"""Propagate monthly source-cell ET reduction fractions through fixed UTrack flows.

Required CSV columns are cellid, month (1-12), and --et-change-column.
A change of 0.01 reduces current ERA5_ET by 1%; negative changes increase ET.
Explicit missing changes (empty, NA, NaN) exclude that source/month with a warning;
absent rows are errors. All configured polygons must have a row for every month.

Each monthly global NetCDF contains signed scenario-minus-current precipitation
changes in mm and fractional prc (-0.01 means a 1% decrease), plus volume change
and receiving-cell area. Relative changes are undefined where ERA5_TP is zero
or missing. Transport is fixed, with no secondary precipitation/ET feedback.
Areas use spherical coordinate-midpoint bounds and radius 6,371,000 m.
UTrack encoded-log sentinel 255 denotes an unavailable link, excluded with a
warning. Other invalid flows and invalid source ET volumes are errors.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
from pathlib import Path
import warnings

import numpy as np
import shapefile
import xarray as xr


ENGINE_PATH = Path(__file__).resolve().with_name("UTrack-Polygon_01_MonthlyP+Eshed.py")
EARTH_RADIUS_M = 6_371_000.0
MISSING_CHANGES = {"", "na", "nan"}


def load_engine():
    spec = importlib.util.spec_from_file_location("utrack_polygon_engine", ENGINE_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load polygon calculation engine: {ENGINE_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario-table", type=Path, required=True)
    parser.add_argument("--et-change-column", required=True, help="Exact CSV column name containing ET reduction fractions.")
    parser.add_argument("--polygon-grid-file", type=Path, required=True)
    parser.add_argument("--utrack-directory", type=Path, required=True)
    parser.add_argument("--era5-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--flow-encoding", choices=("encoded-log", "fraction"), default="encoded-log")
    parser.add_argument("--batch-size", type=int, default=16)
    arguments = parser.parse_args()
    if arguments.batch_size < 1:
        parser.error("--batch-size must be positive.")
    if arguments.et_change_column in {"cellid", "month"}:
        parser.error("--et-change-column must be distinct from cellid and month.")
    return arguments


def read_scenario_table(path: Path, change_column: str) -> dict[tuple[int, int], float | None]:
    if change_column in {"cellid", "month"}:
        raise ValueError("The ET change column must be distinct from cellid and month.")
    changes: dict[tuple[int, int], float | None] = {}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames or []
        required = {"cellid", "month", change_column}
        if not required.issubset(fields):
            raise ValueError(f"{path} is missing required columns: {sorted(required - set(fields))}")
        if len(fields) != len(set(fields)):
            raise ValueError(f"{path} contains duplicate column names.")
        for line, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}, line {line}: malformed CSV row.")
            try:
                cellid = int(row["cellid"])
                month = int(row["month"])
            except ValueError as error:
                raise ValueError(f"{path}, line {line}: cellid and month must be integers.") from error
            if not 1 <= month <= 12:
                raise ValueError(f"{path}, line {line}: month must be between 1 and 12.")
            key = (cellid, month)
            if key in changes:
                raise ValueError(f"{path}, line {line}: duplicate cellid/month {key}.")
            text = row[change_column].strip()
            if text.lower() in MISSING_CHANGES:
                changes[key] = None
                continue
            try:
                fraction = float(text)
            except ValueError as error:
                raise ValueError(f"{path}, line {line}: invalid ET change {text!r}.") from error
            if not np.isfinite(fraction) or fraction > 1:
                raise ValueError(f"{path}, line {line}: ET reduction must be finite and <= 1; got {text!r}.")
            changes[key] = fraction
    if not changes:
        raise ValueError(f"No scenario rows found in {path}.")
    return changes


def read_source_polygons(path: Path) -> tuple[np.ndarray, np.ndarray]:
    if path.suffix.lower() != ".shp" or not path.is_file():
        raise FileNotFoundError(f"Polygon grid shapefile not found: {path}")
    with shapefile.Reader(str(path)) as reader:
        if "cellid" not in [field[0] for field in reader.fields[1:]]:
            raise ValueError(f"{path} must contain a cellid field.")
        ids = []
        centers = []
        for record in reader.iterShapeRecords():
            value = record.record["cellid"]
            if value is None or int(value) != value:
                raise ValueError(f"{path} has a noninteger cellid: {value!r}.")
            ids.append(int(value))
            west, south, east, north = record.shape.bbox
            centers.append(((south + north) / 2, (west + east) / 2))
    if not ids or len(set(ids)) != len(ids):
        raise ValueError(f"{path} must contain nonempty, unique cellid records.")
    return np.asarray(ids), np.asarray(centers)


def validate_coverage(changes: dict[tuple[int, int], float | None], cellids: np.ndarray) -> None:
    expected = {(int(cellid), month) for cellid in cellids for month in range(1, 13)}
    missing = expected - changes.keys()
    extra = changes.keys() - expected
    if missing or extra:
        raise ValueError(
            "Scenario rows must match all configured source cells and months 1-12. "
            f"Missing ({len(missing)}): {sorted(missing)[:10]}; "
            f"unexpected ({len(extra)}): {sorted(extra)[:10]}."
        )


def source_indices(centers: np.ndarray, flow: xr.DataArray, engine) -> tuple[np.ndarray, np.ndarray]:
    lat_indices = engine.nearest_latitude_indices(flow.sourcelat.values, centers[:, 0])
    lon_indices = engine.nearest_longitude_indices(flow.sourcelon.values, centers[:, 1])
    lat_error = np.abs(flow.sourcelat.values[lat_indices] - centers[:, 0])
    lon_error = np.abs((flow.sourcelon.values[lon_indices] - centers[:, 1] + 180) % 360 - 180)
    if not np.isfinite(centers).all() or np.any(lat_error > 1e-6) or np.any(lon_error > 1e-6):
        raise ValueError("Polygon centers must match UTrack source-cell centers (within 1e-6 degrees).")
    if len(np.unique(np.column_stack((lat_indices, lon_indices)), axis=0)) != len(centers):
        raise ValueError("Distinct polygon cellids map to the same UTrack source cell.")
    return lat_indices, lon_indices


def regular_spacing(values: np.ndarray, name: str) -> float:
    if len(values) < 2 or not np.isfinite(values).all():
        raise ValueError(f"{name} must contain at least two finite coordinates.")
    differences = np.diff(values)
    step = float(differences[0])
    if step == 0 or not np.allclose(differences, step, rtol=0, atol=1e-6):
        raise ValueError(f"{name} must be a strictly monotonic, regularly spaced grid.")
    return step


def receiving_cell_areas(flow: xr.DataArray) -> xr.DataArray:
    latitudes = np.asarray(flow.targetlat.values, dtype=float)
    longitudes = np.asarray(flow.targetlon.values, dtype=float)
    lat_step = abs(regular_spacing(latitudes, "targetlat"))
    lon_step = abs(regular_spacing(longitudes, "targetlon"))
    if np.any(np.abs(latitudes) > 90) or lon_step * len(longitudes) > 360 + 1e-6:
        raise ValueError("Receiving coordinates exceed geographic grid bounds.")
    lower = np.clip(latitudes - lat_step / 2, -90, 90)
    upper = np.clip(latitudes + lat_step / 2, -90, 90)
    # Global ERA5 coordinates include 90 but stop at -89.5; close polar caps.
    if latitudes.min() <= -90 + lat_step:
        lower[np.argmin(latitudes)] = -90
    if latitudes.max() >= 90 - lat_step:
        upper[np.argmax(latitudes)] = 90
    bands = EARTH_RADIUS_M**2 * np.deg2rad(lon_step) * (
        np.sin(np.deg2rad(upper)) - np.sin(np.deg2rad(lower))
    )
    return xr.DataArray(
        np.broadcast_to(bands[:, None], (len(latitudes), len(longitudes))).copy(),
        dims=("targetlat", "targetlon"),
        coords={"targetlat": flow.targetlat, "targetlon": flow.targetlon},
        name="receiving_cell_area_m2",
        attrs={"units": "m2", "long_name": "Spherical receiving-cell area"},
    )


def validate_volumes(era5: xr.Dataset) -> None:
    for name in ("ERA5_ET", "ERA5_TP"):
        if name not in era5 or set(era5[name].dims) != {"lat", "lon"}:
            raise ValueError(f"Expected {name} with lat and lon dimensions.")
        units = era5[name].attrs.get("units", "")
        if units not in {"m3", "m^3", "m**3", "m3 month-1", "m^3 month-1"}:
            raise ValueError(f"{name} must contain monthly volumes in m3; got units {units!r}.")


def propagate_month(
    flow: xr.DataArray,
    era5: xr.Dataset,
    lat_indices: np.ndarray,
    lon_indices: np.ndarray,
    reductions: np.ndarray,
    batch_size: int,
    encoding: str,
    engine,
) -> xr.Dataset:
    engine.validate_inputs(flow, era5)
    validate_volumes(era5)
    if batch_size < 1 or encoding not in {"encoded-log", "fraction"}:
        raise ValueError("Expected a positive batch size and supported flow encoding.")
    if not (len(lat_indices) == len(lon_indices) == len(reductions)) or not len(reductions):
        raise ValueError("At least one source is required, with matching indices and reductions.")
    if not np.isfinite(reductions).all() or np.any(reductions > 1):
        raise ValueError("ET reductions must be finite and <= 1.")
    area = receiving_cell_areas(flow)
    total = xr.zeros_like(area)
    unavailable_count = 0
    for start in range(0, len(reductions), batch_size):
        stop = min(start + batch_size, len(reductions))
        lat = xr.DataArray(lat_indices[start:stop], dims="polygon_cell")
        lon = xr.DataArray(lon_indices[start:stop], dims="polygon_cell")
        evaporation = era5["ERA5_ET"].isel(lat=lat, lon=lon).load()
        if not np.isfinite(evaporation.values).all() or np.any(evaporation.values < 0):
            raise ValueError("Selected source ERA5_ET volumes must be finite and nonnegative.")
        selected = flow.isel(sourcelat=lat, sourcelon=lon).load()
        if encoding == "encoded-log":
            unavailable = selected == 255
            unavailable_count += int(unavailable.sum())
            selected = selected.where(~unavailable, 0)
            invalid = (~np.isfinite(selected)) | (selected < 0) | (selected > 254)
        else:
            invalid = (~np.isfinite(selected)) | (selected < 0) | (selected > 1)
        if bool(invalid.any()):
            raise ValueError("Selected UTrack flows contain invalid values.")
        delta_et = -evaporation * xr.DataArray(reductions[start:stop], dims="polygon_cell")
        contribution = (engine.decode_flow(selected, encoding) * delta_et).sum("polygon_cell", skipna=False)
        total = total + contribution
    if unavailable_count:
        warnings.warn(
            f"Excluded {unavailable_count} unavailable UTrack source-to-target links (encoded value 255).",
            UserWarning,
            stacklevel=2,
        )
    precipitation = era5["ERA5_TP"].rename({"lat": "targetlat", "lon": "targetlon"}).load()
    if np.isinf(precipitation.values).any() or bool((precipitation < 0).any()):
        raise ValueError("ERA5_TP must be nonnegative and cannot contain infinite values.")
    undefined = precipitation.isnull() | (precipitation == 0)
    undefined_count = int(undefined.sum())
    if undefined_count:
        warnings.warn(
            f"Relative precipitation change is undefined in {undefined_count} zero/missing ERA5_TP cells.",
            UserWarning,
            stacklevel=2,
        )
    if not np.isfinite(total.values).all():
        raise ValueError("Propagated precipitation volume change contains nonfinite values.")
    result = xr.Dataset(
        {
            "precipitation_change_m3": total,
            "precipitation_change_mm": total / area * 1000,
            "precipitation_change_prc": total / precipitation.where(~undefined),
            "receiving_cell_area_m2": area,
        },
        attrs={
            "sign_convention": "scenario minus current; positive input fractions reduce ET",
            "transport_assumption": "Fixed UTrack transport; no secondary precipitation/ET feedback",
            "et_change_formula": "delta_ET_m3 = -input_fraction * ERA5_ET",
            "precipitation_change_formula": "delta_P_m3 = sum_source(flow_fraction * delta_ET_m3)",
            "mm_formula": "1000 * delta_P_m3 / receiving_cell_area_m2",
            "prc_formula": "delta_P_m3 / ERA5_TP; fractional, not percentage units",
            "area_method": "Spherical coordinate-midpoint bounds, closed polar caps",
            "earth_radius_m": EARTH_RADIUS_M,
            "unavailable_flow_link_count": unavailable_count,
            "undefined_relative_change_cell_count": undefined_count,
        },
    )
    for suffix, units in (("m3", "m3 month-1"), ("mm", "mm month-1"), ("prc", "1")):
        result[f"precipitation_change_{suffix}"].attrs = {
            "units": units,
            "long_name": f"Monthly scenario-minus-current precipitation change ({suffix})",
        }
    return result


def main() -> None:
    arguments = parse_arguments()
    engine = load_engine()
    changes = read_scenario_table(arguments.scenario_table, arguments.et_change_column)
    cellids, centers = read_source_polygons(arguments.polygon_grid_file)
    validate_coverage(changes, cellids)
    # Validate coverage and monthly availability before writing any output.
    for month in range(1, 13):
        engine.monthly_paths(arguments.utrack_directory, arguments.era5_directory, month)
        if all(changes[(int(cellid), month)] is None for cellid in cellids):
            raise ValueError(f"Month {month:02d}: all source ET changes are missing.")
    arguments.output_directory.mkdir(parents=True, exist_ok=True)
    for month in range(1, 13):
        included = np.array([changes[(int(cellid), month)] is not None for cellid in cellids])
        excluded = cellids[~included]
        if len(excluded):
            warnings.warn(
                f"Month {month:02d}: excluding {len(excluded)} source cells with missing ET changes: {excluded.tolist()}",
                UserWarning,
            )
        utrack_path, era5_path = engine.monthly_paths(arguments.utrack_directory, arguments.era5_directory, month)
        with xr.open_dataset(utrack_path) as utrack, xr.open_dataset(era5_path) as era5:
            flow = utrack["moisture_flow"]
            lat, lon = source_indices(centers, flow, engine)
            reductions = np.array([changes[(int(cellid), month)] for cellid in cellids[included]], dtype=float)
            result = propagate_month(
                flow, era5, lat[included], lon[included], reductions,
                arguments.batch_size, arguments.flow_encoding, engine,
            )
            result.attrs.update(
                title="Monthly precipitation change from source-cell ET scenario",
                month=month,
                scenario_table=str(arguments.scenario_table),
                et_change_column=arguments.et_change_column,
                polygon_grid_file=str(arguments.polygon_grid_file),
                utrack_file=str(utrack_path),
                era5_file=str(era5_path),
                flow_encoding=arguments.flow_encoding,
                included_source_cell_count=int(included.sum()),
                excluded_source_cell_count=len(excluded),
                excluded_source_cellids=",".join(str(int(cellid)) for cellid in excluded),
            )
        output_path = arguments.output_directory / f"utrack_polygon_ETchanges_month_{month:02d}.nc"
        result.to_netcdf(output_path)
        print(f"Month {month:02d}: wrote {output_path} ({int(included.sum())} sources, {len(excluded)} excluded).")


if __name__ == "__main__":
    main()
