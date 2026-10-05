"""Calculate monthly UTrack polygon precipitationsheds and evaporationsheds.

The polygon is represented by grid-cell polygons, such as
data\\ERA5_grid\\ERA5_grid_HPA.shp. Every grid cell in the shapefile is used.
For each month, the script writes an upwind precipitationshed and a downwind
evaporationshed without creating a global four-dimensional intermediate array.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import shapefile
import xarray as xr


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polygon-grid-file", type=Path, required=True)
    parser.add_argument("--utrack-directory", type=Path, required=True)
    parser.add_argument("--era5-directory", type=Path, required=True)
    parser.add_argument("--precipitationshed-output-directory", type=Path, required=True)
    parser.add_argument("--evaporationshed-output-directory", type=Path, required=True)
    parser.add_argument("--flow-encoding", choices=("encoded-log", "fraction"), default="encoded-log")
    parser.add_argument("--batch-size", type=int, default=16)
    arguments = parser.parse_args()
    if arguments.batch_size < 1:
        parser.error("--batch-size must be positive.")
    return arguments


def monthly_paths(utrack_directory: Path, era5_directory: Path, month: int) -> tuple[Path, Path]:
    utrack_path = utrack_directory / f"utrack_climatology_0.5_{month:02d}.nc"
    era5_path = era5_directory / f"ERA5_{month:02d}_0.5_volumes_corrected.nc"
    missing = [path for path in (utrack_path, era5_path) if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing input file(s):\n" + "\n".join(f"  - {path}" for path in missing))
    return utrack_path, era5_path


def nearest_latitude_indices(grid: np.ndarray, values: np.ndarray) -> np.ndarray:
    return np.abs(grid[:, None] - values).argmin(axis=0)


def nearest_longitude_indices(grid: np.ndarray, values: np.ndarray) -> np.ndarray:
    grid_wgs84 = (grid + 180) % 360 - 180
    values_wgs84 = (values + 180) % 360 - 180
    distance = np.abs((grid_wgs84[:, None] - values_wgs84 + 180) % 360 - 180)
    return distance.argmin(axis=0)


def polygon_cell_indices(polygon_grid_file: Path, moisture_flow: xr.DataArray) -> tuple[np.ndarray, np.ndarray]:
    if polygon_grid_file.suffix.lower() != ".shp" or not polygon_grid_file.is_file():
        raise FileNotFoundError(f"Polygon grid shapefile not found: {polygon_grid_file}")
    shapes = shapefile.Reader(str(polygon_grid_file)).shapes()
    if not shapes:
        raise ValueError(f"No grid-cell polygons found in {polygon_grid_file}.")
    centers = np.array([((shape.bbox[1] + shape.bbox[3]) / 2, (shape.bbox[0] + shape.bbox[2]) / 2) for shape in shapes])
    latitude_indices = nearest_latitude_indices(moisture_flow.sourcelat.values, centers[:, 0])
    longitude_indices = nearest_longitude_indices(moisture_flow.sourcelon.values, centers[:, 1])
    unique = np.unique(np.column_stack((latitude_indices, longitude_indices)), axis=0)
    if unique.size == 0:
        raise ValueError("The polygon grid did not select any UTrack cells.")
    return unique[:, 0], unique[:, 1]


def validate_inputs(moisture_flow: xr.DataArray, era5_volumes: xr.Dataset) -> None:
    expected_dims = {"sourcelat", "sourcelon", "targetlat", "targetlon"}
    if not expected_dims.issubset(moisture_flow.dims) or "ERA5_ET" not in era5_volumes:
        raise ValueError("Expected UTrack moisture_flow and ERA5_ET input variables.")
    for latitudes, longitudes in ((moisture_flow.sourcelat, moisture_flow.sourcelon), (moisture_flow.targetlat, moisture_flow.targetlon)):
        if not np.array_equal(latitudes, era5_volumes.lat) or not np.array_equal(longitudes, era5_volumes.lon):
            raise ValueError("UTrack and ERA5 grids do not match.")


def decode_flow(values: xr.DataArray, encoding: str) -> xr.DataArray:
    return values.where(values > 0, 0) if encoding == "fraction" else xr.where(values == 0, 0, np.exp(-0.1 * values))


def aggregate_precipitationshed(flow: xr.DataArray, evaporation: xr.DataArray, lat_indices: np.ndarray, lon_indices: np.ndarray, batch_size: int, encoding: str) -> xr.DataArray:
    total = None
    for start in range(0, len(lat_indices), batch_size):
        indices = slice(start, min(start + batch_size, len(lat_indices)))
        selected = flow.isel(
            targetlat=xr.DataArray(lat_indices[indices], dims="polygon_cell"),
            targetlon=xr.DataArray(lon_indices[indices], dims="polygon_cell"),
        ).load()
        contribution = (decode_flow(selected, encoding) * evaporation.rename({"lat": "sourcelat", "lon": "sourcelon"})).sum("polygon_cell", skipna=True)
        total = contribution if total is None else total + contribution
    return total.rename("evaporation_contribution_m3")


def aggregate_evaporationshed(flow: xr.DataArray, evaporation: xr.DataArray, lat_indices: np.ndarray, lon_indices: np.ndarray, batch_size: int, encoding: str) -> xr.DataArray:
    total = None
    for start in range(0, len(lat_indices), batch_size):
        indices = slice(start, min(start + batch_size, len(lat_indices)))
        source_cell = xr.DataArray(np.arange(start, min(start + batch_size, len(lat_indices))), dims="polygon_cell")
        selected = flow.isel(
            sourcelat=xr.DataArray(lat_indices[indices], dims="polygon_cell"),
            sourcelon=xr.DataArray(lon_indices[indices], dims="polygon_cell"),
        ).load()
        source_evaporation = evaporation.isel(
            lat=xr.DataArray(lat_indices[indices], dims="polygon_cell"),
            lon=xr.DataArray(lon_indices[indices], dims="polygon_cell"),
        ).assign_coords(polygon_cell=source_cell)
        contribution = (decode_flow(selected, encoding) * source_evaporation).sum("polygon_cell", skipna=True)
        total = contribution if total is None else total + contribution
    return total.rename("precipitation_contribution_m3")


def main() -> None:
    arguments = parse_arguments()
    arguments.precipitationshed_output_directory.mkdir(parents=True, exist_ok=True)
    arguments.evaporationshed_output_directory.mkdir(parents=True, exist_ok=True)
    selected_count = None
    for month in range(1, 13):
        utrack_path, era5_path = monthly_paths(arguments.utrack_directory, arguments.era5_directory, month)
        with xr.open_dataset(utrack_path) as utrack, xr.open_dataset(era5_path) as era5:
            flow = utrack["moisture_flow"]
            validate_inputs(flow, era5)
            lat_indices, lon_indices = polygon_cell_indices(arguments.polygon_grid_file, flow)
            selected_count = len(lat_indices)
            precipitationshed = aggregate_precipitationshed(flow, era5["ERA5_ET"], lat_indices, lon_indices, arguments.batch_size, arguments.flow_encoding)
            evaporationshed = aggregate_evaporationshed(flow, era5["ERA5_ET"], lat_indices, lon_indices, arguments.batch_size, arguments.flow_encoding)
            attributes = {"polygon_grid_file": str(arguments.polygon_grid_file), "selected_grid_cell_count": selected_count, "flow_encoding": arguments.flow_encoding}
            precipitationshed = precipitationshed.to_dataset()
            evaporationshed = evaporationshed.to_dataset()
            precipitationshed.attrs = {"title": "Monthly UTrack polygon precipitationshed", **attributes}
            evaporationshed.attrs = {"title": "Monthly UTrack polygon evaporationshed", **attributes}
        precipitationshed.to_netcdf(arguments.precipitationshed_output_directory / f"utrack_polygon_precipitationshed_month_{month:02d}.nc")
        evaporationshed.to_netcdf(arguments.evaporationshed_output_directory / f"utrack_polygon_evaporationshed_month_{month:02d}.nc")
        print(f"Month {month:02d}: wrote polygon sheds for {selected_count} grid cells.")


if __name__ == "__main__":
    main()
