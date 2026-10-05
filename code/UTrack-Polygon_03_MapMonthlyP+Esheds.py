"""Map monthly polygon sheds, including evaporationshed precipitation fractions."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import numpy as np
import shapefile
import xarray as xr


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
DEFAULT_STATE_BOUNDARY_FILE = SCRIPT_DIRECTORY.parent / "data" / "boundaries" / "CONUS-States_TIGRIS.shp"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument(
        "--evaporationshed-input-directory",
        type=Path,
        required=True,
        help="Directory containing utrack_polygon_evaporationshed_month_*.nc.",
    )
    parser.add_argument(
        "--era5-directory",
        type=Path,
        required=True,
        help="Directory containing ERA5_01_0.5_volumes_corrected.nc through _12_0.5_volumes_corrected.nc.",
    )
    parser.add_argument(
        "--evaporationshed-fraction-output-directory",
        type=Path,
        required=True,
        help="Directory for maps of evaporationshed contributions as receiving-cell precipitation percentages.",
    )
    parser.add_argument(
        "--evaporationshed-map-output-directory",
        type=Path,
        required=True,
        help="Directory for volumetric evaporationshed maps.",
    )
    parser.add_argument(
        "--state-boundary-file",
        type=Path,
        default=DEFAULT_STATE_BOUNDARY_FILE,
        help=f"State-boundary shapefile to overlay (default: {DEFAULT_STATE_BOUNDARY_FILE}).",
    )
    return parser.parse_args()


def load_state_boundaries(boundary_file: Path) -> tuple[list[np.ndarray], tuple[float, float, float, float]]:
    if boundary_file.suffix.lower() != ".shp" or not boundary_file.is_file():
        raise FileNotFoundError(f"State-boundary shapefile not found: {boundary_file}")
    reader = shapefile.Reader(str(boundary_file))
    lines = []
    for shape in reader.shapes():
        part_starts = list(shape.parts) + [len(shape.points)]
        for start, stop in zip(part_starts[:-1], part_starts[1:]):
            points = np.asarray(shape.points[start:stop])
            if points.size:
                points[:, 0] = (points[:, 0] + 180) % 360 - 180
                lines.append(points)
    if not lines:
        raise ValueError(f"No state boundaries found in {boundary_file}")
    west, south, east, north = reader.bbox
    return lines, ((west + 180) % 360 - 180, south, (east + 180) % 360 - 180, north)


def draw_state_boundaries(axis, state_boundaries: list[np.ndarray]) -> None:
    for boundary in state_boundaries:
        axis.plot(boundary[:, 0], boundary[:, 1], color="black", linewidth=0.4, zorder=2)


def apply_conus_extent(axis, bounds: tuple[float, float, float, float]) -> None:
    west, south, east, north = bounds
    axis.set(xlim=(west, east), ylim=(south, north), aspect="equal")


def main() -> None:
    arguments = parse_arguments()
    arguments.output_directory.mkdir(parents=True, exist_ok=True)
    arguments.evaporationshed_fraction_output_directory.mkdir(parents=True, exist_ok=True)
    arguments.evaporationshed_map_output_directory.mkdir(parents=True, exist_ok=True)
    state_boundaries, state_bounds = load_state_boundaries(arguments.state_boundary_file)
    monthly_fields = []
    for month in range(1, 13):
        input_path = arguments.input_directory / f"utrack_polygon_precipitationshed_month_{month:02d}.nc"
        if not input_path.is_file():
            raise FileNotFoundError(f"Expected monthly precipitationshed output not found: {input_path}")
        with xr.open_dataset(input_path) as dataset:
            if "evaporation_contribution_m3" not in dataset:
                raise KeyError(f"{input_path} does not contain evaporation_contribution_m3.")
            monthly_fields.append((month, dataset["evaporation_contribution_m3"].load()))

    positive = [field.values[np.isfinite(field.values) & (field.values > 0)] for _, field in monthly_fields]
    positive = np.concatenate([values for values in positive if values.size])
    if not positive.size:
        raise ValueError("No positive precipitationshed contributions are available to map.")
    limits = (float(positive.min()), float(positive.max()))

    for month, field in monthly_fields:
        longitudes = (field.sourcelon.values + 180) % 360 - 180
        longitude_order = np.argsort(longitudes)
        values = field.isel(sourcelon=longitude_order).values
        values = np.ma.masked_where(~np.isfinite(values) | (values <= 0), values)
        figure, axis = plt.subplots(figsize=(12, 6), constrained_layout=True)
        image = axis.pcolormesh(
            longitudes[longitude_order],
            field.sourcelat.values,
            values,
            cmap="viridis",
            norm=LogNorm(*limits),
            shading="auto",
        )
        axis.set(
            title=f"UTrack polygon precipitationshed: month {month:02d}",
            xlabel="Longitude",
            ylabel="Latitude",
            xlim=(-180, 180),
            ylim=(-90, 90),
            aspect="equal",
        )
        draw_state_boundaries(axis, state_boundaries)
        figure.colorbar(image, ax=axis, label="Evaporation contribution (m3 month-1)")
        output_path = arguments.output_directory / f"utrack_polygon_precipitationshed_month_{month:02d}.png"
        figure.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(figure)
        print(f"Month {month:02d}: wrote {output_path}")

    monthly_evaporationsheds = []
    monthly_fractions = []
    for month in range(1, 13):
        evaporationshed_path = (
            arguments.evaporationshed_input_directory
            / f"utrack_polygon_evaporationshed_month_{month:02d}.nc"
        )
        era5_path = arguments.era5_directory / f"ERA5_{month:02d}_0.5_volumes_corrected.nc"
        missing = [path for path in (evaporationshed_path, era5_path) if not path.is_file()]
        if missing:
            raise FileNotFoundError("Missing input file(s):\n" + "\n".join(f"  - {path}" for path in missing))
        with xr.open_dataset(evaporationshed_path) as evaporationshed, xr.open_dataset(era5_path) as era5:
            if "precipitation_contribution_m3" not in evaporationshed:
                raise KeyError(f"{evaporationshed_path} does not contain precipitation_contribution_m3.")
            if "ERA5_TP" not in era5:
                raise KeyError(f"{era5_path} does not contain ERA5_TP.")
            contribution = evaporationshed["precipitation_contribution_m3"].load()
            monthly_evaporationsheds.append((month, contribution))
            total_precipitation = era5["ERA5_TP"].rename({"lat": "targetlat", "lon": "targetlon"}).sel(
                targetlat=contribution.targetlat,
                targetlon=contribution.targetlon,
            ).load()
            monthly_fractions.append(
                (month, xr.where(total_precipitation > 0, contribution / total_precipitation * 100, np.nan))
            )

    positive_evaporationshed_values = [
        contribution.values[np.isfinite(contribution.values) & (contribution.values > 0)]
        for _, contribution in monthly_evaporationsheds
    ]
    positive_evaporationshed_values = [values for values in positive_evaporationshed_values if values.size]
    if not positive_evaporationshed_values:
        raise ValueError("No positive evaporationshed contributions are available to map.")
    evaporationshed_limits = (
        float(np.concatenate(positive_evaporationshed_values).min()),
        float(np.concatenate(positive_evaporationshed_values).max()),
    )

    for month, contribution in monthly_evaporationsheds:
        longitudes = (contribution.targetlon.values + 180) % 360 - 180
        longitude_order = np.argsort(longitudes)
        values = contribution.isel(targetlon=longitude_order).values
        values = np.ma.masked_where(~np.isfinite(values) | (values <= 0), values)
        figure, axis = plt.subplots(figsize=(12, 6), constrained_layout=True)
        image = axis.pcolormesh(
            longitudes[longitude_order],
            contribution.targetlat.values,
            values,
            cmap="viridis",
            norm=LogNorm(*evaporationshed_limits),
            shading="auto",
        )
        axis.set(
            title=f"UTrack polygon evaporationshed: month {month:02d}",
            xlabel="Longitude",
            ylabel="Latitude",
        )
        draw_state_boundaries(axis, state_boundaries)
        apply_conus_extent(axis, state_bounds)
        figure.colorbar(image, ax=axis, label="Precipitation contribution (m3 month-1)")
        output_path = arguments.evaporationshed_map_output_directory / f"utrack_polygon_evaporationshed_month_{month:02d}.png"
        figure.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(figure)
        print(f"Month {month:02d}: wrote {output_path}")

    if not any(np.any(np.isfinite(fraction.values) & (fraction.values > 0)) for _, fraction in monthly_fractions):
        raise ValueError("No positive evaporationshed precipitation fractions are available to map.")

    for month, fraction in monthly_fractions:
        longitudes = (fraction.targetlon.values + 180) % 360 - 180
        longitude_order = np.argsort(longitudes)
        values = fraction.isel(targetlon=longitude_order).values
        values = np.ma.masked_where(~np.isfinite(values) | (values <= 0), values)
        figure, axis = plt.subplots(figsize=(12, 6), constrained_layout=True)
        image = axis.pcolormesh(
            longitudes[longitude_order],
            fraction.targetlat.values,
            values,
            cmap="viridis",
            norm=Normalize(vmin=0, vmax=100),
            shading="auto",
        )
        axis.set(
            title=f"UTrack polygon evaporationshed fraction of receiving precipitation: month {month:02d}",
            xlabel="Longitude",
            ylabel="Latitude",
        )
        draw_state_boundaries(axis, state_boundaries)
        apply_conus_extent(axis, state_bounds)
        figure.colorbar(image, ax=axis, label="Contribution to receiving-cell precipitation (%)")
        output_path = (
            arguments.evaporationshed_fraction_output_directory
            / f"utrack_polygon_evaporationshed_precipitation_fraction_month_{month:02d}.png"
        )
        figure.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(figure)
        print(f"Month {month:02d}: wrote {output_path}")


if __name__ == "__main__":
    main()
