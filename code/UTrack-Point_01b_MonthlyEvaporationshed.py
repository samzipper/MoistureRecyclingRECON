"""Calculate monthly UTrack evaporationsheds for one source coordinate.

For each monthly UTrack file, this script selects the UTrack source cell nearest
to the requested latitude/longitude. It converts that source cell's
source-to-target flows to fractions, weights them by the source cell's monthly
ERA5_ET volume, and writes one NetCDF file containing the downwind
precipitation contribution at every target cell.

Example:
    python code\\UTrack_05_MonthlyEvaporationshed.py ^
      --latitude 37.78228045733923 --longitude -100.90621369000884 ^
      --utrack-directory "C:\\data\\UTrack" ^
      --era5-directory "C:\\data\\RECON" ^
      --output-directory "figures+tables\\UTrack_monthly_evaporationshed"
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr


SOURCE_DIMS = ("sourcelat", "sourcelon")
TARGET_DIMS = ("targetlat", "targetlon")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--latitude", type=float, required=True, help="Source latitude in degrees north.")
    parser.add_argument("--longitude", type=float, required=True, help="Source longitude in degrees east or west.")
    parser.add_argument("--utrack-directory", type=Path, required=True)
    parser.add_argument("--era5-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--flow-encoding", choices=("encoded-log", "fraction"), default="encoded-log")
    return parser.parse_args()


def monthly_paths(utrack_directory: Path, era5_directory: Path, month: int) -> tuple[Path, Path]:
    utrack_path = utrack_directory / f"utrack_climatology_0.5_{month:02d}.nc"
    era5_path = era5_directory / f"ERA5_{month:02d}_0.5_volumes_corrected.nc"
    missing_paths = [path for path in (utrack_path, era5_path) if not path.is_file()]
    if missing_paths:
        missing_list = "\n".join(f"  - {path}" for path in missing_paths)
        raise FileNotFoundError(f"Month {month:02d} is missing input file(s):\n{missing_list}")
    return utrack_path, era5_path


def nearest_latitude_index(values: np.ndarray, requested_latitude: float) -> int:
    if not -90 <= requested_latitude <= 90:
        raise ValueError(f"Latitude must be between -90 and 90 degrees; got {requested_latitude}.")
    return int(np.abs(values - requested_latitude).argmin())


def nearest_longitude_index(values: np.ndarray, requested_longitude: float) -> int:
    normalized_longitude = (requested_longitude + 180) % 360 - 180
    values_wgs84 = (values + 180) % 360 - 180
    distances = np.abs((values_wgs84 - normalized_longitude + 180) % 360 - 180)
    return int(distances.argmin())


def validate_inputs(moisture_flow: xr.DataArray, era5_volumes: xr.Dataset, utrack_path: Path, era5_path: Path) -> None:
    if not set(SOURCE_DIMS + TARGET_DIMS).issubset(moisture_flow.dims):
        raise ValueError(f"{utrack_path.name} has unexpected moisture_flow dimensions: {moisture_flow.dims}.")
    if "ERA5_ET" not in era5_volumes or not {"lat", "lon"}.issubset(era5_volumes["ERA5_ET"].dims):
        raise KeyError(f"{era5_path.name} must contain ERA5_ET with lat and lon dimensions.")
    for latitudes, longitudes, grid_name in (
        (moisture_flow.sourcelat, moisture_flow.sourcelon, "source"),
        (moisture_flow.targetlat, moisture_flow.targetlon, "target"),
    ):
        if not np.array_equal(latitudes.values, era5_volumes.lat.values) or not np.array_equal(
            longitudes.values, era5_volumes.lon.values
        ):
            raise ValueError(f"The ERA5 grid in {era5_path.name} does not match the UTrack {grid_name} grid.")


def decode_flow(encoded_flow: xr.DataArray, flow_encoding: str) -> xr.DataArray:
    if flow_encoding == "fraction":
        return encoded_flow.where(encoded_flow > 0, 0)
    return xr.where(encoded_flow == 0, 0, np.exp(-0.1 * encoded_flow))


def calculate_monthly_evaporationshed(
    moisture_flow: xr.DataArray,
    evaporation_volume: xr.DataArray,
    latitude: float,
    longitude: float,
    flow_encoding: str,
) -> xr.Dataset:
    source_latitude_index = nearest_latitude_index(moisture_flow.sourcelat.values, latitude)
    source_longitude_index = nearest_longitude_index(moisture_flow.sourcelon.values, longitude)
    encoded_source_flow = moisture_flow.isel(
        sourcelat=source_latitude_index, sourcelon=source_longitude_index
    ).load()
    flow_fraction = decode_flow(encoded_source_flow, flow_encoding)
    source_evaporation = evaporation_volume.isel(lat=source_latitude_index, lon=source_longitude_index).load()
    contribution = (flow_fraction * source_evaporation).rename("precipitation_contribution_m3")
    source_latitude = float(moisture_flow.sourcelat.isel(sourcelat=source_latitude_index))
    source_longitude = float(moisture_flow.sourcelon.isel(sourcelon=source_longitude_index))

    result = xr.Dataset(
        {
            "precipitation_contribution_m3": contribution,
            "moisture_flow_fraction": flow_fraction.rename("moisture_flow_fraction"),
        },
        attrs={
            "title": "Monthly UTrack evaporationshed",
            "requested_source_latitude_degrees_north": latitude,
            "requested_source_longitude_degrees_east": (longitude + 180) % 360 - 180,
            "selected_source_latitude_degrees_north": source_latitude,
            "selected_source_longitude_degrees_east": (source_longitude + 180) % 360 - 180,
            "source_evaporation_m3_month-1": float(source_evaporation),
            "flow_encoding": flow_encoding,
            "description": "Downwind precipitation contributions from evaporation at the selected source cell.",
        },
    )
    result["precipitation_contribution_m3"].attrs.update(
        long_name="Monthly precipitation contribution from selected source evaporation",
        units="m3 month-1",
    )
    result["moisture_flow_fraction"].attrs.update(
        long_name="UTrack source-to-target moisture flow fraction",
        units="1",
    )
    return result


def main() -> None:
    arguments = parse_arguments()
    arguments.output_directory.mkdir(parents=True, exist_ok=True)
    for month in range(1, 13):
        utrack_path, era5_path = monthly_paths(arguments.utrack_directory, arguments.era5_directory, month)
        with xr.open_dataset(utrack_path) as utrack_dataset, xr.open_dataset(era5_path) as era5_volumes:
            validate_inputs(utrack_dataset["moisture_flow"], era5_volumes, utrack_path, era5_path)
            evaporationshed = calculate_monthly_evaporationshed(
                utrack_dataset["moisture_flow"],
                era5_volumes["ERA5_ET"],
                arguments.latitude,
                arguments.longitude,
                arguments.flow_encoding,
            )
        output_path = arguments.output_directory / f"utrack_evaporationshed_month_{month:02d}.nc"
        evaporationshed.to_netcdf(output_path)
        total = float(evaporationshed["precipitation_contribution_m3"].sum(skipna=True))
        print(f"Month {month:02d}: wrote {output_path} ({total:.3e} m3 month-1)")


if __name__ == "__main__":
    main()
