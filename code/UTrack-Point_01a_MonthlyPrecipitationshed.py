"""Calculate monthly UTrack precipitationsheds for one target coordinate.

For each monthly UTrack file, this script selects the UTrack target cell nearest
to the requested latitude/longitude.  It converts that target cell's
source-to-target flows to fractions, weights them by monthly ERA5_ET source
volumes, and writes one NetCDF file containing the upwind evaporation
contribution of every source cell. This is a precipitationshed: the source
evaporation that contributes to precipitation at the selected target.

The UTrack files used by sifa4152/utrack_dataset_processing store
``moisture_flow`` as -10 times the natural log of the flow fraction, with zero
representing no flow.  Use ``--flow-encoding fraction`` only for files whose
``moisture_flow`` is already a flow fraction.

Example:
    python code\\UTrack_03_MonthlyPrecipitationshed.py ^
      --latitude 39.0 --longitude -98.0 ^
      --utrack-directory "C:\\data\\UTrack" ^
      --era5-directory "C:\\data\\RECON" ^
      --output-directory "figures+tables\\UTrack_monthly_precipitationshed"
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
    parser.add_argument("--latitude", type=float, required=True, help="Target latitude in degrees north.")
    parser.add_argument(
        "--longitude",
        type=float,
        required=True,
        help="Target longitude in degrees east or west; values are normalized to [-180, 180).",
    )
    parser.add_argument(
        "--utrack-directory",
        type=Path,
        required=True,
        help="Directory containing utrack_climatology_0.5_01.nc through _12.nc.",
    )
    parser.add_argument(
        "--era5-directory",
        type=Path,
        required=True,
        help="Directory containing ERA5_01_0.5_volumes_corrected.nc through _12_0.5_volumes_corrected.nc.",
    )
    parser.add_argument("--output-directory", type=Path, required=True, help="Directory for monthly NetCDF outputs.")
    parser.add_argument(
        "--flow-encoding",
        choices=("encoded-log", "fraction"),
        default="encoded-log",
        help="Storage convention of moisture_flow (default: encoded-log).",
    )
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
    circular_distance = np.abs((values_wgs84 - normalized_longitude + 180) % 360 - 180)
    return int(circular_distance.argmin())


def validate_inputs(moisture_flow: xr.DataArray, era5_volumes: xr.Dataset, utrack_path: Path, era5_path: Path) -> None:
    expected_dims = set(SOURCE_DIMS + TARGET_DIMS)
    if not expected_dims.issubset(moisture_flow.dims):
        raise ValueError(f"{utrack_path.name} has unexpected moisture_flow dimensions: {moisture_flow.dims}.")
    if "ERA5_ET" not in era5_volumes:
        raise KeyError(f"{era5_path.name} must contain an ERA5_ET variable.")
    if not {"lat", "lon"}.issubset(era5_volumes["ERA5_ET"].dims):
        raise ValueError(f"{era5_path.name} ERA5_ET must have lat and lon dimensions.")
    for flow_latitudes, flow_longitudes, grid_name in (
        (moisture_flow.sourcelat, moisture_flow.sourcelon, "source"),
        (moisture_flow.targetlat, moisture_flow.targetlon, "target"),
    ):
        if not np.array_equal(flow_latitudes.values, era5_volumes.lat.values) or not np.array_equal(
            flow_longitudes.values, era5_volumes.lon.values
        ):
            raise ValueError(f"The ERA5 grid in {era5_path.name} does not match the UTrack {grid_name} grid.")


def decode_flow(encoded_flow: xr.DataArray, flow_encoding: str) -> xr.DataArray:
    if flow_encoding == "fraction":
        return encoded_flow.where(encoded_flow > 0, 0)
    return xr.where(encoded_flow == 0, 0, np.exp(-0.1 * encoded_flow))


def calculate_monthly_precipitationshed(
    moisture_flow: xr.DataArray,
    evaporation_volume: xr.DataArray,
    latitude: float,
    longitude: float,
    flow_encoding: str,
) -> xr.Dataset:
    target_latitude_index = nearest_latitude_index(moisture_flow.targetlat.values, latitude)
    target_longitude_index = nearest_longitude_index(moisture_flow.targetlon.values, longitude)
    encoded_target_flow = moisture_flow.isel(
        targetlat=target_latitude_index, targetlon=target_longitude_index
    ).load()
    flow_fraction = decode_flow(encoded_target_flow, flow_encoding)
    source_evaporation = evaporation_volume.rename({"lat": "sourcelat", "lon": "sourcelon"}).load()
    contribution = (flow_fraction * source_evaporation).rename("evaporation_contribution_m3")
    contribution_total = contribution.sum(skipna=True)
    contribution_share = xr.where(
        contribution_total > 0,
        contribution / contribution_total,
        np.nan,
    ).rename("fraction_of_utrack_target_precipitation")

    target_latitude = float(moisture_flow.targetlat.isel(targetlat=target_latitude_index))
    target_longitude = float(moisture_flow.targetlon.isel(targetlon=target_longitude_index))
    target_longitude_wgs84 = (target_longitude + 180) % 360 - 180
    result = xr.Dataset(
        data_vars={
            "evaporation_contribution_m3": contribution,
            "fraction_of_utrack_target_precipitation": contribution_share,
            "moisture_flow_fraction": flow_fraction.rename("moisture_flow_fraction"),
            "source_evaporation_m3": source_evaporation.rename("source_evaporation_m3"),
        },
        attrs={
            "title": "Monthly UTrack precipitationshed",
            "requested_target_latitude_degrees_north": latitude,
            "requested_target_longitude_degrees_east": (longitude + 180) % 360 - 180,
            "selected_target_latitude_degrees_north": target_latitude,
            "selected_target_longitude_degrees_east": target_longitude_wgs84,
            "flow_encoding": flow_encoding,
            "description": (
                "Evaporation contributions from source cells to the selected target cell. "
                "The contribution equals UTrack source-to-target flow fraction times monthly ERA5_ET volume."
            ),
        },
    )
    result["evaporation_contribution_m3"].attrs.update(
        long_name="Monthly source evaporation contribution to selected target precipitation",
        units="m3 month-1",
    )
    result["fraction_of_utrack_target_precipitation"].attrs.update(
        long_name="Source contribution share of UTrack-derived precipitation at selected target",
        units="1",
    )
    result["moisture_flow_fraction"].attrs.update(
        long_name="UTrack source-to-target moisture flow fraction",
        units="1",
    )
    result["source_evaporation_m3"].attrs.update(long_name="Monthly ERA5 source evapotranspiration volume", units="m3 month-1")
    return result


def main() -> None:
    arguments = parse_arguments()
    arguments.output_directory.mkdir(parents=True, exist_ok=True)

    for month in range(1, 13):
        utrack_path, era5_path = monthly_paths(arguments.utrack_directory, arguments.era5_directory, month)
        with xr.open_dataset(utrack_path) as utrack_dataset, xr.open_dataset(era5_path) as era5_volumes:
            moisture_flow = utrack_dataset["moisture_flow"]
            validate_inputs(moisture_flow, era5_volumes, utrack_path, era5_path)
            precipitationshed = calculate_monthly_precipitationshed(
                moisture_flow,
                era5_volumes["ERA5_ET"],
                arguments.latitude,
                arguments.longitude,
                arguments.flow_encoding,
            )

        output_path = arguments.output_directory / f"utrack_precipitationshed_month_{month:02d}.nc"
        precipitationshed.to_netcdf(output_path)
        total_contribution = float(precipitationshed["evaporation_contribution_m3"].sum(skipna=True))
        print(f"Month {month:02d}: wrote {output_path} ({total_contribution:.3e} m3 month-1)")


if __name__ == "__main__":
    main()
