"""Calculate monthly UTrack precipitationsheds for all cells in a polygon grid."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
ENGINE_PATH = SCRIPT_DIRECTORY / "UTrack-Polygon_01_MonthlyP+Eshed.py"


def load_engine():
    spec = importlib.util.spec_from_file_location("utrack_polygon_engine", ENGINE_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load polygon calculation engine: {ENGINE_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polygon-grid-file", type=Path, required=True)
    parser.add_argument("--utrack-directory", type=Path, required=True)
    parser.add_argument("--era5-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--flow-encoding", choices=("encoded-log", "fraction"), default="encoded-log")
    parser.add_argument("--batch-size", type=int, default=16)
    arguments = parser.parse_args()
    if arguments.batch_size < 1:
        parser.error("--batch-size must be positive.")
    return arguments


def main() -> None:
    arguments = parse_arguments()
    engine = load_engine()
    arguments.output_directory.mkdir(parents=True, exist_ok=True)
    for month in range(1, 13):
        utrack_path, era5_path = engine.monthly_paths(arguments.utrack_directory, arguments.era5_directory, month)
        with engine.xr.open_dataset(utrack_path) as utrack, engine.xr.open_dataset(era5_path) as era5:
            flow = utrack["moisture_flow"]
            engine.validate_inputs(flow, era5)
            latitude_indices, longitude_indices = engine.polygon_cell_indices(arguments.polygon_grid_file, flow)
            precipitationshed = engine.aggregate_precipitationshed(
                flow,
                era5["ERA5_ET"],
                latitude_indices,
                longitude_indices,
                arguments.batch_size,
                arguments.flow_encoding,
            ).to_dataset()
            precipitationshed.attrs = {
                "title": "Monthly UTrack polygon precipitationshed",
                "polygon_grid_file": str(arguments.polygon_grid_file),
                "selected_grid_cell_count": len(latitude_indices),
                "flow_encoding": arguments.flow_encoding,
            }
        output_path = arguments.output_directory / f"utrack_polygon_precipitationshed_month_{month:02d}.nc"
        precipitationshed.to_netcdf(output_path)
        print(f"Month {month:02d}: wrote {output_path} for {len(latitude_indices)} grid cells.")


if __name__ == "__main__":
    main()
