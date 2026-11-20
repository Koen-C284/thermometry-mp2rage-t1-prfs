"""
Stage s02, in vivo only: motion correction of the dynamic series and alignment of the B1/T1 maps (elastix).

Usage:  python python/s02_registration.py <input folder> <output folder> [--grid-px 5] [--samples 2048]
The input folder holds the stage s01 outputs (mp2rage_image.nii.gz + .json, t1_map[_image].nii.gz,
b1_map[_image].nii.gz). Requires itk-elastix. Defaults follow the paper (5-pixel B-spline grid, 2048 samples).
"""

import argparse
from pathlib import Path

from t1thermo.registration import main

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Stage s02: in vivo motion correction")
    parser.add_argument("src", type=Path)
    parser.add_argument("dst", type=Path)
    parser.add_argument("--grid-px", type=float, default=5.0, help="final B-spline grid spacing in pixels")
    parser.add_argument("--samples", type=int, default=2048, help="number of random samples for the metric")
    args = parser.parse_args()
    main(args.src, args.dst, args.grid_px, args.samples)
