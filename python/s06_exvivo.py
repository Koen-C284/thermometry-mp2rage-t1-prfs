"""
Stages s06-s07, ex vivo: probe calibration, validation against both fiber-optic probes, precision and the
paper figures (Figs. 5, 6, 9A, Supplementary Fig. 2).

Usage (from the repository root):  python python/s06_exvivo.py [config/exvivo.json]
Outputs: results/exvivo/
"""

import json
import sys
import warnings

from t1thermo import exvivo
from t1thermo.io import REPO_ROOT, load_config

if __name__ == "__main__":
    warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN background voxels outside the muscle mask
    cfg = load_config(sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "config" / "exvivo.json")
    print(json.dumps(exvivo.run(cfg), indent=2))
    print(f"Outputs written to {cfg['output_dir']}")
