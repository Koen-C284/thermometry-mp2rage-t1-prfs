"""
Stages s06-s07, in vivo: PRFS-referenced calibration, transfer of the reference coefficient, cross-subject
transferability and the paper figures (Figs. 7, 8, 9B). In vivo data are not distributed; see config/invivo.json.

Usage (from the repository root):  python python/s06_invivo.py [config/invivo.json]
Outputs: results/invivo/
"""

import json
import sys
import warnings

from t1thermo import invivo
from t1thermo.io import REPO_ROOT, load_config

if __name__ == "__main__":
    warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN background voxels outside the muscle mask
    cfg = load_config(sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "config" / "invivo.json")
    print(json.dumps(invivo.run(cfg), indent=2))
    print(f"Outputs written to {cfg['output_dir']}")
