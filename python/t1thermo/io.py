"""Loading of configuration files, NIfTI volumes and the fiber-optic probe trace."""

import json
from pathlib import Path

import nibabel as nib
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]


def load_config(path: str | Path) -> dict:
    """Read a JSON configuration; relative paths inside it are resolved against the repository root."""
    cfg = json.loads(Path(path).read_text(encoding="utf-8"))
    for key in ("data_dir", "output_dir"):
        if key in cfg:
            cfg[key] = (REPO_ROOT / cfg[key]).resolve()
    return cfg


def load_nifti(path: str | Path) -> np.ndarray:
    return np.asarray(nib.load(str(path)).dataobj, dtype=float)


def save_nifti(data: np.ndarray, path: str | Path, reference: str | Path | None = None) -> None:
    affine = nib.load(str(reference)).affine if reference is not None else np.eye(4)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(np.asarray(data, dtype=np.float32), affine), str(path))


def load_signal_change(path: str | Path) -> np.ndarray:
    """Relative T1-weighted signal change dS/S0 from stage s05, returned as (X, Y, T) of the first slice."""
    data = load_nifti(path)
    if data.ndim == 4:  # (X, Y, slice, T)
        data = data[:, :, 0, :]
    return data


def load_prfs(path: str | Path) -> np.ndarray:
    """PRFS temperature change from stage s04 (°C), returned as (X, Y, T) of the first slice."""
    data = load_nifti(path)
    return data[:, :, 0, :] if data.ndim == 4 else data


def load_muscle_mask(path: str | Path, shape_xy: tuple[int, int]) -> np.ndarray:
    """Muscle mask from stage s04 (masks.nii.gz stacks body, muscle and fat masks along the last axis)."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Muscle mask not found: {path}")
    m = load_nifti(path)
    if m.ndim == 3 and m.shape[2] > 1:
        m = m[:, :, 1]
    m = np.squeeze(m)
    if m.shape != tuple(shape_xy):
        raise ValueError(f"Mask shape {m.shape} does not match image shape {shape_xy}")
    return np.isfinite(m) & (m > 0)


def load_probe_trace(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Fiber-optic probe file: header lines starting with '#' give the sampling period and the probe
    sample at which the first MRI dynamic started. Returns time relative to MRI start (s), probe 1, probe 2.
    """
    meta: dict[str, float] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("#") and ":" in line:
            key, value = line[1:].split(":", 1)
            try:
                meta[key.strip()] = float(value.split()[0])
            except (ValueError, IndexError):
                pass
    data = np.loadtxt(path, comments="#")
    dt, start = meta["sampling_period_s"], meta["mri_start_sample"]
    time_s = (np.arange(data.shape[0]) - start) * dt
    return time_s, data[:, 0], data[:, 1]
