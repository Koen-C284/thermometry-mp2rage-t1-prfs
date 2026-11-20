"""
PRFS-referenced calibration of the T1-based signal change (in vivo procedure of the paper).

Model per voxel: dT = alpha * (dS / S0), with zero intercept. A voxel joins the calibration pool if its
PRFS temperature rises by at least `min_heating_c` and |Pearson r| between signal change and PRFS
is at least `min_abs_r`. The applied coefficient is the R^2-weighted mean of the voxel slopes.
"""

from dataclasses import dataclass

import numpy as np

from .regression import OriginFit, error_metrics, fit_through_origin


@dataclass
class Calibration:
    alpha: float  # applied coefficient (R^2-weighted mean of voxel slopes), °C per unit dS/S0
    pooled: OriginFit  # single fit through all calibration samples (Fig. 9 line)
    slope_map: np.ndarray  # (X, Y) voxelwise zero-intercept slope
    r_map: np.ndarray  # (X, Y) voxelwise Pearson r
    pool_mask: np.ndarray  # (X, Y) calibration voxels
    n_pool: int
    n_muscle: int


def voxelwise_fit(signal: np.ndarray, prfs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Zero-intercept slope and Pearson r per voxel over the time axis. Inputs are (X, Y, T)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        sxx = np.nansum(signal**2, axis=-1)
        slope = np.where(sxx > 0, np.nansum(signal * prfs, axis=-1) / sxx, np.nan)
        ds = signal - np.nanmean(signal, axis=-1, keepdims=True)
        dp = prfs - np.nanmean(prfs, axis=-1, keepdims=True)
        vs, vp = np.nansum(ds**2, axis=-1), np.nansum(dp**2, axis=-1)
        r = np.where((vs > 0) & (vp > 0), np.nansum(ds * dp, axis=-1) / np.sqrt(vs * vp), np.nan)
    return slope, r


def calibrate_against_prfs(
    signal: np.ndarray,
    prfs: np.ndarray,
    muscle: np.ndarray,
    min_heating_c: float = 0.5,
    min_abs_r: float = 0.9,
    pooled_exclude_first: int = 0,
) -> Calibration:
    """
    Calibrate dS/S0 against PRFS. `signal` and `prfs` are (X, Y, T) over the dynamics used for calibration.
    For the pooled fit, the first `pooled_exclude_first` dynamics (baseline) and samples with
    PRFS < `min_heating_c` are left out, as in Fig. 9.
    """
    slope, r = voxelwise_fit(signal, prfs)
    with np.errstate(invalid="ignore"):
        heated = np.nanmax(prfs, axis=-1) >= min_heating_c
    pool = muscle & heated & (np.abs(r) >= min_abs_r) & np.isfinite(slope)

    weights = np.where(pool, r**2, 0.0)
    alpha = float(np.sum(slope[pool] * weights[pool]) / np.sum(weights[pool])) if pool.any() else float("nan")

    x = signal[pool][:, pooled_exclude_first:].ravel()
    y = prfs[pool][:, pooled_exclude_first:].ravel()
    keep = y >= min_heating_c
    pooled = fit_through_origin(x[keep], y[keep])
    return Calibration(alpha, pooled, slope, r, pool, int(pool.sum()), int(muscle.sum()))


def apply_calibration(signal: np.ndarray, alpha: float, baseline_dynamics: int, mask: np.ndarray | None = None) -> np.ndarray:
    """Temperature change map alpha * dS/S0, zero-aligned to the mean of the first `baseline_dynamics` dynamics."""
    temp = alpha * signal
    n = min(max(1, baseline_dynamics), temp.shape[-1])
    with np.errstate(invalid="ignore"):
        temp = temp - np.nanmean(temp[..., :n], axis=-1, keepdims=True)
    if mask is not None:
        temp = np.where(mask[..., None], temp, np.nan)
    return temp


def transfer_matrix(
    subjects: dict[str, dict],
    baseline_dynamics: int,
) -> dict[tuple[str, str], dict[str, dict[str, float]]]:
    """
    Apply every subject's coefficient to every subject and score against that subject's PRFS after the
    baseline period, on its calibration voxels and on all its muscle voxels.

    `subjects[name]` holds 'alpha', 'signal' (X, Y, T, full series), 'prfs' (X, Y, T), 'muscle', 'pool'.
    """
    scores = {}
    for src, s in subjects.items():
        for tgt, t in subjects.items():
            temp = apply_calibration(t["signal"], s["alpha"], baseline_dynamics, t["muscle"])
            prfs = t["prfs"][..., : temp.shape[-1]]
            scores[(src, tgt)] = {
                region: error_metrics(temp[mask][:, baseline_dynamics:], prfs[mask][:, baseline_dynamics:])
                for region, mask in (("calibration_voxels", t["pool"]), ("muscle", t["muscle"]))
            }
    return scores
