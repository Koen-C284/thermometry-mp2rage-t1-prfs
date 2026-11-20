"""
Ex vivo validation against fiber-optic probes (Figs. 5, 6, 9A and Supplementary Fig. 2).

The T1-based signal change at the probe-1 position is calibrated against probe 1 (zero intercept);
probe 2 is an independent position used only for validation. Probe, PRFS and T1-based series are all
zero-referenced to the same baseline dynamics and aligned on the acquisition time of each dynamic.
"""

import json
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy import ndimage

from . import figures
from .calibration import calibrate_against_prfs
from .io import load_muscle_mask, load_prfs, load_probe_trace, load_signal_change, save_nifti
from .regression import error_metrics, fit_through_origin, quadratic_precision


def zero_reference(x: np.ndarray, n: int) -> np.ndarray:
    """Subtract the mean of the first n samples along the last axis."""
    with np.errstate(invalid="ignore"):
        return x - np.nanmean(x[..., :n], axis=-1, keepdims=True)


def smooth_in_mask(images: np.ndarray, mask: np.ndarray, sigma: float, size: int = 5) -> np.ndarray:
    """Gaussian smoothing (size x size kernel) restricted to the mask, ignoring voxels outside it."""
    ax = np.arange(size) - size // 2
    k1 = np.exp(-(ax**2) / (2 * sigma**2))
    kernel = np.outer(k1, k1)
    out = np.empty_like(images)
    for t in range(images.shape[-1]):
        inside = mask & np.isfinite(images[..., t])
        weight = ndimage.convolve(inside.astype(float), kernel, mode="constant")
        total = ndimage.convolve(np.where(inside, images[..., t], 0.0), kernel, mode="constant")
        with np.errstate(invalid="ignore", divide="ignore"):
            out[..., t] = np.where(mask & (weight > 0), total / weight, np.nan)
    return out


def run(cfg: dict) -> dict:
    data, out = Path(cfg["data_dir"]), Path(cfg["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    skip, n_base = cfg["discard_dynamics"], cfg["baseline_dynamics"]
    cal_xy, val_xy = tuple(cfg["probe_voxels"]["calibration"]), tuple(cfg["probe_voxels"]["validation"])

    sidecar = json.loads((data / "mp2rage_image.json").read_text(encoding="utf-8"))
    dyn_t = np.asarray(sidecar["dynamic_start_times_s"], dtype=float)[skip:]
    time_s = dyn_t - dyn_t[0]

    probe_t, p1, p2 = load_probe_trace(data / cfg["probe_file"])
    probes = {"calibration": zero_reference(np.interp(dyn_t, probe_t, p1), n_base),
              "validation": zero_reference(np.interp(dyn_t, probe_t, p2), n_base)}

    prfs = zero_reference(load_prfs(data / "derived" / "temp_prfs.nii.gz")[..., skip:], n_base)
    muscle = load_muscle_mask(data / "derived" / "masks.nii.gz", prfs.shape[:2])
    voxels = {"calibration": cal_xy, "validation": val_xy}

    results: dict = {"settings": {k: cfg[k] for k in ("discard_dynamics", "baseline_dynamics", "probe_voxels")}}
    prfs_ts = {name: prfs[xy] for name, xy in voxels.items()}
    results["prfs"] = {name: error_metrics(prfs_ts[name][n_base:], probes[name][n_base:]) for name in voxels}
    results["prfs"]["precision_sd"] = {name: quadratic_precision(time_s, prfs_ts[name])[0] for name in voxels}

    temps = {}
    for tag in ("img1", "img2"):
        signal = zero_reference(load_signal_change(data / "derived" / f"signal_change_{tag}.nii.gz")[..., skip:], n_base)
        fit = fit_through_origin(signal[cal_xy][n_base:], probes["calibration"][n_base:])
        temp = np.where(muscle[..., None], fit.alpha * signal, np.nan)
        temps[tag] = temp

        smoothed = fit.alpha * smooth_in_mask(signal, muscle, cfg["filtered_sigma_px"])
        res = {
            "alpha": fit.alpha,
            "fit_r2_uncentered": fit.r2_uncentered,
            "calibration_probe": error_metrics(temp[cal_xy][n_base:], probes["calibration"][n_base:]),
            "validation_probe": error_metrics(temp[val_xy][n_base:], probes["validation"][n_base:]),
            "precision_sd": {name: quadratic_precision(time_s, temp[xy])[0] for name, xy in voxels.items()},
            "precision_sd_filtered": {name: quadratic_precision(time_s, smoothed[xy])[0] for name, xy in voxels.items()},
        }

        # Same PRFS-referenced voxel selection as in vivo, for Fig. 9A (not used for the probe calibration)
        prfs_cal = calibrate_against_prfs(signal, prfs, muscle, cfg["min_heating_c"], cfg["min_abs_r"], n_base)
        res["prfs_referenced"] = {"alpha_weighted": prfs_cal.alpha, "alpha_pooled": prfs_cal.pooled.alpha,
                                  "pooled_r2_uncentered": prfs_cal.pooled.r2_uncentered,
                                  "n_voxels": prfs_cal.n_pool, "n_muscle": prfs_cal.n_muscle}
        figures.calibration_scatter(signal, prfs, prfs_cal, muscle, n_base, cfg["min_heating_c"],
                                    out / f"fig9a_calibration_{tag}.png", f"Ex vivo {tag.upper().replace('IMG', 'GRE ')}")
        save_nifti(temp[:, :, None, :], out / f"temperature_t1_{tag}.nii.gz", data / "derived" / "temp_prfs.nii.gz")
        results[tag] = res

    # The paper's T1-based curves and maps (Figs. 5 and 6) are the mean of the GRE 1 and GRE 2 temperatures
    temp_mean = np.nanmean(np.stack([temps["img1"], temps["img2"]]), axis=0)
    temps["mean"] = temp_mean
    results["mean_gre1_gre2"] = {
        "calibration_probe": error_metrics(temp_mean[cal_xy][n_base:], probes["calibration"][n_base:]),
        "validation_probe": error_metrics(temp_mean[val_xy][n_base:], probes["validation"][n_base:]),
        "precision_sd": {name: quadratic_precision(time_s, temp_mean[xy])[0] for name, xy in voxels.items()},
    }
    save_nifti(temp_mean[:, :, None, :], out / "temperature_t1_mean.nii.gz", data / "derived" / "temp_prfs.nii.gz")

    underlay = np.asarray(nib.load(str(data / "mp2rage_image.nii.gz")).dataobj[:, :, 0, 0, :10, 1, 0], dtype=float).mean(-1)
    figures.probe_timecourses(time_s, probes, prfs_ts, {n: temp_mean[xy] for n, xy in voxels.items()}, voxels,
                              underlay, muscle, out / "fig5_probe_timecourses.png")
    figures.period_maps(prfs, temp_mean, muscle, time_s, cfg["map_periods_s"], voxels, out / "fig6_temperature_maps.png",
                        underlay=underlay)
    figures.precision_plot(time_s, {f"T1 {t}": {n: temps[t][xy] for n, xy in voxels.items()} for t in ("img1", "img2")} |
                           {"PRFS": prfs_ts}, out / "suppfig2_precision.png")

    (out / "exvivo_results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    return results
