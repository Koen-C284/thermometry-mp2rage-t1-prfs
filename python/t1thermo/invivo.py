"""
In vivo calibration against PRFS (Figs. 7, 8, 9B) and cross-subject transferability.

The reference subject is calibrated against its own PRFS; its coefficient is applied unchanged to the
other subjects. In addition, every subject is calibrated independently and each coefficient is scored on
every subject (RMSE against PRFS after the baseline period).
"""

import json
from pathlib import Path

import numpy as np

from . import figures
from .calibration import apply_calibration, calibrate_against_prfs, transfer_matrix
from .io import load_muscle_mask, load_prfs, load_signal_change, save_nifti


def load_subject(folder: Path, tag: str) -> dict:
    prfs = load_prfs(folder / "temp_prfs.nii.gz")
    signal = load_signal_change(folder / f"signal_change_{tag}.nii.gz")
    n = min(prfs.shape[-1], signal.shape[-1])
    return {"signal": signal[..., :n], "prfs": prfs[..., :n],
            "muscle": load_muscle_mask(folder / "masks.nii.gz", prfs.shape[:2])}


def run(cfg: dict) -> dict:
    data, out = Path(cfg["data_dir"]), Path(cfg["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    skip, n_base, dt = cfg["discard_dynamics"], cfg["baseline_dynamics"], cfg["dynamic_interval_s"]
    reference = cfg["reference_subject"]
    results: dict = {"reference_subject": reference}

    for tag in ("img1", "img2"):
        subjects = {name: load_subject(data / name, tag) for name in cfg["subjects"]}
        for name, s in subjects.items():
            cal = calibrate_against_prfs(s["signal"][..., skip:], s["prfs"][..., skip:], s["muscle"],
                                         cfg["min_heating_c"], cfg["min_abs_r"], max(0, n_base - skip))
            s.update(alpha=cal.alpha, pool=cal.pool_mask, cal=cal)
            figures.calibration_scatter(s["signal"][..., skip:], s["prfs"][..., skip:], cal, s["muscle"],
                                        max(0, n_base - skip), cfg["min_heating_c"],
                                        out / f"fig9b_calibration_{name}_{tag}.png", f"{name} {tag}")

        alpha_ref = subjects[reference]["alpha"]
        for name, s in subjects.items():
            time_s = np.arange(s["signal"].shape[-1]) * dt
            temp = apply_calibration(s["signal"], alpha_ref, n_base, s["muscle"])
            save_nifti(temp[:, :, None, :], out / f"temperature_t1_{name}_{tag}.nii.gz", data / name / "temp_prfs.nii.gz")
            spec = cfg["subjects"][name]
            positions = {f"P{k + 1}": tuple(v) for k, v in enumerate(spec.get("voxels", []))}
            if positions:
                figures.voxel_timecourses(time_s, s["prfs"], temp, positions, cfg["phases_s"],
                                          out / f"fig7_timecourses_{name}_{tag}.png", f"{name} ({tag})")
            figures.period_maps(s["prfs"], temp, s["muscle"], time_s, {"Peak heating": spec["peak_window_s"]},
                                positions, out / f"fig8_peak_maps_{name}_{tag}.png")

        scores = transfer_matrix(subjects, n_base)
        results[tag] = {
            "coefficients": {n: {"alpha_weighted": s["alpha"], "alpha_pooled": s["cal"].pooled.alpha,
                                 "pooled_r2_uncentered": s["cal"].pooled.r2_uncentered,
                                 "n_calibration_voxels": s["cal"].n_pool, "n_muscle_voxels": s["cal"].n_muscle}
                             for n, s in subjects.items()},
            "transfer_rmse": {f"{src} -> {tgt}": v for (src, tgt), v in scores.items()},
        }

    (out / "invivo_results.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    return results
