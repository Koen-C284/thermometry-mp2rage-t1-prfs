"""Figures of the paper (ex vivo Figs. 5, 6, 9A, Supp. Fig. 2; in vivo Figs. 7, 8, 9B)."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .regression import quadratic_precision

COLORS = {"calibration": "#0072B2", "validation": "#D55E00"}
POSITION_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def calibration_scatter(signal, prfs, cal, muscle, exclude_first, min_heating_c, path, title):
    """Calibration voxels (blue) vs all other muscle samples (grey) with the zero-intercept fit (Fig. 9)."""
    x_c = signal[cal.pool_mask][:, exclude_first:].ravel()
    y_c = prfs[cal.pool_mask][:, exclude_first:].ravel()
    keep = y_c >= min_heating_c
    other = muscle & ~cal.pool_mask
    x_o = np.concatenate([signal[other].ravel(), x_c[~keep], signal[cal.pool_mask][:, :exclude_first].ravel()])
    y_o = np.concatenate([prfs[other].ravel(), y_c[~keep], prfs[cal.pool_mask][:, :exclude_first].ravel()])

    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.scatter(x_o, y_o, s=4, c="#9AA0A6", alpha=0.3, lw=0, label="Other muscle samples")
    ax.scatter(x_c[keep], y_c[keep], s=6, c="#0072B2", alpha=0.6, lw=0, label=f"Calibration voxels (n={cal.n_pool})")
    finite = np.isfinite(x_o)
    xs = np.linspace(np.nanpercentile(np.r_[x_o[finite], x_c], 1), np.nanpercentile(np.r_[x_o[finite], x_c], 99), 50)
    ax.plot(xs, cal.pooled.alpha * xs, c="#E69F00", lw=2.5,
            label=f"dT = {cal.pooled.alpha:.1f} dS/S0 (R² = {cal.pooled.r2_uncentered:.2f})")
    ax.set(xlabel="Relative signal change dS/S0", ylabel="PRFS temperature change (°C)", title=title)
    ax.grid(ls="--", alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    _save(fig, path)


PROBE_STYLE = {  # (probe, T1-based, PRFS) colours per position, as in Fig. 5
    0: ("#1F6FD6", "#0B3D91", "#4BA3F2"),
    1: ("#E8590C", "#B23A14", "#FFA94D"),
}


def _crop(images, mask, pad=2):
    """Bounding box of the mask (with padding), applied to every image."""
    rows, cols = np.where(mask)
    r0, r1 = max(rows.min() - pad, 0), min(rows.max() + pad + 1, mask.shape[0])
    c0, c1 = max(cols.min() - pad, 0), min(cols.max() + pad + 1, mask.shape[1])
    return [im[r0:r1, c0:c1] for im in images], (r0, c0)


def probe_timecourses(time_s, probes, prfs_ts, t1_ts, positions, underlay, mask, path):
    """Voxel locations and probe, T1-based and PRFS temperature at P1 and P2 (Fig. 5)."""
    fig = plt.figure(figsize=(13.5, 5))
    ax_img = fig.add_axes([0.0, 0.15, 0.22, 0.7])
    ax = fig.add_axes([0.28, 0.12, 0.7, 0.76])

    (img,), (r0, c0) = _crop([underlay], mask)
    ax_img.imshow(img, cmap="gray")
    ax_img.set_title("Selected Voxel Locations", fontsize=12, fontweight="bold")
    ax_img.axis("off")

    for k, name in enumerate(probes):
        c_probe, c_t1, c_prfs = PROBE_STYLE[k % 2]
        x, y = positions[name]
        ax_img.text(y - c0, x - r0, f" P{k + 1} ", color="white", fontsize=11, fontweight="bold", ha="center",
                    va="center", bbox=dict(boxstyle="round,pad=0.2", fc=c_t1, ec="k", lw=0.8))
        ax.plot(time_s, probes[name], c=c_probe, lw=2.2, label=f"Fiber probe P{k + 1}")
        ax.plot(time_s, t1_ts[name], c=c_t1, lw=1.6, label=f"P{k + 1} T$_1$-Based")
        ax.plot(time_s, prfs_ts[name], c=c_prfs, lw=1.6, label=f"P{k + 1} PRFS")

    dt = float(np.median(np.diff(time_s)))
    ax.set_xlabel("Acquisition Dynamic", fontsize=13, fontweight="bold")
    ax.set_ylabel("ΔT (°C)", fontsize=13, fontweight="bold")
    ticks = np.arange(0, len(time_s) + 50, 50)
    ax.set_xticks(ticks * dt, [str(t) for t in ticks])
    top = ax.secondary_xaxis("top")
    top.set_xlabel("Time (s)", fontsize=13, fontweight="bold")
    ax.grid(ls="--", alpha=0.4)
    handles, labels = ax.get_legend_handles_labels()
    order = [0, 3, 1, 4, 2, 5] if len(handles) == 6 else range(len(handles))
    ax.legend([handles[i] for i in order], [labels[i] for i in order], ncol=3, fontsize=10, loc="lower right")
    _save(fig, path)


def period_maps(prfs, t1, muscle, time_s, periods, positions, path, vmin=0.0, vmax=4.0, underlay=None):
    """Temperature maps averaged over time windows: PRFS (top) and T1-based (bottom) (Figs. 6 and 8)."""
    fig, axes = plt.subplots(2, len(periods), figsize=(3.7 * len(periods), 6.4), squeeze=False)
    bg = np.zeros(muscle.shape) if underlay is None else underlay
    for j, (label, (t0, t1_end)) in enumerate(periods.items()):
        sel = (time_s >= t0) & (time_s < t1_end)
        for i, (name, vol) in enumerate((("PRFS", prfs), ("T$_1$-Based", t1))):
            temp = np.where(muscle, np.nanmean(vol[..., sel], axis=-1), np.nan)
            (b, t), (r0, c0) = _crop([bg, temp], muscle)
            ax = axes[i, j]
            ax.imshow(b, cmap="gray")
            im = ax.imshow(t, cmap="afmhot", vmin=vmin, vmax=vmax)
            for k, (x, y) in enumerate(positions.values()):
                c = PROBE_STYLE[k % 2][0]
                ax.plot(y - c0, x - r0, "s", mfc="none", mec=c, ms=5, mew=1.5)
                ax.text(y - c0 + 1.5, x - r0, f"P{k + 1}", color=c, fontsize=11, fontweight="bold", va="center")
            ax.axis("off")
            if i == 0:
                ax.set_title(label, fontsize=13, fontweight="bold")
            if j == 0:
                ax.text(-0.08, 0.5, name, transform=ax.transAxes, rotation=90, fontsize=14,
                        fontweight="bold", ha="center", va="center")
    cbar = fig.colorbar(im, ax=axes, orientation="horizontal", fraction=0.05, pad=0.06, shrink=0.5)
    cbar.set_label("ΔT (°C)", fontsize=12, fontweight="bold")
    _save(fig, path)


def precision_plot(time_s, series_by_method, path):
    """Quadratic fit to each temperature trajectory with the ±1 SD residual band (Supplementary Fig. 2)."""
    methods = list(series_by_method)
    positions = list(next(iter(series_by_method.values())))
    fig, axes = plt.subplots(len(positions), len(methods), figsize=(4 * len(methods), 3.2 * len(positions)), squeeze=False)
    for i, pos in enumerate(positions):
        for j, method in enumerate(methods):
            y = series_by_method[method][pos]
            sd, fit = quadratic_precision(time_s, y)
            ax = axes[i, j]
            ax.plot(time_s, y, ".", ms=2, c="#555555")
            ax.plot(time_s, fit, c="#D55E00", lw=1.5)
            ax.fill_between(time_s, fit - sd, fit + sd, color="#D55E00", alpha=0.2)
            ax.set_title(f"{method}, {pos}: SD = {sd:.3f} °C", fontsize=9)
            ax.grid(ls="--", alpha=0.3)
    _save(fig, path)


def voxel_timecourses(time_s, prfs, t1, positions, phases_s, path, title=""):
    """PRFS (solid) and T1-based (dashed) temperature at selected voxels (Fig. 7)."""
    fig, ax = plt.subplots(figsize=(8, 4))
    for name, (start, end) in phases_s.items():
        ax.axvspan(start, end, color="#F2F2F2" if name != "Heating" else "#FDE7D9", zorder=0)
    for k, (name, (x, y)) in enumerate(positions.items()):
        c = POSITION_COLORS[k % 4]
        ax.plot(time_s, prfs[x, y], c=c, lw=1.4, label=f"{name} PRFS")
        ax.plot(time_s, t1[x, y], c=c, lw=1.4, ls="--", label=f"{name} T1-based")
    ax.set(xlabel="Time (s)", ylabel="Temperature change (°C)", title=title)
    ax.grid(ls="--", alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    _save(fig, path)
