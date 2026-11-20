"""Zero-intercept (regression through the origin) fit used for every calibration in the paper."""

from dataclasses import dataclass

import numpy as np


@dataclass
class OriginFit:
    """Result of fitting y = alpha * x."""

    alpha: float
    r2_uncentered: float
    rmse: float
    n: int


def fit_through_origin(x: np.ndarray, y: np.ndarray) -> OriginFit:
    """
    Least-squares fit of y = alpha * x, ignoring non-finite pairs.

    The uncentered R^2 = 1 - SS_res / sum(y^2) is the appropriate goodness of fit for a model
    without intercept (it is bounded in [0, 1]).
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    sxx = float(np.sum(x * x))
    if x.size < 2 or sxx == 0.0:
        return OriginFit(float("nan"), float("nan"), float("nan"), int(x.size))

    alpha = float(np.sum(x * y) / sxx)
    residual = y - alpha * x
    syy = float(np.sum(y * y))
    r2 = 1.0 - float(np.sum(residual**2)) / syy if syy > 0 else float("nan")
    return OriginFit(alpha, r2, float(np.sqrt(np.mean(residual**2))), int(x.size))


def error_metrics(estimate: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    """Mean error (bias), RMSE and SD of the error of an estimate against a reference (°C)."""
    e = np.asarray(estimate, dtype=float) - np.asarray(reference, dtype=float)
    e = e[np.isfinite(e)]
    if e.size == 0:
        return {"mean_error": float("nan"), "rmse": float("nan"), "sd_error": float("nan"), "n": 0}
    return {
        "mean_error": float(np.mean(e)),
        "rmse": float(np.sqrt(np.mean(e**2))),
        "sd_error": float(np.std(e, ddof=1)) if e.size > 1 else float("nan"),
        "n": int(e.size),
    }


def quadratic_precision(time_s: np.ndarray, temperature: np.ndarray) -> tuple[float, np.ndarray]:
    """
    Precision of a smoothly varying temperature trajectory: SD of the residuals of a quadratic fit
    in time (Supplementary Fig. 2). Returns the SD and the fitted curve.
    """
    t = np.asarray(time_s, dtype=float)
    y = np.asarray(temperature, dtype=float)
    ok = np.isfinite(t) & np.isfinite(y)
    fitted = np.full_like(y, np.nan)
    if ok.sum() < 4:
        return float("nan"), fitted
    coef = np.polyfit(t[ok], y[ok], 2)
    fitted = np.polyval(coef, t)
    return float(np.std(y[ok] - fitted[ok], ddof=3)), fitted
