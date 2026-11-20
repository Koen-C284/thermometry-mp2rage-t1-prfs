"""Unit tests for the calibration and validation maths (run: python -m pytest python/tests)."""

import numpy as np
import pytest

from t1thermo.calibration import apply_calibration, calibrate_against_prfs, transfer_matrix
from t1thermo.exvivo import smooth_in_mask, zero_reference
from t1thermo.regression import error_metrics, fit_through_origin, quadratic_precision


def synthetic(alpha=-50.0, n_dyn=120, n_base=20, noise=0.0, seed=0):
    rng = np.random.default_rng(seed)
    x = np.zeros(n_dyn)
    x[n_base:] = np.linspace(-0.005, -0.06, n_dyn - n_base)
    signal = np.broadcast_to(x, (6, 6, n_dyn)).copy() + noise * rng.standard_normal((6, 6, n_dyn))
    prfs = alpha * np.broadcast_to(x, (6, 6, n_dyn)).copy()
    return signal, prfs, np.ones((6, 6), bool)


def test_fit_through_origin_recovers_slope():
    x = np.linspace(-1, 1, 50)
    fit = fit_through_origin(x, 3.0 * x)
    assert fit.alpha == pytest.approx(3.0)
    assert fit.r2_uncentered == pytest.approx(1.0)
    assert fit.rmse == pytest.approx(0.0, abs=1e-12)


def test_error_metrics_and_precision():
    m = error_metrics(np.array([1.0, 2.0, 3.0]), np.array([0.5, 1.5, 2.5]))
    assert m["mean_error"] == pytest.approx(0.5) and m["rmse"] == pytest.approx(0.5)
    t = np.linspace(0, 100, 200)
    sd, fit = quadratic_precision(t, 1e-4 * t**2 + 0.01 * t)
    assert sd == pytest.approx(0.0, abs=1e-9)


def test_calibration_selects_heated_correlated_voxels():
    signal, prfs, muscle = synthetic()
    prfs[0, 0] = 0.0  # unheated voxel must be excluded
    cal = calibrate_against_prfs(signal, prfs, muscle, 0.5, 0.9, pooled_exclude_first=20)
    assert cal.alpha == pytest.approx(-50.0)
    assert cal.pooled.alpha == pytest.approx(-50.0)
    assert not cal.pool_mask[0, 0] and cal.n_pool == 35


def test_apply_calibration_zero_aligns_baseline():
    signal, _, muscle = synthetic()
    temp = apply_calibration(signal + 0.01, -50.0, 20, muscle)
    assert np.allclose(np.nanmean(temp[..., :20], axis=-1), 0.0)


def test_transfer_matrix_scales_with_coefficient_ratio():
    sig_a, prfs_a, m = synthetic(-50.0)
    sig_b, prfs_b, _ = synthetic(-60.0)
    subjects = {n: {"alpha": a, "signal": s, "prfs": p, "muscle": m, "pool": m}
                for n, a, s, p in (("a", -50.0, sig_a, prfs_a), ("b", -60.0, sig_b, prfs_b))}
    scores = transfer_matrix(subjects, 20)
    assert scores[("b", "b")]["muscle"]["rmse"] == pytest.approx(0.0, abs=1e-12)
    ref = prfs_b[0, 0, 20:]
    expected = abs(-50.0 / -60.0 - 1) * np.sqrt(np.mean(ref**2))
    assert scores[("a", "b")]["muscle"]["rmse"] == pytest.approx(expected)


def test_smooth_in_mask_matches_gaussian_and_ignores_outside():
    img = np.zeros((9, 9, 1))
    img[4, 4, 0] = 1.0
    mask = np.ones((9, 9), bool)
    out = smooth_in_mask(img, mask, sigma=1.0)
    k = np.exp(-np.arange(-2, 3) ** 2 / 2.0)
    expected_centre = 1.0 / np.outer(k, k).sum()
    assert out[4, 4, 0] == pytest.approx(expected_centre)
    mask[0, 0] = False
    assert np.isnan(smooth_in_mask(img, mask, 1.0)[0, 0, 0])


def test_zero_reference():
    x = np.array([1.0, 3.0, 5.0, 7.0])
    assert np.allclose(zero_reference(x, 2), [-1.0, 1.0, 3.0, 5.0])
