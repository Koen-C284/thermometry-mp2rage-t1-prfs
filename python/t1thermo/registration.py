"""
In vivo motion correction with elastix (stage s02; requires `itk-elastix`).

Each dynamic is registered with a B-spline transform to the mean of the first `n_reference` dynamics,
using the echo- and block-averaged magnitude. The transform is applied to every echo and GRE block: the
magnitude is resampled directly and the phase is taken from the resampled real and imaginary parts
(resampling wrapped phase directly would corrupt it near 2*pi wraps). Baseline B1 and T1 maps are aligned to the same reference with an affine transform.
Entry point: python/s02_registration.py
"""

import shutil
from pathlib import Path

import itk
import nibabel as nib
import numpy as np


def _parameters(kind: str, pixel_size_mm: float, grid_px: float = 5.0, samples: int = 2048):
    """elastix settings. B-spline: final grid of `grid_px` pixels and `samples` random samples (paper: 5 px / 2048)."""
    params = itk.ParameterObject.New()
    pmap = params.GetDefaultParameterMap(kind)
    pmap["Metric"] = ["AdvancedMattesMutualInformation"]
    pmap["ImageSampler"] = ["RandomCoordinate"]
    pmap["NumberOfSpatialSamples"] = [str(samples if kind == "bspline" else 2048)]
    pmap["RandomSeed"] = ["121212"]  # fixed sampler seed: reproducible transforms
    if kind == "bspline":
        pmap["NumberOfResolutions"] = ["2"]
        pmap["FinalGridSpacingInPhysicalUnits"] = [str(grid_px * pixel_size_mm)] * 2
        pmap["GridSpacingSchedule"] = ["4.0", "4.0", "2.0", "2.0"]
        pmap["Optimizer"] = ["AdaptiveStochasticGradientDescent"]
        pmap["MaximumNumberOfIterations"] = ["500"]
    else:
        pmap["AutomaticTransformInitialization"] = ["true"]
        pmap["AutomaticTransformInitializationMethod"] = ["GeometricalCenter"]
    params.AddParameterMap(pmap)
    return params


def _itk(img2d: np.ndarray, spacing: float):
    image = itk.GetImageFromArray(np.ascontiguousarray(img2d.T, dtype=np.float32))
    image.SetSpacing([spacing, spacing])
    return image


def _apply(img2d: np.ndarray, transform, spacing: float) -> np.ndarray:
    return itk.GetArrayFromImage(itk.transformix_filter(_itk(img2d, spacing), transform)).T


def register(scan: np.ndarray, spacing: float, n_reference: int = 20, grid_px: float = 5.0, samples: int = 2048):
    """scan: [x, y, slice, echo, dynamic, block, (magnitude, phase)]; returns the registered array and reference image."""
    magnitude = scan[..., 0].mean(axis=(3, 5))[:, :, 0, :]  # [x, y, dynamic]
    fixed = _itk(magnitude[..., :n_reference].mean(axis=-1), spacing)
    params = _parameters("bspline", spacing, grid_px, samples)
    complex_scan = scan[..., 0] * np.exp(1j * scan[..., 1])
    out = np.empty_like(scan)
    for t in range(scan.shape[4]):
        _, transform = itk.elastix_registration_method(fixed, _itk(magnitude[..., t], spacing),
                                                        parameter_object=params, log_to_console=False)
        for z in range(scan.shape[2]):
            for e in range(scan.shape[3]):
                for b in range(scan.shape[5]):
                    # magnitude is resampled directly (complex resampling would lose signal where the phase
                    # varies within a voxel); the phase is taken from the resampled complex signal
                    c = complex_scan[:, :, z, e, t, b]
                    moved = _apply(c.real, transform, spacing) + 1j * _apply(c.imag, transform, spacing)
                    out[:, :, z, e, t, b, 0] = _apply(scan[:, :, z, e, t, b, 0], transform, spacing)
                    out[:, :, z, e, t, b, 1] = np.angle(moved)
    return out, fixed


def align_map(fixed, anatomical: np.ndarray, quantitative: np.ndarray, spacing: float) -> np.ndarray:
    """Affine registration of a map's anatomical image to the reference; the transform is applied to the map."""
    _, transform = itk.elastix_registration_method(fixed, _itk(anatomical, spacing),
                                                    parameter_object=_parameters("affine", spacing), log_to_console=False)
    return _apply(quantitative, transform, spacing)


def main(src: Path, dst: Path, grid_px: float = 5.0, samples: int = 2048) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    img = nib.load(src / "mp2rage_image.nii.gz")
    scan = np.asarray(img.dataobj, dtype=np.float32)
    spacing = float(img.header.get_zooms()[0]) or 1.0
    registered, fixed = register(scan, spacing, grid_px=grid_px, samples=samples)
    nib.save(nib.Nifti1Image(registered, img.affine), dst / "mp2rage_image.nii.gz")
    for name in ("t1_map", "b1_map"):
        anatomical = np.squeeze(np.asarray(nib.load(src / f"{name}_image.nii.gz").dataobj, dtype=np.float32))
        quantitative = nib.load(src / f"{name}.nii.gz")
        aligned = align_map(fixed, anatomical, np.squeeze(np.asarray(quantitative.dataobj, dtype=np.float32)), spacing)
        nib.save(nib.Nifti1Image(aligned.astype(np.float32), quantitative.affine), dst / f"{name}.nii.gz")
    for sidecar in src.glob("*.json"):
        shutil.copy2(sidecar, dst / sidecar.name)
