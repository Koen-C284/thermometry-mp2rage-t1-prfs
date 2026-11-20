# Simultaneous T1-based and PRFS thermometry with a multi-echo MP2RAGE sequence

Code and ex vivo data for the paper *"Simultaneous T1-based and PRFS Thermometry with a Multi-Echo MP2RAGE
sequence"*. The T1-weighted magnitude change of both MP2RAGE readout blocks is converted to temperature with a
zero-intercept calibration, `dT = alpha * dS/S0`, derived from concurrent PRFS thermometry (in vivo) or from a
fiber-optic probe (ex vivo).

## Quick start: ex vivo results

```bash
pip install -r requirements.txt
python python/s06_exvivo.py
```

This reproduces the ex vivo analysis (Figs. 5, 6, 9A and Supplementary Fig. 2) from the shipped dataset and writes
figures, calibrated temperature maps and `exvivo_results.json` to `results/exvivo/`.

To also re-run the MATLAB stages on the ex vivo data:

```matlab
cd matlab
run_exvivo                                   % s05 from the shipped intermediates (this repository's code only)
run_exvivo('path/to/external_code', true)    % also recompute s04 (PRFS) and s05 step 1 (needs the Fat-Water Toolbox)
```

Unit tests: `python -m pytest python/tests`.

## Repository layout

```
config/            sequence_parameters.json (EPG timing), exvivo.json, invivo.json (analysis settings)
data/exvivo/       ex vivo dataset: inputs and the intermediate results of every stage (see its README)
matlab/            s01_compute_maps.m, s03_epg_dictionary.m, s04_prfs.m, s05_t1_thermometry.m, run_exvivo.m
  functions/       helpers of the MATLAB stages (T1/B1 maps, EPG dictionary)
  third_party/     EPG-X (MIT) and HIMM PRFS thermometry (CC BY-NC 4.0)
python/            s02_registration.py, s06_exvivo.py, s06_invivo.py
  t1thermo/        library used by the Python stages
  tests/
results/           created when the stages are run
```

## Pipeline

| Stage | Entry point | What it does | Needs external code | Ex vivo |
|---|---|---|---|---|
| s01 | `matlab/s01_compute_maps.m` | PAR/REC to NIfTI; Look-Locker T1 map and AFI B1 map | PAR/REC reader | output shipped (`mp2rage_image`) |
| s02 | `python/s02_registration.py` | In vivo motion correction (elastix B-spline, complex-valued resampling) | – | not needed |
| s03 | `matlab/s03_epg_dictionary.m` | EPG dictionary and B1/T1-dependent relative T1 sensitivity map | – | not used |
| s04 | `matlab/s04_prfs.m` | Fat-referenced, drift-corrected multi-echo PRFS temperature (bundled HIMM code) | Fat-Water Toolbox | output shipped (`temp_prfs`, `masks`) |
| s05 | `matlab/s05_t1_thermometry.m` | Step 1: shared R2* map and muscle mask (water-fat separation); steps 2–4: S0 per dynamic, smoothing, dS/S0, sensitivity correction | step 1 only: Fat-Water Toolbox | step 1 output shipped; steps 2–4 run |
| s06–s07 | `python/s06_exvivo.py`, `python/s06_invivo.py` | Calibration, validation, transferability and figures | – | runs |

Every stage can be run on its own from the shipped intermediates with only this repository's code. The external
code (see [THIRD_PARTY.md](THIRD_PARTY.md)) is needed only to recompute s04 and s05 step 1 or to convert new
PAR/REC data. Pass its folder to `setup_paths` / `run_exvivo`, or set the environment variable `T1THERMO_EXTERNAL`.

Sequence timing for the EPG simulations is hard-coded in `config/sequence_parameters.json`, so no scanner log
files are needed. Stage settings (baseline length, smoothing, sensitivity correction, voxel selection) are in
`config/exvivo.json` and `config/invivo.json`.

## Methods implemented

- **S0 (s05):** a voxelwise R2* map is estimated once from the baseline-averaged GRE 2 data with graph-cut
  water-fat separation; with R2* fixed, S0 = Σ S·exp(−R2*·TE) / Σ exp(−2R2*·TE) for each dynamic and block.
  (`"s0_fit": "wls"` instead fits S0 and R2* per dynamic with weighted log-linear least squares.)
- **Ex vivo:** the first 6 dynamics are discarded (approach to steady state) and the next 6 form the common
  baseline of the probe, PRFS and T1-based series. The probe trace is aligned on the acquisition time of every
  dynamic. alpha is fitted at the probe-1 voxel against probe 1, per GRE block. Probe 2 is used only for
  validation. Figs. 5 and 6 show the mean of the GRE 1 and GRE 2 temperatures. Precision is the SD of the
  residuals of a quadratic fit in time, without and with 5x5 Gaussian filtering (sigma = 1 pixel).
- **In vivo:** dS/S0 (5x5 Gaussian, sigma = 1 pixel, divided by the EPG relative sensitivity map) is calibrated
  against PRFS in muscle voxels with a PRFS rise of at least 0.5 °C and |r| of at least 0.9. alpha is the
  R²-weighted mean of the voxel slopes. The reference subject's alpha is applied to the other subjects, and every
  subject is also calibrated independently to quantify transferability (RMSE against PRFS after the 75-dynamic
  baseline).

## Data

`data/exvivo/` contains the ex vivo porcine muscle phantom dataset (see its README). The in vivo volunteer data
are not shared. `config/invivo.json` documents the expected layout, so the in vivo analysis can be run on
equivalent data.

## Licence

Code: MIT ([LICENSE](LICENSE)); `matlab/third_party/` keeps its own licences (EPG-X: MIT; HIMM: CC BY-NC 4.0,
non-commercial use only). Data in `data/`: CC BY 4.0 ([data/LICENSE](data/LICENSE)).

## Citation

Please cite the paper when using this code or data (reference to be added upon publication).
