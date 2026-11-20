# HIMM multi-echo PRFS thermometry (body version)

This folder contains the HIMM (harmonic initialized model-based multi-echo) PRFS thermometry code
used for the paper, including the driver script `s04_prfs_thermometry.m` and the water-fat helper functions
used in stage s05 (`perform_water_fat_separation.m`, `get_masks.m`).

- Creator: Mathijs Kikken (University Medical Center Utrecht), included with his permission.
- Based on the multi-echo water/fat separated thermometry code of Megan Poorman and William Grissom
  (Vanderbilt University Institute of Imaging Science):
  https://github.com/poormanme/waterFatSeparated_MRThermometry
- Public release of HIMM: https://github.com/MathijsKikken/MR_Thermometry_HIMM

Licence: Creative Commons Attribution-NonCommercial 4.0 International (see `LICENSE`), the licence of the public
HIMM and Poorman/Grissom repositories. This folder is **not** covered by the MIT licence of the rest of the
repository; commercial use is not permitted.

Modifications with respect to the upstream HIMM repository: `s04_prfs_thermometry.m` was adapted for the
MP2RAGE NIfTI inputs/outputs of this pipeline, and commented-out local debug paths and unused `addpath`
statements were removed; the core HIMM algorithm is unchanged. Older backup copies and editor autosave files
were omitted.

Notes for running `s04_prfs_thermometry.m`: set `scan_data_dir` to a folder containing
`02_interim/03_registered/mp2rage_image.nii.gz` and `.json`; results are written to
`03_processed/prfs_thermometry/`. The script also saves extra copies under `\Optimization_multiref_poly3\`
at the root of the current drive. It requires the ISMRM Fat-Water Toolbox (see `THIRD_PARTY.md`).
On Linux or macOS, compile `MeanAdjacent/MeanAdjacentWithMask.c` and
`MeanAdjacent/MeanAdjacentWithMaskIncChange.c` using MATLAB's `mex` command (only Windows `.mexw64` binaries
are precompiled).

Please cite:
- Kikken MWI, Steensma BR, van den Berg CAT, Raaijmakers AJE. Multi-echo MR thermometry in the upper leg at 7 T
  using near-harmonic 2D reconstruction for initialization. Magn Reson Med. 2023;89(6):2347-2360. doi:10.1002/mrm.29591
- Poorman ME, Braškutė I, Bartels LW, Grissom WA. Multi-echo MR thermometry using iterative separation of baseline
  water and fat images. Magn Reson Med. 2019;81(4):2385-2398. doi:10.1002/mrm.27567
- Salomir R, et al. Reference-free PRFS MR-thermometry using near-harmonic 2-D reconstruction of the background
  phase. IEEE Trans Med Imaging. 2012;31(2):287-301. doi:10.1109/TMI.2011.2168421
