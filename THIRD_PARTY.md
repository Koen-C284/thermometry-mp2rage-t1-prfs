# Third-party code and acknowledgements

## Included in this repository

| Code | Location | Licence | Notes |
|---|---|---|---|
| EPG-X (Extended Phase Graphs), S. Malik | `matlab/third_party/EPG-X` | MIT (see `LICENSE` in that folder) | `EPG_GRE2.m`/`EPG_GRE3.m` include extensions by S. Dokumaci (initial state and M0 options), documented in source. |
| HIMM multi-echo PRFS thermometry, body version (M. Kikken), based on the code of M. Poorman and W. Grissom | `matlab/third_party/HIMM` | CC BY-NC 4.0 (see `LICENSE` and `NOTICE.md` in that folder); non-commercial use only | Included with the author's permission; `s04_prfs_thermometry.m` is adapted for MP2RAGE NIfTI I/O (algorithm unchanged). On Linux/macOS, compile `MeanAdjacent/*.c` with MATLAB's `mex` command before running s04 (only Windows `.mexw64` binaries are precompiled). |
| Reduced-dimension nonlinear least-squares T1 fitting (`rdNls` / `rdNlsPr`), J. K. Barral et al. | `matlab/functions/compute_T1.m` | See Barral et al. (2010) | Helper functions `buildNlsStruct`, `rdNls`, and `rdNlsPr` adapted from the Stanford T1 mapping implementation (Barral et al., 2010). |

## Required for some stages, not included

Every stage can be run from the intermediate results shipped in `data/exvivo/` (see the README), so this code
is only needed to recompute those intermediates or to process new data. Place it in one folder and pass that
folder to `setup_paths` (or set the environment variable `T1THERMO_EXTERNAL`).

| Code | Licence | Expected subfolder | Used in | Shipped intermediate that replaces it |
|---|---|---|---|---|
| ISMRM Fat-Water Toolbox (graph-cut fat-water separation, D. Hernando et al.) | No licence published; download from the ISMRM Fat-Water workshop page | `fatwater_toolbox/` | s04 (PRFS), s05 step 1 (R2* map, muscle mask) | `derived/temp_prfs.nii.gz`, `derived/masks.nii.gz`, `derived/r2star_map.nii.gz`, `derived/muscle_mask_s05.nii.gz` |
| Philips PAR/REC reader `loadPARREC` | No licence published | `parrec/` | s01 | `mp2rage_image.nii.gz` + `.json` |

On Windows, keep the external folder at a short path (for example `C:\t1thermo_external`): the toolbox's
compiled graph-cut file (`max_flow_mex.mexw64`) does not load when its full path exceeds 260 characters.

Python dependencies (numpy, scipy, nibabel, matplotlib, itk-elastix) are installed from PyPI; see `requirements.txt`.

## Acknowledgements

The PRFS thermometry in this work uses the HIMM (harmonic initialized model-based multi-echo) method and
implementation by Mathijs Kikken, which builds on the multi-echo thermometry algorithm and code of
Megan Poorman and William Grissom. Please cite:

- Kikken MWI, Steensma BR, van den Berg CAT, Raaijmakers AJE. Multi-echo MR thermometry in the upper leg at
  7 T using near-harmonic 2D reconstruction for initialization. Magn Reson Med. 2023;89(6):2347-2360.
  doi:10.1002/mrm.29591
- Poorman ME, Braškutė I, Bartels LW, Grissom WA. Multi-echo MR thermometry using iterative separation of
  baseline water and fat images. Magn Reson Med. 2019;81(4):2385-2398. doi:10.1002/mrm.27567
- Salomir R, et al. Reference-free PRFS MR-thermometry using near-harmonic 2-D reconstruction of the background
  phase. IEEE Trans Med Imaging. 2012;31(2):287-301. doi:10.1109/TMI.2011.2168421
- Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large field
  inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010;63(1):79-90. doi:10.1002/mrm.22177
- Malik SJ, Teixeira RPAG, Hajnal JV. Extended phase graph formalism for systems with magnetization transfer
  and exchange. Magn Reson Med. 2018;80(2):767-779. (EPG-X)
- Barral JK, Gudmundson E, Stikov N, Etezadi-Amoli M, Stoica P, Nishimura DG. A robust methodology for in vivo
  T1 mapping. Magn Reson Med. 2010;64(4):1057-1067. doi:10.1002/mrm.22497
- Klein S, Staring M, Murphy K, Viergever MA, Pluim JPW. elastix: a toolbox for intensity-based medical image
  registration. IEEE Trans Med Imaging. 2010;29(1):196-205.
