# Ex vivo porcine muscle phantom dataset

Heated porcine muscle (neck tissue) measured with the 2D multi-echo MP2RAGE sequence on a 3T scanner,
with two fiber-optic temperature probes. No human data and no scanner identifiers are included.

| File | Content |
|---|---|
| `mp2rage_image.nii.gz` | Multi-echo MP2RAGE series `[x, y, slice, echo, dynamic, GRE block (1, 2), (magnitude, phase)]`, 96x96, 11 echoes, 250 dynamics (~4.01 s each) |
| `mp2rage_image.json` | Echo times (ms), voxel size, slice orientation, water-fat shift and the start time of every dynamic (s) |
| `probe_temperatures.txt` | Probe 1 (calibration position) and probe 2 (validation position) in °C, sampled every 0.7 s; the header gives the probe sample at which the first MRI dynamic started |
| `derived/temp_prfs.nii.gz` | PRFS temperature change (°C) from stage s04 (GRE 2) |
| `derived/masks.nii.gz` | Body, muscle and fat masks from stage s04 (last axis) |
| `derived/r2star_map.nii.gz` | Shared R2* map (1/s) from the water-fat separation in stage s05 step 1 |
| `derived/muscle_mask_s05.nii.gz` | Muscle mask from stage s05 step 1 |
| `derived/signal_change_img1.nii.gz`, `..._img2.nii.gz` | Relative signal change dS/S0 of GRE 1 and GRE 2 from stage s05 (no smoothing, no B1/T1 correction, S0 = mean of the first 10 dynamics) |

Each stage can be re-run on its own from these files: s05 reuses `r2star_map` and `muscle_mask_s05` when
present (so it needs only this repository's code), and the Python analysis starts from `derived/`.

Probe positions in the image (0-based array indices `[x, y]`): calibration `[51, 48]`, validation `[46, 40]`.
