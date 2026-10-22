# KCWI DRP: sparse-matrix cube construction — staged implementation plan

Repo: `~/drp/kcwi/KCWI_DRP` (at `d0b4d7b`, 2026-10-13). Reference implementations
consulted: SCALES DRP (`~/drp/scales/SCALES-DRP`) and PypeIt
(`~/drp/MyPypeIt/PypeIt`, `pypeit/core/datacube.py::subpixellate`).

## 1. Goal

Replace the chain of sequential interpolations that turns a 2D KCWI frame into a
3D cube with **one sparse linear operator `W` (n_voxels × n_pixels) per frame**:

```
cube = W @ sci          var = (W ∘ W) @ var_pix          (optional) Cov = W diag(var_pix) Wᵀ
```

`W` encodes the slice geometry, wavelength solution, air→vacuum,
heliocentric/barycentric correction and DAR, with exact pixel/voxel footprint
overlap weights. Then, optionally, use the same operators to fit dithered frames
jointly.

## 2. Why: current behaviour and its problems

Science path (`kcwidrp/pipelines/kcwi_pipeline.py`, ~L256-271):
`SubtractSky → MakeCube → WavelengthCorrections → CorrectDar → MakeInvsens/FluxCalibrate`.

| Step | File | What it does | Problem |
|---|---|---|---|
| MakeCube | `primitives/MakeCube.py` L19-103, L224-228 | `skimage.transform.warp(order=3)` per slice (24 slices, `multiprocessing.Pool`) using per-slice `asympolynomial` transforms from `geom.pkl`; same warp applied to **σ**, mask (order 0), flags (order 0), noskysub, obj, sky, delta-wave | Bicubic interpolation of σ is not error propagation: ringing, no covariance. Bicubic interpolation of data rings around CRs and bad pixels |
| WavelengthCorrections | `primitives/WavelengthCorrections.py` L46-90 (`heliocentric` L198, `air2vac` L92) | Per-spaxel `interp1d(kind='cubic')` resampling of the cube, applied once for heliocentric and once for air→vac (both on by default, `configs/kcwi.cfg` L15-17) | Two more resamplings along λ |
| CorrectDar | `primitives/CorrectDar.py` L82-, loop ~L242-260 | Pads the cube, then `scipy.ndimage.shift` (spline) on every λ plane of data, **stddev**, mask and flags | Another resampling, and σ is interpolated again |
| GenerateMaps | `primitives/GenerateMaps.py` L38-80 | Makes wavemap/slicemap/posmap with Python loops over every pixel (`for ix…: for iy…`) | Slow, though the transform is already vectorized (`itrf(coords)` accepts (N,2)) |

So data are resampled up to 4 times, the noise is smoothed and correlated each
time, and the reported σ is wrong (no covariance, and interpolating σ is not
propagation). KCWI cube noise is known to be underestimated, and downstream tools
apply empirical covariance corrections for this reason.

What we keep: the geometry calibration (`SolveGeom`: 5 bars/slice, polynomial
order (2,4); `geom.pkl` with `tform`, `invtf`, `xl0`, `xl1`, `xsize`, `ysize`,
`wave0out`, `dwout`, `pxscl`, `slscl`). The matrix approach changes how the
calibration is **applied**, not how accurate it is.

### Approaches compared (decision already made: matrix)

- **KCWI DRP today:** bicubic warp, then λ interpolation ×2, then DAR shift. About 10 s/frame. Problems as above.
- **PypeIt `subpixellate`:** each pixel is split into 5×5 sub-pixels and each slice into 5 sub-slices (125 points/pixel). Points are mapped to voxels, with DAR and dither offsets included, and histogrammed. Only one resampling, but no covariance, about 1-10 min/frame, and the variance weighting looks inconsistent (`subpix_wght**2` vs `**3`, datacube.py ~L3140).
- **SCALES:** explicit sparse matrix R from monochromator PSFlets, Horne/χ² extraction. The concept transfers. The **build code must not be copied**: it uses Python loops with appends, which would take about 8 min at KCWI scale versus about 5 s vectorized (measured).

## 3. Sizes and measured costs (user's Mac, single process)

KCWI 2×2 binned: detector 2048², cube ≈ 24 × ~70 × ~2000 ≈ 3.4M voxels.

| Operation | Measured |
|---|---|
| Build CSR from index arrays, 2×10⁷ nnz (float32) | 0.4 s, 160 MB |
| `W @ d` | 0.02 s |
| `(W∘W) @ var` | 0.03 s |
| `W diag(v) Wᵀ` (full covariance) | 6 s, ~10⁸ nnz |
| `tf.warp(order=3)` one slice, one image | 0.06 s |
| LSQR, 5 stacked frames, 50 iterations | 19 s (random-column worst case) |
| Python loop build (SCALES style), 20k voxels × 100 px | 2.8 s (vectorized: 0.03 s, identical result) |

Benchmark scripts lived in a temporary scratchpad. Re-create them if needed.

## 4. Core technique: building `W` vectorized

**Rule: no Python loops over pixels or voxels.** Produce `rows, cols, vals` as
NumPy arrays and call `scipy.sparse.csr_matrix` once (per slice, then
`sp.vstack`/`sp.block_diag`).

Work **pixel → voxel**: each detector pixel touches only a few voxels.

1. **Pixel edge coordinates** (cache per `geom.pkl`). For each slice, take all
   on-slice pixels (`np.nonzero(slicemap == isl)`) and call
   `invtf[isl](np.column_stack([ix - xl0 ∓ 0.5, iy ∓ 0.5]))` once to get
   `u_lo, u_hi` (along slice, output px) and `lam_lo, lam_hi` (air Å). Slice
   coordinate: `s_lo, s_hi = isl ∓ 0.5`.
2. **Per-frame corrections** (elementwise, under 1 s):
   - `lam = helio(a2v(lam))`. Reuse the math in `WavelengthCorrections.a2v_conversion` and the velocity computation in `heliocentric`.
   - DAR: get the offset vs λ from `CorrectDar.atm_disper` and the same projection onto the slice and cross-slice axes as `CorrectDar` (uses `PXSCL`, `SLSCL`, PA). Add the along-slice part to `u` and the cross-slice part to `s` (in fractional slices).
   - RED channel: the 180° rotation becomes an index relabel of the output grid.
3. **Separable overlap.** Treat each pixel footprint as an axis-aligned box in
   (u, λ, s). The overlap with the regular voxel grid is the product of three 1D
   interval overlaps. If pixel ≤ voxel there are ≤ 2 candidates per axis, so 8
   per pixel:

```python
def overlap_1d(lo, hi, edge0, step, ncand=2):
    i0 = np.floor((lo - edge0) / step).astype(int)
    i  = i0[:, None] + np.arange(ncand)
    a  = np.maximum(lo[:, None], edge0 + i*step)
    b  = np.minimum(hi[:, None], edge0 + (i+1)*step)
    return i, np.clip(b - a, 0, None) / (hi - lo)[:, None]
# w = fu[:,:,None,None]*fl[:,None,:,None]*fs[:,None,None,:]  -> (Npix,2,2,2)
# vox = ravel_multi_index(...); keep w>0 and in-bounds; build csr once
```

   `ncand` must cover the largest footprint/voxel ratio. Compute it per axis
   from the data, not a hard-coded 2. Without DAR, the s-axis overlap is
   trivially 1 for the pixel's own slice.
4. **Mask then normalize.**

```python
O = O @ sp.diags(gpm.ravel().astype(np.float32))      # bad px/CR -> zero column
norm = np.asarray(O.sum(axis=1)).ravel()               # coverage (cf. PypeIt normcube)
W = sp.diags(1/np.where(norm > 0, norm, np.inf)) @ O   # row-normalized
```

   Row normalization gives a weighted mean, i.e. the same units as today's cube
   (interpolated counts per detector pixel). That keeps `MakeInvsens` and
   `FluxCalibrate` consistent. Keep `norm` as a coverage cube.
5. **Memory.** About 4M px × 8 gives about 3×10⁷ temporaries. Build per slice to
   bound peak memory. Use float32 vals and int32 indices.

If the box approximation is too coarse (check the footprint tilt from the
transform Jacobian), split each pixel into 2×2 parts (still vectorized) or clip
the true quadrilateral against the grid.

**Correction to an earlier claim:** `tf.warp(order=3)` uses spline prefiltering,
which is not a finite-support kernel. A bicubic warp therefore **cannot** be
reproduced exactly as a sparse matrix. For bit-level validation of the matrix
machinery, compare against `tf.warp(order=1)` (bilinear), which has finite
support, using a bilinear-weight variant of `W`.

## 5. Stages

Each stage is mergeable on its own, sits behind a config switch until validated,
and keeps the old path for comparison.

### Stage 0: Baseline and harness

- Pick a reference dataset: one BL or BM night with a standard star, a
  science frame with a point source, arcs, contbars and flats. Include a
  nod-and-shuffle frame if available, and RED-channel (KCRM) data.
- Run the current pipeline and archive `icube`, `icubew`, `icubed`,
  `invsens`, `icubes`, and per-primitive timings.
- Write `kcwidrp/tests/test_cube_matrix.py` (fast, synthetic) and a separate
  regression script for the real data (not in CI).
- Metrics:
  - arc-line centroid and width in the cube vs `geom['avwvsig']`
  - contbar positions along the slice
  - standard-star total flux and invsens curve
  - σ calibration: in sky-only regions compare the reported σ with the empirical
    scatter, including binned/aperture sums (the covariance test)
  - runtime

### Stage 1: Vectorized coordinate maps (quick win, independent)

- Rewrite the `GenerateMaps._perform` loops (L64-80) to call `itrf` once per slice
  on all pixels. Output must match the current maps exactly (same pixel-centre
  convention, same `0 <= x <= xsize` test).
- Add a helper (e.g. `kcwidrp/core/cube_matrix.py`) that returns per-pixel
  **edge** coordinates (`u_lo/hi`, `lam_lo/hi`, slice id, flat pixel index)
  from `geom`, cached next to `*_geom.pkl` (e.g. `*_pixcoords.npz`).
- Acceptance: wavemap/slicemap/posmap identical to the current output.
  `MakeMasterSky` (which reads them) is unchanged. Large speedup.

### Stage 2: `W` for geometry only, replacing the warp in MakeCube

- `cube_matrix.build_W(geom, pixcoords, gpm, out_grid, mode='overlap'|'bilinear')`.
- Config: `cube_method = "warp" | "matrix"` in `kcwi.cfg` (default `warp`
  until validated).
- In `MakeCube`, if `matrix`: build W, then compute each output with one matvec:
  - `icube` data: `W @ sci`
  - uncertainty: `sqrt((W∘W) @ σ²)`
  - noskysub, obj (`_objf`), sky (`_skyf`), delta-wave map: `W @ …`
- **Mask and flags:** flags are bitmasks, so don't average them. For each bit
  plane `b`: `((O_bool @ bitplane_b) > 0)`, then OR the planes together, where
  `O_bool` is the structure of the unmasked overlap matrix. The mask cube is
  `coverage == 0` (or below a threshold), plus the flag-derived mask.
- Keep the 180° RED rotation (it can stay as `np.rot90` on the result, or be
  folded into the index mapping) and all the WCS/header code unchanged.
- Validate the machinery first with `mode='bilinear'` against
  `tf.warp(order=1)`: match to ~1e-6. Then switch to `overlap`.
- Acceptance (overlap mode vs current warp):
  - arc centroids agree within a small fraction of a pixel
  - standard-star flux agrees within ~1%
  - in sky regions the ratio of reported σ to empirical scatter moves toward 1
  - MakeCube time ≤ current

### Stage 3: Fold the wavelength corrections into `W`

- Compute the radial-velocity correction and the air→vac mapping *before*
  building the cube, and apply them to `lam_lo/hi`. The output grid is defined
  in the corrected frame (vacuum, heliocentric).
- In `WavelengthCorrections`: if the cube was built by the matrix method and
  already includes the corrections (set a header flag, e.g. `WAVCORR=T`), skip
  the resampling but still write `VCORR`, `VCORRTYP`, `CTYPE3='WAVE'`, etc., and
  still produce the `icubew` output for pipeline compatibility.
- Pipeline order: the correction parameters only need header information
  (time, RA/Dec, Keck location), so compute them in a small helper called from
  MakeCube. No event-table reordering is needed.
- Acceptance: sky-line and arc-line positions in `icubew` match the current
  output. Line profiles are narrower or equal, since there are fewer cubic
  resamplings.

### Stage 4: Fold DAR into `W`

- Compute the per-λ (Δalong-slice, Δcross-slice) offsets from the same
  `atm_disper` plus projection code as `CorrectDar`. Add Δu to `u` and Δs to `s`.
  Use the padded output grid that CorrectDar uses (`padding_x/y` per grating,
  CorrectDar ~L100-111).
- Cross-slice DAR moves flux by a fraction of a slice. With one frame the
  overlap weights split a pixel between two adjacent slice voxels. This is
  equivalent in resolution to today's `ndimage.shift`, but done once and with
  correct variance. Document that true cross-slice resolution needs Stage 6.
- In `CorrectDar`: if the cube already includes DAR (`DARCOR=T`), skip the
  shifts and keep the header updates, padding flags (value 128) and `icubed`
  output.
- Acceptance: standard-star centroid vs λ is flat, at least as well as the
  current path. Standard-star flux agrees. **Re-derive invsens with the new path**
  and compare it to the old invsens.

### Stage 5: Covariance products (optional output)

- Option A (cheap): for each voxel, store the correlation with its ±1 neighbours
  along each axis (6 extra planes, or a compact summary). These come from
  `W diag(var) Wᵀ`, computed per slice.
- Option B: store `W` itself (`.npz`, about 160 MB) with the cube, so users can
  compute exact errors for any aperture `a`: `var(aᵀcube) = (Wᵀa)ᵀ diag(var) (Wᵀa)`.
- Option C: compute a binned-aperture noise-scaling curve, which is what CWITools
  users fit empirically today, and write it to the header.
- Decide with maintainers. A is probably the default and B opt-in.

### Stage 6 (research / optional): Joint inversion of dithered frames

- Separate script or primitive (in the style of `kcwi_combine`), not in the
  per-frame pipeline.
- Model: `d_k = A_k f + noise`, where `A_k` (pixels × voxels) maps a common
  sky grid (finer than the slice width) to frame k, including the per-frame
  DAR, offsets and position angle. Solve weighted least squares with LSQR/CG and
  Tikhonov regularization. Optionally add shared sky B-spline terms (replacing
  MakeMasterSky for those frames).
- Gains: sub-slice spatial sampling from dithers or multiple position angles,
  and a single resampling for the combined cube.
- Cost: about 1-10 min for 5 frames, 100-300 iterations. Avoid bounded/NNLS
  (`lsq_linear`) as the default: it takes minutes to hours.
- Open issues: how to choose the regularization, how to report a deconvolved vs
  re-convolved cube, and how to obtain the posterior covariance (approximate).

## 6. Expected downstream effects

- Cube values: same units (row-normalized). Small changes in effective spatial
  and spectral resolution, because there are fewer resamplings and overlap
  weights replace bicubic (no ringing, slightly different smoothing kernel).
- Reported σ grows to its true value, and aperture-summed errors become
  correct, especially once covariance is available.
- `invsens` and `FluxCalibrate`: must be re-derived using the same cube method
  as the science frames. Never mix methods.
- `icubew` and `icubed` keep their names. Their content comes from the single
  resampling, and header flags mark which method was used.
- Runtime: MakeCube + WavelengthCorrections + CorrectDar together go from ~10+ s to
  ~1-5 s per frame. GenerateMaps goes from seconds-to-minutes to under 1 s.
- MakeCube runs for four event-table entries: `arc_make_cube`,
  `flat_make_cube`, `object_make_cube` and `nandshuff_make_cube`
  (kcwi_pipeline.py L158, L207, L259, L308). All four should use the same
  method, which also covers the arc `dcube`.

## 7. Open questions to settle before or while implementing

1. Is KCRM's detector and binning handled identically by `geom` (xsize/ysize,
   rotation)? Test on RED data early.
2. Nod-and-shuffle: confirm the obj/sky frame handling and the shuffled-row
   limits (`wavensgood*`) map cleanly onto the coverage/mask logic.
3. Does the header WCS (CRPIX/CD written in MakeCube ~L420-470) stay correct
   with a padded, DAR-folded grid? CorrectDar currently adds the padding to
   CRPIX1/2 (CorrectDar.py L114-117). Replicate that.
4. Choice of flag-combination semantics (OR of contributing pixels vs
   weight-thresholded).
5. Upstream acceptance: KCWI_DRP is a Keck project. Propose Stages 1-2 as a
   first PR behind the config switch.

## 8. Related notes from the review (not in scope)

- PypeIt `subpixellate` variance weighting looks inconsistent with its flux
  weighting (`subpix_wght**2` per sub-pixel, giving n³ instead of n⁴ per pixel).
  This is worth reporting upstream after an empirical check.
- SCALES DRP has a separate drafted issue: vectorize its R-matrix build and
  stop `crop_sparse_vals` from mutating `ims_cal` in place.
