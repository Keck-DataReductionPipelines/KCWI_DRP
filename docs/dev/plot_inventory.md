# Plot and prompt inventory (task S.5 / bokeh-plan 0.2)

Taken at commit `a360685` (`git rev-parse --short HEAD`, branch `upgrade/pr0-tests`). **All line numbers in this file refer to that commit** and will move after Track A. The primitives and scripts inspected were unmodified in the working tree at that time.

Conventions:

- `pl` is `self.config.instrument.plot_level`. "pl>=1" means `plot_level >= 1`.
- `STD` means `figure(title=..., x_axis_label=..., y_axis_label=..., plot_width=<plot_width>, plot_height=<plot_height>)`. Any other `figure()` keyword is spelled out.
- `plab` is the string from `plotlabel(self.action.args)` (`"[ Img # %d (%s) Slicer: %s, Grating: %s, CWave: %d[, Filter: %s ] "`), prepended to most titles. `kcwi_file_primitives.plotlabel`, `self.action.args.plotlabel` and `self.action.args.stdlabel` are used as given in the source.
- Every figure is pushed with `bokeh_plot(p, self.context.bokeh_session)` straight after it is built, so it is the figure on screen until the next `bokeh_plot` call. "Live-only" means no `save_plot` is applied to it.
- Figure IDs `F01`-`F48` are in file order (alphabetical by primitive, then by line); prompt IDs `P01`-`P51` and `S01`-`S08` are likewise.
- Sections: 1 figure table, 2 feature summary, 3 prompt table, 4 sleep list, 5 counts and verification, 6 easy-to-miss items.

## 1. Figure table (48 rows, 16 primitives)

`StartBokeh.py:41` (`p = figure()`, the placeholder root of the Bokeh session) is not a row.

| ID | File | Line | Var | Title string (as in source) | Figure guard | Saved? | `save_plot` file name and its guard | Bokeh features |
|---|---|---|---|---|---|---|---|---|
| F01 | ArcOffsets.py | 79 | `p` | `plab + "BAR OFFSET for Arc: %d Slice: %d = %d" % (arc_number, int(arc_number/5), offset)` | `do_plot` (pl>=2) and `arc_number == next_bar_to_plot`; inside the per-arc loop | live-only | - | STD; `line(color, legend_label)` x2 |
| F02 | ArcOffsets.py | 108 | `p` | `plab + "BAR OFFSETS "` | pl>=1 | saved (L125) | `pfname + '_baroffs.png'`, `pfname = "arc_%05d_%s_%s_%s_ta%03d"`; pl>=1 (same block; save is after the prompt/sleep) | STD; `diamond(size)`; `xgrid.grid_line_color`; `oplot_slices`; `set_plot_lims` |
| F03 | ExtractArcs.py | 140 | `p` | `plab + "ARC # %d" % len(arcs)` | `do_plot` (pl>=3); branch `NBARS < 60` | live-only | - | STD; `line(legend_label, color)` x3 |
| F04 | ExtractArcs.py | 207 | `p` | `plab + "ARC # %d" % len(arcs)` | `do_plot` (pl>=3); branch `NBARS >= 60` | live-only | - | STD; `line(legend_label, color)` x2 |
| F05 | FindBars.py | 85 | `p` | `plab + "BARS TRACE at = %d, TRY %d" % (middle_y_row, n_tries)` | pl>=3; inside the peak-finding retry loop | live-only | - | STD; `line(color, legend_label)`; `scatter(marker='x', color, legend_label)`; `line(color, line_dash='dashed')`; `legend.location` |
| F06 | FindBars.py | 124 | `p` | `plab + "BARS MID TRACE Thresh = %.2f" % bar_thresh` | pl>=1; only when `peaks_found == NBARS` | live-only | - | same as F05 (`legend_label="MidTrace"`) |
| F07 | FindBars.py | 154 | `p` | `plab + "BAR %d CENTRD = %.2f" % (ip, xc)` | `do_inter` (pl>=3, cleared by answering q); per bar | live-only | - | STD; `line(color, legend_label)`; `circle(color, legend_label)`; `line(color, legend_label)` |
| F08 | FitCenter.py | 324 | `p` | `plab + "CENTRAL DISPERSION FIT for Bar: %d Slice: %d" % (b, int(b / 5))` | `do_inter` (pl>=2) and `ir == next_bar_to_plot`; per bar. Note `p` was a spawn `Pool` at L302-304 | live-only | - | STD (kw order differs); `scatter(color, legend_label)` **no `marker`**; `line(color, legend_label)` x3 |
| F09 | FitCenter.py | 353 | `p` | `plab + "CENTRAL VALUES"` (y label "Central Wavelength (A)") | pl>=1 | saved (L378) | `pfname + '_cwaves.png'`, `pfname = "arc_%05d_%s_%s_%s"`; pl>=1 (same block; after prompt/sleep) | STD; `scatter(marker='x', legend_label)`; `line(color, legend_label)`; `xgrid.grid_line_color`; `oplot_slices`; `legend.location`; `set_plot_lims` |
| F10 | FitCenter.py | 380 | `p` | `plab + "CENTRAL VALUES"` (y label "Central Dispersion (A)") | pl>=1 | saved (L401) | `pfname + '_cdisp.png'`; pl>=1 | as F09 |
| F11 | GetAtlasLines.py | 520 | `p` | `plotlabel(self.action.args) + "ATLAS LINES Ngood = %d, Nrej = %d" % (len(refws), nrej)` | pl>=1 | saved (L565) | `atfnam + ".png"`, `atfnam = "arc_%05d_%s_%s_%s_atlines"`; pl>=1 (same block; after prompt/sleep) | STD; `line(legend_label, color)`; `diamond(legend_label, color, size)` up to 6 calls (size 8, final "Kept" size 16); `x_range = Range1d(a, b)`; `y_range = Range1d(-0.04, 1.04)` |
| F12 | MakeCube.py | 269 | `p` | `self.action.args.plotlabel + "WARPED Slice %d" % isl` | pl>=3; loop over 24 slices. Note `p` was a `Pool` at L244-246 | live-only | - | `figure(tooltips=[("x","$x"),("y","$y"),("value","@image")], title, x_axis_label, y_axis_label, plot_width, plot_height)`; `x_range.range_padding = y_range.range_padding = 0`; `image([warped], x, y, dw, dh, palette="Spectral11", level="image")` |
| F13 | MakeInvsens.py | 186 | `p` | `self.action.args.stdlabel + ' Std Slice %d (DAR padded)' % i` | `do_plots` (pl>=3); per slice (loop) | live-only | - | STD; `scatter(px, py, marker='x')` |
| F14 | MakeInvsens.py | 348 | `p` | `self.action.args.stdlabel + ' Obs Spec'` | pl>1 (i.e. >=2); `while not done` loop | live-only | - | `figure(tooltips=[("x","@x{0,0.0}"),("y","@y{0,0.0}")], STD)`; `ColumnDataSource(data=dict(x=w, y=obsspec))` (L345); `line('x','y', line_color, source=source)`; `line(line_color, legend_label)`; `set_plot_lims` |
| F15 | MakeInvsens.py | 446 | `p` | `self.action.args.stdlabel + ' Obs Spec'` | pl>=2; `while not done` loop | live-only | - | `figure(tooltips=[("x","@x{0.0}"),("y","@y{0.0}")], STD)` (no `source` is ever supplied, so the `@x`/`@y` fields have no data source); `line(w, obsspec, line_color)`; `line(line_color, legend_label)`; `line(line_color, line_dash='dashed')` per mask; `set_plot_lims` |
| F16 | MakeInvsens.py | 564 | `peff` | `self.action.args.stdlabel + ' Efficiency'` | pl>=1; `while not done` loop | saved (L699) | `pfname + '_eff.png'`, `pfname = "std_%05d_%s_%s_%s_%d"`; pl>=1 (after the loop) | STD; `line(line_color, legend_label)` x2; `scatter(marker='x', legend_label)`; `line(line_color, legend_label)`; `line(line_color)`; `line(line_color, line_dash='dashed'/'dashdot')`; `set_plot_lims` |
| F17 | MakeInvsens.py | 590 | `pivs` | `self.action.args.stdlabel + ' Inverse sensitivity'` | pl>=1; loop | saved (L700) | `pfname + '_invsens.png'`; pl>=1 | `figure(..., y_axis_type='log', STD)`; `line(line_color, legend_label)` x3; `scatter(marker='x', legend_label)`; `line(line_color)`; `set_plot_lims` |
| F18 | MakeInvsens.py | 611 | `pcal` | `self.action.args.stdlabel + ' Calibrated'` | pl>=1; loop | saved (L702) | `pfname + '_cal.png'`; pl>=1 | STD; `line(line_color, legend_label)` x3; `line(line_color)`; `set_plot_lims` |
| F19 | MakeInvsens.py | 630 | `prsd` | `self.action.args.stdlabel + ' Residuals = %.1f +- %.1f (%%)' % (rsd_mean, rsd_stdv)` | pl>=1; loop | saved (L701) | `pfname + '_resid.png'`; pl>=1 | STD; `line(line_color, legend_label)`; `scatter(marker='x', legend_label)` only if `len(not_used) > 0`; `line`; `line(line_color, line_dash='dashed')` x2; `set_plot_lims` |
| F20 | MakeMasterBias.py | 139 | `p` | `plabel + 'BIAS NOISE amp %d = %.3f' % (ia, bias_rn)` (`plabel` is built locally from header keywords) | pl>=1 | saved (L157) | `biasfnam + ".png"`, `biasfnam = "bias_%05d_%s_amp%d_rdnoise"`; pl>=1 (after prompt/sleep) | STD; `quad(top, bottom, left, right, fill_color, line_color, alpha)`; `line(line_color, line_width, alpha, legend_label)`; `line(color, legend_label)`; `line(color)`; `y_range.start = 0` |
| F21 | MakeMasterFlat.py | 260 | `p` | `plab + ' WAVE SLOPE FIT'` | pl>=1; inside `if internal:` (SFLAT) | live-only | - | STD; `circle(legend_label)`; `line(line_color, line_width, legend_label)` |
| F22 | MakeMasterFlat.py | 309 | `p` | `plab + ' Vignetting'` | pl>=1; inside `if internal:` | live-only | - | STD; `circle(legend_label)`; `line(line_color, legend_label)` x2; `line(line_color)` x7; `y_range = Range1d(ycmin, ycmax)` |
| F23 | MakeMasterFlat.py | 366 | `p` | `plab + ' Buffer Region'` | pl>=1; inside `if internal:` | live-only | - | STD; `circle()`; `line(line_color)` |
| F24 | MakeMasterFlat.py | 441 | `p` | `plab + ' Ledge'` | **no `plot_level` guard**; inside `if 'BM' in grating and not twiflat:` and `wavegood0 <= ledge_wave <= wavegood1`. Its `bokeh_plot` (L450) and sleep branch (L454) also run at pl 0 | live-only | - | STD; `circle(fill_color)`; `line()`; `line(color)` |
| F25 | MakeMasterFlat.py | 455 | `p` | `plab + ' Deriv'` | **no `plot_level` guard** (same BM block); `bokeh_plot` at L466 also unguarded | live-only | - | STD; `circle()`; `line(color)` x1 plus one per peak |
| F26 | MakeMasterFlat.py | 483 | `p` | `plab + ' Peak of ledge'` | pl>=3 (BM block) | live-only | - | STD; `circle(legend_label)`; `line(line_color, legend_label)` |
| F27 | MakeMasterFlat.py | 518 | `p` | `plab + ' BM Ledge Region'` | pl>=1 (BM block) | live-only | - | STD; `circle(fill_color, legend_label)` x2; `line(line_color, legend_label)`; `line(line_color)` x6; `line(line_color, legend_label)`; `y_range = Range1d(ylmin, ylmax)`; `legend.location = 'top_left'` |
| F28 | MakeMasterFlat.py | 643 | `p` | `plab + ' Blue/Red fits'` | pl>=1 (block starts L620) | saved (L661) | `rbfnam + '.png'`, `rbfnam = "redblue_%05d_%s_%s_%s"`; pl>=1 | `figure(title, x_axis_label, y_axis_label, y_range=ylim, plot_width, plot_height)` with `ylim` a 2-list; `line(line_color, legend_label)` x3; `circle(size=1, line_alpha=0., fill_color)` x3 |
| F29 | MakeMasterFlat.py | 713 | `p` | `plab + ' Blue fits'` | pl>=1 and `nqb > 1` | live-only | - | STD; `line(line_color, legend_label)` x2; `circle(fill_color)` x2 |
| F30 | MakeMasterFlat.py | 731 | `p` | `plab + ' Blue ratios'` | pl>=1 and `nqb > 1` | live-only | - | as F29 (y label "Ratio") |
| F31 | MakeMasterFlat.py | 752 | `p` | `plab + ' Red fits'` | pl>=1 and `nqr > 1` | live-only | - | as F29 |
| F32 | MakeMasterFlat.py | 770 | `p` | `plab + ' Red ratios'` | pl>=1 and `nqr > 1` | live-only | - | as F29 |
| F33 | MakeMasterFlat.py | 871 | `p` | `plab + ' Master Illumination'` | pl>=1 (block starts L856) | saved (L891) | `fltfnam + ".png"`, `fltfnam = "flat_%05d_%s_%s_%s"`; pl>=1 | `figure(..., y_range=yran, ...)` with `yran` a 2-list; `circle(size=1, line_alpha=0., fill_color, legend_label)`; `line(line_color, legend_label, line_width=2)`; `line(line_color, legend_label)`; `line(line_color)`; `legend.location = "top_left"` |
| F34 | MakeMasterSky.py | 331 | `p` | `self.action.args.plotlabel + ' Std max sl %d' % std_sl_max` | pl>=1; inside `if self.action.args.stdname is not None:` | live-only | - | STD; `circle(size=1, line_alpha=0., fill_color='purple', legend_label)`; `line(color)` x4 |
| F35 | MakeMasterSky.py | 413 | `p` | `self.action.args.plotlabel + ' Cont source max sl %d' % con_sl_max` | pl>=1; inside `elif self.action.args.contsky:` and `cont_source_pos is None` | live-only | - | STD; `circle(size=1, line_alpha=0., fill_color, legend_label)`; `line(color)` x6 |
| F36 | MakeMasterSky.py | 592 | `p` | `self.action.args.plotlabel + ' Master Sky'` | pl>=1 | saved (L610) | `skyfnam + ".png"`, `skyfnam = "sky_%05d_%s_%s_%s"`; pl>=1 (after prompt/sleep) | STD; `circle(size=1, line_alpha=0., fill_color, legend_label)`; `line(line_color, legend_label)`; `line(line_color)` x2; `y_range.start = yrng[0]`; `y_range.end = yrng[1]` |
| F37 | ReadAtlas.py | 178 | `p` | `plab + "ATLAS OFFSET = %d px" % offset_pix` | pl>=1 | live-only | - | STD; `line(legend_label)`; `line(color, legend_label)`; `line(color, legend_label)` only if `atoff != 0` |
| F38 | ReadAtlas.py | 202 | `p` | `plab + "ATLAS OFFSET = %.1f Ang (%d px)" % (offset_wav, offset_pix)` | pl>=1; rebuilt on every pass of `while q:` | saved (L240), the last pass only | `pfname + '_atoff.png'`, `pfname = "arc_%05d_%s_%s_%s"`; pl>=1 (after the loop) | STD; `line(legend_label)`; `line(color, legend_label)`; `x_range = Range1d(minwav, maxwav)` (positional); `line(color, legend_label, line_dash="dashdot")` |
| F39 | SolveAIT.py | 155 | `p` | `plotlabel(self.action.args) + "ARC # %d" % ib` | `do_plot` (pl>=1); per bar | saved (L170) | `pfname + '_bar%d.png' % ib`, `pfname = "arc_%05d_%s_%s_%s_tf%02d"`; `do_plot` (pl>=1); **save is before the prompt** | STD; `line(legend_label, color)` x4 (Arc, Bar, Good WLs, All WLs); `line(color)` x2 (the second line of each limit pair); `set_plot_lims` |
| F40 | SolveArcs.py | 226 | `p` | `plab + "ATLAS/ARC LINE FITS" + ptitle` | `do_inter` (pl>=3) and `ib == next_bar_to_plot`; per line, inside `try` | live-only | - | STD; `line(color, legend_label)`; `circle(color, legend_label)`; `x_range = Range1d(start=, end=)`; `extra_x_ranges = {"pix": Range1d(start=, end=)}`; `add_layout(LinearAxis(x_range_name="pix", axis_label="CCD Y pix"), 'above')`; `line(..., x_range_name="pix")`; `circle(..., x_range_name="pix")`; `line_dash='dashed'/'dashdot'`; `y_range.start = 0` |
| F41 | SolveArcs.py | 389 | `p` | `plab + "RESIDUALS" + ptitle` | `master_inter` (pl>=2) and `ib == next_bar_to_plot`; per bar | live-only | - | STD; `diamond(color, legend_label, size=8/6)` x3; `line(color, line_dash='dotted'/'dashdot')` x3; `line(legend_label, color, line_dash)`; `set_plot_lims` |
| F42 | SolveArcs.py | 415 | `p` | `plab + "ATLAS/ARC FIT" + ptitle` | same as F41 | live-only | - | STD; `line(color, legend_label)` x2; `line(color, line_dash, legend_label)`; `diamond(color, size)`; `diamond(legend_label, color, size)` x2 (the second and third in `if` blocks, `color` + `legend_label` + `size`) |
| F43 | SolveArcs.py | 480 | `p` | `ptitle` (`plab + "COEF %d VALUES <C%d> = ..."`) | pl>=1; per coefficient | saved (L497) | `pfname + '_coef%d.png' % cn`, `pfname = "arc_%05d_%s_%s_%s_tf%02d"`; pl>=1 (after prompt/sleep) | STD; `diamond(size)`; `xgrid.grid_line_color`; `oplot_slices`; `set_plot_lims` |
| F44 | SolveArcs.py | 505 | `p` | `ptitle` (`plab + "FIT STATS <Nlns> = %.1f +- %.1f"`) | **built at every `plot_level`**; its `bokeh_plot` is pl>=1 (L530) | saved (L537) | `pfname + '_nlines.png'`; **no guard: saved at pl 0 too** | STD; `diamond(size)`; `line(color)`; `line(color, line_dash='dashed')` x2; `xgrid.grid_line_color`; `oplot_slices`; `set_plot_lims` |
| F45 | SolveArcs.py | 549 | `p` | `ptitle` (`plab + "FIT STATS <RMS> = %.3f +- %.3f"`) | built at every `plot_level`; `bokeh_plot` is pl>=1 (L570) | saved (L578) | `pfname + '_resid.png'`; **no guard: saved at pl 0 too** | as F44 |
| F46 | SubtractOverscan.py | 103 | `p` | `'Img # %05d OSCAN [%d:%d, %d:%d] amp %d, noise: %.3f e-/px' % (frameno, x0, x1, y0, y1, ia, sdrs)` | pl>=2; per amp | **no individual save**; saved only inside the `gridplot` (L135) | `"oscan_%05d.png" % frameno`; pl>=3 and `len(plts) > 0` | STD; `circle(legend_label)`; `line(line_color, line_width, legend_label)`; then `gridplot(plts, ncols, plot_width, plot_height, toolbar_location=None)` |
| F47 | SubtractScatteredLight.py | 96 | `p` | `self.action.args.plotlabel + "SCATTERED LIGHT FIT"` | pl>=1 | saved (L110) | `scfnam + ".png"`, `scfnam = "scat_%05d_%s_%s_%s"`; pl>=1 (same block; after prompt/sleep) | STD; `circle(legend_label)`; `line(color, line_width, legend_label)`; `legend.location = "top_left"` |
| F48 | TraceBars.py | 135 | `p` | `plab + 'SPATIAL CONTROL POINTS'` | `do_plot` (pl>=1); inside the "bars found" branch | saved (L147) | `trcfnam + ".png"`, `trcfnam = "bars_%05d_%s_%s_%s"`; `do_plot` (pl>=1; after prompt/sleep) | STD; `scatter(marker='x', size=2, color)`; `scatter(color)` **no `marker`** |

Per primitive: ArcOffsets 2, ExtractArcs 2, FindBars 3, FitCenter 3, GetAtlasLines 1, MakeCube 1, MakeInvsens 7, MakeMasterBias 1, MakeMasterFlat 13, MakeMasterSky 3, ReadAtlas 2, SolveAIT 1, SolveArcs 6, SubtractOverscan 1, SubtractScatteredLight 1, TraceBars 1 = 48 rows in 16 primitives.

Saved rows: there are 20 `save_plot(` calls (section 5). Nineteen save a figure directly: F02, F09, F10, F11, F16, F17, F18, F19, F20, F28, F33, F36, F38, F39, F43, F44, F45, F47, F48. The twentieth (SubtractOverscan.py:135) saves a `gridplot` that wraps the F46 figures. So 29 rows have no direct save; F46 is saved only through the gridplot, and the other 28 are never saved (live-only).

Notes on figure/save behaviour that differ from what the guard suggests:

- `MakeMasterFlat.py:441` and `:455` (and their `bokeh_plot` calls at `:450` and `:466`) have no `plot_level` guard. With `enable_bokeh = False` and a BM-grating flat they raise at `plot_level` 0 as well (bug B9 is stated only for `plot_level >= 1`). With a live session they push and sleep at `plot_level` 0.
- `SolveArcs.py:537` and `:578` call `save_plot` at any `plot_level`, so PNG export (and its browser dependency) runs at `plot_level` 0.
- `SubtractOverscan.py:113` and `:131` push the same figure objects twice; `:131` and `:135` build two separate `gridplot`s from the same `plts` list.
- In `SolveAIT.py` and `SubtractOverscan.py` the save happens before the prompt; in the others it happens after the prompt or sleep (so a pause adds delay before the PNG).

## 2. Feature summary

This is the specification for task 1.1 and the expected set for the coverage scan in 1.2. It was built by reading each primitive and cross-checked with an AST walk (a throwaway script, not kept) over `kcwidrp/primitives/*.py`. The AST walk lists calls whose receiver is a name bound to `figure(...)`; the names `map` and `close` are excluded on purpose (see the note at the end of 2.3).

### 2.1 Keywords passed to `figure()`

| Keyword | Where |
|---|---|
| `title` | all 48 |
| `x_axis_label` | all 48 |
| `y_axis_label` | all 48 |
| `plot_width` | all 48 (value `self.config.instrument.plot_width`) |
| `plot_height` | all 48 (value `self.config.instrument.plot_height`) |
| `tooltips` | MakeCube.py:269 `[("x","$x"),("y","$y"),("value","@image")]`; MakeInvsens.py:349 `[("x","@x{0,0.0}"),("y","@y{0,0.0}")]`; MakeInvsens.py:447 `[("x","@x{0.0}"),("y","@y{0.0}")]` |
| `y_axis_type` | MakeInvsens.py:594, value `'log'` |
| `y_range` | MakeMasterFlat.py:647 (`ylim`, a 2-list) and :875 (`yran`, a 2-list) |

No other keyword is passed. `StartBokeh.py:41` calls `figure()` with none.

### 2.2 Glyph methods and their keywords

Positional arguments are always data (x, y) except for `image`, where the first positional is `[warped]` (a list holding a 2-D array). `quad` takes only keywords.

| Method | Keywords used (union) | Files |
|---|---|---|
| `line` | `color`, `line_color`, `line_width`, `line_dash`, `legend_label`, `alpha`, `source`, `x_range_name` | ArcOffsets, ExtractArcs, FindBars, FitCenter, GetAtlasLines, MakeInvsens, MakeMasterBias, MakeMasterFlat, MakeMasterSky, ReadAtlas, SolveAIT, SolveArcs, SubtractOverscan, SubtractScatteredLight; also in `kcwi_plotting.oplot_slices` (`color`, `line_dash`) |
| `circle` | `size`, `color`, `fill_color`, `line_alpha`, `legend_label`, `x_range_name` | FindBars, MakeMasterFlat, MakeMasterSky, SolveArcs, SubtractOverscan, SubtractScatteredLight |
| `scatter` | `marker`, `color`, `legend_label`, `size` | FindBars, FitCenter, MakeInvsens, TraceBars |
| `diamond` | `color`, `size`, `legend_label` | ArcOffsets, GetAtlasLines, SolveArcs |
| `quad` | `top`, `bottom`, `left`, `right`, `fill_color`, `line_color`, `alpha` | MakeMasterBias |
| `image` | `x`, `y`, `dw`, `dh`, `palette`, `level` | MakeCube |

Per-keyword details:

- `line(source=source)` is used once, with string column names: `p.line('x', 'y', line_color='black', source=source)` (MakeInvsens.py:355).
- `line(alpha=0.7)` is used once (MakeMasterBias.py:146) together with `line_width=4`.
- `line(x_range_name="pix")` and `circle(x_range_name="pix")` appear only in SolveArcs.py:245-258.
- `circle(line_alpha=0.)` appears at MakeMasterSky.py:339, :421, :598 and MakeMasterFlat.py:651, :653, :655, :878, always with `size=1, fill_color=...`.
- `scatter` has no `marker` at FitCenter.py:331 and TraceBars.py:140 (the second `scatter` of that figure). All other `scatter` calls pass `marker='x'`.

### 2.3 Attributes read or assigned on a figure

| Attribute | Use | Where |
|---|---|---|
| `x_range` | assigned a `Range1d` | ReadAtlas.py:214 (positional), GetAtlasLines.py:557 (positional), SolveArcs.py:239 (`start=`, `end=`) |
| `y_range` | assigned a `Range1d` (positional) | GetAtlasLines.py:559, MakeMasterFlat.py:328, :547 |
| `x_range.range_padding` | assigned `0` (chained: `p.x_range.range_padding = p.y_range.range_padding = 0`) | MakeCube.py:275 |
| `y_range.range_padding` | assigned `0` | MakeCube.py:275 |
| `y_range.start` | assigned | MakeMasterBias.py:151 (`0`), MakeMasterSky.py:603, SolveArcs.py:259 (`0`); also through `set_plot_lims` |
| `y_range.end` | assigned | MakeMasterSky.py:604; also through `set_plot_lims` (`ylim[1]`) |
| `x_range.start`, `x_range.end` | assigned only through `set_plot_lims` | see 2.5 |
| `legend.location` | assigned | `"top_left"`: SubtractScatteredLight.py:104, MakeMasterFlat.py:548 (single quotes), :885; `"bottom_center"`: FindBars.py:97, :135, FitCenter.py:394; `"top_center"`: FitCenter.py:366. Every assignment follows glyph calls that passed `legend_label`, so a legend exists |
| `xgrid.grid_line_color` | assigned `None` | ArcOffsets.py:116, FitCenter.py:364, :392, SolveArcs.py:488, :527, :567 (`oplot_slices` does not touch the grid) |
| `extra_x_ranges` | assigned a dict `{"pix": Range1d(start=, end=)}` | SolveArcs.py:240 |
| `add_layout(LinearAxis(x_range_name="pix", axis_label="CCD Y pix"), 'above')` | method | SolveArcs.py:242 |

Nothing is read from a figure (no attribute is used as an rvalue except the chained reads implicit in the assignments above).

Do not list `map` or `close` as figure methods: in `MakeCube.py:244-246` (`p = Pool()`; `p.map(...)`; `p.close()`) and `FitCenter.py:302-304` (`p = get_context("spawn").Pool()`) the name `p` is first bound to a multiprocessing pool, and the figure bound to `p` later in the same function is a different object.

### 2.4 Other names imported from `bokeh`, with keywords used

| Name | Module | Usage |
|---|---|---|
| `figure` | `bokeh.plotting` (all 16 primitives); `bokeh.plotting.figure` (StartBokeh.py, core) | see 2.1 |
| `ColumnDataSource` | `bokeh.plotting` | MakeInvsens.py:345: `ColumnDataSource(data=dict(x=w, y=obsspec))` (keyword `data`). A commented copy is at :434 |
| `Range1d` | `bokeh.models` | positional `Range1d(a, b)`: ReadAtlas.py:214, GetAtlasLines.py:557, :559, MakeMasterFlat.py:328, :547. **Keyword form** `Range1d(start=, end=)`: SolveArcs.py:239, :240 |
| `LinearAxis` | `bokeh.models` | SolveArcs.py:242: `LinearAxis(x_range_name="pix", axis_label="CCD Y pix")` |
| `gridplot` | `bokeh.layouts` | SubtractOverscan.py:131, :135: `gridplot(plts, ncols=(2 if namps > 2 else 1), plot_width=500, plot_height=300, toolbar_location=None)` |
| `basicConfig` | `bokeh.util.logconfig` | FindBars.py:26: `basicConfig(level=logging.ERROR)` (called inside `__init__`; keyword `level`) |
| `pull_session` | `bokeh.client` | StartBokeh.py:31: `pull_session()` (no arguments; a commented call with `session_id='kcwi', url='http://localhost:5006'` is at :30) |
| `column` | `bokeh.layouts` | StartBokeh.py:42: `column(children=[p])`; core/bokeh_plotting.py:38 |

Session and document calls (StartBokeh and `core/bokeh_plotting.py`, outside the primitives' figure code): `session.document.clear()`, `session.document.add_root(c)`, `session.show(c)`, `session.push()`, `session.document.select_one(selector=dict(type=Figure))`, `select_one(selector=dict(type=Column))`, `layout.children.remove(...)`, `layout.children.insert(0, plot)`. Also in core: `from bokeh.plotting import output_file, save` (`bokeh_save`, unused by primitives), `from bokeh.models import Column`, `from bokeh.plotting.figure import Figure, figure`, and in `core/kcwi_plotting.py` `from bokeh.io import export_png` with `export_png(fig, filename=fnam)` and `export_png(fig, filename=fnam, webdriver=driver)`.

### 2.5 What `set_plot_lims` and `oplot_slices` (`core/kcwi_plotting.py`) do to a figure

- `set_plot_lims(fig, xlim=None, ylim=None)` (L108-115). If `xlim` is truthy: assigns `fig.x_range.start = xlim[0]` and `fig.x_range.end = xlim[1]`. If `ylim` is truthy: assigns `fig.y_range.start = ylim[0]` and `fig.y_range.end = ylim[1]`. It never replaces the range object, so it needs the default `DataRange1d`. Called from SolveAIT (F39), MakeInvsens (F14, F15, F16-F19), SolveArcs (F41, F43, F44, F45), ArcOffsets (F02), FitCenter (F09, F10). It is **not** called on a figure whose range was replaced by `Range1d` (SolveArcs F40 uses `Range1d` directly).
- `oplot_slices(fig, yrange)` (L101-105). For `ix` in 1..23 it calls `fig.line([sx, sx], yrange, color='black', line_dash='dashdot')` with `sx = ix * 5 - 0.5`: 23 vertical dash-dot black lines at x = 4.5, 9.5, ... 114.5, with no `legend_label`. Called from SolveArcs (F43, F44, F45), ArcOffsets (F02), FitCenter (F09, F10).
- `get_plot_lims(data, padding=0.05, clip=True)` (L91-98) does not touch a figure; it returns `(min - pad*range, max + pad*range)` (range forced to 1.0 if below 1e-4 and `clip`).
- `save_plot(fig, filename=None)` writes `plots/<filename>` through `export_png` (default file `plots/kcwi_drp_plot.png`), via a cached headless Firefox `webdriver` when `configure_plot_driver` was used, otherwise Bokeh's own browser detection.

### 2.6 Literal values

| Keyword | Literal values used |
|---|---|
| `marker=` | `'x'` only (FindBars.py:93, :131, FitCenter.py:359, :386, MakeInvsens.py:193, :574, :599, :640, TraceBars.py:139) |
| `line_dash=` | `'dashed'`, `'dashdot'`, `'dotted'` (`'dotted'` only at SolveArcs.py:404; `'dashdot'` also inside `oplot_slices`) |
| `palette=` | `"Spectral11"` (MakeCube.py:277) |
| `level=` (on `image`) | `"image"` (MakeCube.py:277) |
| `y_axis_type=` | `'log'` (MakeInvsens.py:594) |
| `legend.location` | `"top_left"`, `"bottom_center"`, `"top_center"` (the lines are listed in 2.3) |
| `add_layout` position | `'above'` (SolveArcs.py:244) |
| `toolbar_location=` (on `gridplot`) | `None` |

### 2.7 Prompts (full list in section 3)

Answers that change the reduction: P14 (MakeInvsens.py:365), P15 (:469), P19 (:655), P25 (MakeMasterFlat.py:471/473), P37 and P38 (NandshuffSubtractSky.py:185, :269), P40 (ReadAtlas.py:224). The eight prompts in the two scripts change which frames are ingested (S01-S08). All other prompts only pace the display or choose which figures are shown.

## 3. Prompt table

Reading guide. "Guard" is the condition under which the prompt is reached, beyond the figure guard. "Default" is the effect of an empty reply. "Pause" is the matching `time.sleep(...plot_pause...)` branch (section 4). Unless stated, any non-matching text is accepted as `<cr>`; closed standard input (EOF) raises `EOFError`, which is not caught anywhere. "B5" means the site tests `'Q' in q.upper()` (any text containing the letter q or Q quits, bug B5). "Quit" always means "stop showing further figures in this loop"; the reduction continues. Section 3.1 covers the primitives (51 logical prompts from 52 `input(` lines), 3.2 covers the two scripts (8 lines).

### 3.1 Primitives

| ID | File:line | Guard | Exact prompt string | Figure on screen | Accepted answers | Default `<cr>` | Changes reduction? | Pause branch | Misbehaviour / bug |
|---|---|---|---|---|---|---|---|---|---|
| P01 | ArcOffsets.py:92 | `do_plot` (pl>=2), `arc_number == next_bar_to_plot` | `Next? <int> or <cr>, q to quit: ` | F01 | `<cr>`: next bar; `<int> n`: jump to bar n; `q`: stop plotting (`do_plot = False`) | next bar (`int('')` raises ValueError, caught) | No (display only) | none | B4 (n <= current bar, negative, or >= NBARS never matches, review ends silently), B5; non-integer text such as `2.5` acts as `<cr>` |
| P02 | ArcOffsets.py:121 | pl>=2 (else branch of pl>=1 block) | `Next? <cr>: ` | F02 | any text = next | next | No | ArcOffsets.py:123 (pl==1) | none |
| P03 | ExtractArcs.py:149 | `do_plot` (pl>=3), `NBARS < 60` branch | `Next? <cr>, q to quit: ` | F03 | `<cr>`: next; `q`: `do_plot = False` | next | No | none | B5 |
| P04 | ExtractArcs.py:215 | `do_plot` (pl>=3), `NBARS >= 60` branch | `Next? <cr>, q to quit: ` | F04 | as P03 | next | No | none | B5 |
| P05 | FindBars.py:99 | pl>=3 (unconditional inside that block) | `Next? <cr>: ` | F05 | any text = next | next | No | none | none |
| P06 | FindBars.py:138 | pl>=2 inside pl>=1 block | `Next? <cr>: ` | F06 | any text = next | next | No | FindBars.py:140 (pl==1) | none |
| P07 | FindBars.py:165 | `do_inter` (pl>=3) | `Next? <cr>, q - quit: ` | F07 | `<cr>`: next bar; `q`: `do_inter = False`, no further centroid plots | next | No | FindBars.py:169 exists but is **dead code** (`else` of `if do_inter` nested inside `if do_inter`) | B5 |
| P08 | FitCenter.py:339 | `do_inter` (pl>=2), `ir == next_bar_to_plot` | `Next? <int> or <cr>, q to quit: ` | F08 | `<cr>`: next; `<int> n`: jump to bar index n; `q`: stop | next (ValueError caught) | No | none | B4, B5 |
| P09 | FitCenter.py:370 | pl>=2 | `Next? <cr>: ` | F09 | any = next | next | No | FitCenter.py:372 (pl==1) | none |
| P10 | FitCenter.py:398 | pl>=2 | `Next? <cr>: ` | F10 | any = next | next | No | FitCenter.py:400 (pl==1) | none |
| P11 | GetAtlasLines.py:562 | pl>=2 | `Next? <cr>: ` | F11 | any = next | next | No | GetAtlasLines.py:564 (pl==1) | none |
| P12 | MakeCube.py:280 | `do_inter` (pl>=3), per slice (24) | `Next? <cr>, q to quit: ` | F12 | `<cr>`: next slice; `q`: `do_inter = False` | next | No | MakeCube.py:284: **live after q** (the remaining slices are still plotted, each followed by a `plot_pause` sleep) | B5 |
| P13 | MakeInvsens.py:196 | `do_plots` (pl>=3), per slice | `Next? <cr>, q - quit: ` | F13 | `<cr>`: next; `q`: `do_plots = False` | next | No | none | B5 |
| P14 | MakeInvsens.py:365 | pl>1, `while not done` | `New values? <float> <float> (or <cr> - done): ` | F14 | `<cr>`: accept current limits, loop ends; two floats `lo hi`: new fit limits (valid if `lo <= hi` and within `[wall0, wall1]`), figure is redrawn and asked again; anything else: both limits reset to WAVGOOD, message "range/order error" or "format error, try again", figure redrawn | accept | **Yes**: sets the wavelength region of the fit (`wl_good`, L381) and so the fitted inverse sensitivity and effective area in `*_invsens.fits` and `*_ea.fits`, and header `INVFW0`/`INVFW1` | none | **B6**: any invalid entry resets *both* limits to WAVGOOD, losing limits accepted earlier. One number only: `IndexError` is caught here (reset), so it does not raise at this prompt; `q` is treated as bad format |
| P15 | MakeInvsens.py:469-470 | pl>=2, `while not done` | `Mask line: wavelength start stop (Ang) comment (str) <float> <float> <str> (<cr> - done): ` | F15 | `<cr>`: done; `w0 w1 [comment...]`: add a line mask if both lie strictly inside `(wlm0, wlm1)` (else message "line mask outside wl range"); non-numeric: "bad line", asked again | done | **Yes**: masks become zero weight in the polynomial fits (L509-530), changing `*_invsens.fits`, `*_ea.fits`; also writes/overwrites a local `<stdfile>.lmsk` (L488-499) | none | **B2**: one number, or whitespace only, raises uncaught `IndexError` (`qstr.split()[1]` / `[0]`); only `ValueError` is caught. Also observed, not in section 13: header lines of an existing `.lmsk` are written back without a newline (L493) |
| P16 | MakeInvsens.py:585 | pl>=2 | `Next? <cr>: ` | F16 | any = next | next | No | MakeInvsens.py:587 (`2. * plot_pause`, pl==1) | none |
| P17 | MakeInvsens.py:606 | pl>=2 | `Next? <cr>: ` | F17 | any = next | next | No | MakeInvsens.py:608 (`2. *`) | none |
| P18 | MakeInvsens.py:625 | pl>=2 | `Next? <cr>: ` | F18 | any = next | next | No | MakeInvsens.py:627 (`2. *`) | none |
| P19 | MakeInvsens.py:655-656 | pl>=2, end of each pass of the fit loop | `Current fit order = %d, New fit order? <int>, <cr> - done: ` (with `% ford`) | F19 (the other three figures of the pass F16-F18 were shown just before) | `<cr>`: done; `<int>`: refit with that polynomial order and repeat all four figures | done | **Yes**: the fit order, so `finvsen`, `fearea`, header `INVFORD`, outputs `*_invsens.fits` and `*_ea.fits` | MakeInvsens.py:686 (pl==1): sets `done = True` then sleeps `2. * plot_pause` | **B3**: `ford = int(qstr)` is assigned before `np.polyfit`; if the fit raises `ValueError` (for example a negative order) the message "Bad fit order, try again" is printed but `ford` keeps the bad value while the fit is the old one. No upper bound. `q` gives "Bad fit order". `LinAlgError` (a `ValueError`) is also swallowed; other exceptions (for example `MemoryError` for a huge order) are not caught |
| P20 | MakeMasterBias.py:154 | pl>=2 | `Next? <cr>: ` | F20 | any = next | next | No | MakeMasterBias.py:156 (pl==1) | none |
| P21 | MakeMasterFlat.py:270 | pl>=2 | `Next? <cr>: ` | F21 | any = next | next | No | MakeMasterFlat.py:272 | none |
| P22 | MakeMasterFlat.py:331 | pl>=2 | `Next? <cr>: ` | F22 | any = next | next | No | MakeMasterFlat.py:333 | none |
| P23 | MakeMasterFlat.py:375 | pl>=2 | `Next? <cr>: ` | F23 | any = next | next | No | MakeMasterFlat.py:377 | none |
| P24 | MakeMasterFlat.py:452 | pl>=2 | `Next? <cr>: ` | F24 | any = next | next | No | MakeMasterFlat.py:454 (taken at pl 0 and 1) | the figure and sleep also occur at pl 0 (no pl>=1 guard) |
| P25 | MakeMasterFlat.py:471 and :473 (one logical prompt: the second call repeats the first inside `while not spk.isnumeric()`) | pl>=2 and `len(peaks) != 1`; also printed first: `Please indicate the integer pixel value of the peak` | `Peak? <int>: ` | F25 | a non-negative integer n: `ipk = n`, used as an index into `xvals` (`apk = xvals[ipk]`); anything else (empty, negative, decimal, text) is re-asked | none (an empty reply re-asks; there is no way to accept the automatic choice) | **Yes**: `apk` is the ledge position; it sets `ratio`, the correction `yfr[i] /= ratio`, hence the master flat (written by MakeMasterFlat) | at pl<2 the same decision is made automatically (L475-478, no pause) | **B1c**: no bounds check; `n >= len(xvals)` raises `IndexError`. Also observed: `str.isnumeric()` is true for characters such as `'²'` where `int()` then raises `ValueError`, uncaught. **B1/B1b** apply to the non-interactive `else` branch (L478), not to this prompt |
| P26 | MakeMasterFlat.py:494 | pl>=2 (always true inside the pl>=3 block) | `Next? <cr>: ` | F26 | any = next | next | No | MakeMasterFlat.py:496 is **dead code** (else of `pl>=2` inside `pl>=3`) | none |
| P27 | MakeMasterFlat.py:551 | pl>=2 | `Next? <cr>: ` | F27 | any = next | next | No | MakeMasterFlat.py:553 | none |
| P28 | MakeMasterFlat.py:658 | pl>=2 | `Next? <cr>: ` | F28 | any = next | next | No | MakeMasterFlat.py:660 | none |
| P29 | MakeMasterFlat.py:727 | pl>=2 | `Next? <cr>: ` | F29 | any = next | next | No | MakeMasterFlat.py:729 | none |
| P30 | MakeMasterFlat.py:747 | pl>=2 | `Next? <cr>: ` | F30 | any = next | next | No | MakeMasterFlat.py:749 | none |
| P31 | MakeMasterFlat.py:766 | pl>=2 | `Next? <cr>: ` | F31 | any = next | next | No | MakeMasterFlat.py:768 | none |
| P32 | MakeMasterFlat.py:784 | pl>=2 | `Next? <cr>: ` | F32 | any = next | next | No | MakeMasterFlat.py:786 | none |
| P33 | MakeMasterFlat.py:888 | pl>=2 | `Next? <cr>: ` | F33 | any = next | next | No | MakeMasterFlat.py:890 | none |
| P34 | MakeMasterSky.py:347 | pl>=2 | `Next? <cr>: ` | F34 | any = next | next | No | MakeMasterSky.py:349 | none |
| P35 | MakeMasterSky.py:435 | pl>=2 | `Next? <cr>: ` | F35 | any = next | next | No | MakeMasterSky.py:437 | none |
| P36 | MakeMasterSky.py:607 | pl>=2 | `Next? <cr>: ` | F36 | any = next | next | No | MakeMasterSky.py:609 | none |
| P37 | NandshuffSubtractSky.py:185-186 | pl>=2; BLUE camera, non-standard N&S configuration (`skyrow0 >= 10`), inside `while skyscl < 0.` | `Enter scale factor for sky to match obj (float): ` | none (this primitive makes no figure; whatever figure was pushed last stays up) | a float >= 0: sky scale factor. Negative or non-numeric: warning "Invalid input ... try again" and re-asked | **none**: an empty reply is `float('')`, a `ValueError`, and re-asks. With pl<2 the scale is silently 1.0 | **Yes**: `sky *= skyscl` (L205); changes the nod-and-shuffle subtracted image (`intk`), the sky and object images, and header `NASSCL` | none (no figure; pl<2 takes `skyscl = 1.0`, L194) | **B7**: 0 is accepted (no sky subtraction, no warning). Also observed: `nan` and `inf` are accepted (`nan < 0.` is False) |
| P38 | NandshuffSubtractSky.py:269-270 | as P37 for the RED camera, `objrow0 >= 10` (L268) | `Enter scale factor for sky to match obj (float): ` | none | as P37 | as P37 | **Yes**: as P37 | none (`skyscl = 1.0` at L278) | **B7**, `nan`/`inf` as P37. Also observed (not in section 13, not run): in this RED branch `sky[0, objrow0, :] = 0.` and `obj[0, objrow0, :] = 0.` (L303, L305) index a 2-D array with three indices, which would raise `IndexError` after the prompt is answered |
| P39 | ReadAtlas.py:193 | pl>=2 | `Next? <cr>: ` | F37 | any = next | next | No | ReadAtlas.py:195 (pl==1) | none |
| P40 | ReadAtlas.py:224 | pl>=2, inside `while q:` | `Enter: <cr> - next, new offset (int px): ` | F38 | `<cr>`: accept and leave the loop; integer n: use n px as the arc-atlas offset (`offset_wav = n * refdisp`) and redraw F38; non-integer: prints "Try again: integer pixel values accepted" and asks again | accept the current offset | **Yes**: stored as `args.offset_pix` and `args.offset_wave` (L247-248); `FitCenter.py:239` subtracts `offset_wave` from the starting wavelengths `p0`, so it changes the wavelength solution (`twkcoeff`, used by SolveArcs and GetAtlasLines) and all products derived from it | ReadAtlas.py:233 (pl==1; sets `q = None`) | `q` (letter) is not special here: it is a non-integer and is re-asked. Negative integers and out-of-range values are accepted without check |
| P41 | SolveAIT.py:172 | pl>=2 (`do_plot` block) | `Next? <cr>: ` (assigned to `_`) | F39 | any = next | next | No | **none**: at pl==1 there is no pause | the PNG is saved (L170) before the prompt |
| P42 | SolveArcs.py:262 | `do_inter` (pl>=3), `ib == next_bar_to_plot`, inside `try` | `ptitle + " - Next? <cr>, q to quit: "` where `ptitle = " Bar# %d - line %3d/%3d: xc = %.1f, Wave = %9.2f"` | F40 | `<cr>`: next line; `q`: `do_inter = False` | next | No | none | B5 |
| P43 | SolveArcs.py:412 | `master_inter` (pl>=2), `ib == next_bar_to_plot` | `Next? <cr>: ` | F41 | any = next | next | No | none | none |
| P44 | SolveArcs.py:438 | as P43 | `Next? <int> or <cr>, q - quit: ` | F42 | `<cr>`: next bar; `<int> n`: jump to bar n; `q`: `master_inter = False` | next bar (`int('')` raises ValueError, caught) | No | none | B4, B5 |
| P45 | SolveArcs.py:493 | pl>=2, per coefficient | `Next? <cr>: ` | F43 | any = next | next | No | SolveArcs.py:495 (pl==1) | none |
| P46 | SolveArcs.py:533 | pl>=2 inside pl>=1 | `Next? <cr>: ` | F44 | any = next | next | No | SolveArcs.py:535 | none |
| P47 | SolveArcs.py:573 | pl>=2 inside pl>=1 | `Next? <cr>: ` | F45 | any = next | next | No | SolveArcs.py:575 | none |
| P48 | SubtractOverscan.py:116 | pl>=3, per amp | `Next? <cr>: ` | F46 (the amp just drawn) | any = next | next | No | SubtractOverscan.py:118 (pl==2) | none |
| P49 | SubtractOverscan.py:140 | pl>=2 (always true inside the pl>=3 block) | `Next? <cr>: ` | the `gridplot` of the F46 figures (no row of its own) | any = next | next | No | SubtractOverscan.py:142 is **dead code** | `bokeh_clear(self.context.bokeh_session)` follows at L143 and raises when the session is `None` (**B10**) |
| P50 | SubtractScatteredLight.py:107 | pl>=2 | `Next? <cr>: ` | F47 | any = next | next | No | SubtractScatteredLight.py:109 (pl==1) | none |
| P51 | TraceBars.py:144 | pl>=2 | `Next? <cr>: ` | F48 | any = next | next | No | TraceBars.py:146 (pl==1) | none |

### 3.2 Scripts

None of these has a `plot_level` guard. No figure is on screen (they run before the framework starts). The test in each is a substring match on the file name or path, so a directory name containing `kr` or `kb` triggers it (**B8**).

| ID | File:line | Guard | Exact prompt string | Accepted answers | Default `<cr>` | Changes reduction? | Misbehaviour / bug |
|---|---|---|---|---|---|---|---|
| S01 | scripts/reduce_kcwi.py:437 | `-f/--frames` branch: `args.blue and 'kr' in frame.lower()` | `Proceed? <cr>=yes or Q=quit: ` | `<cr>` or any text without q: keep this frame and go on; text containing q/Q: `frames = []`, `break` (frame list emptied, `ingest_data(None, [], False)` is still called) | proceed | Yes: decides which frames are ingested | B8 (prompt repeats for every offending file and consumes lines of a piped answers file), B5. Behaviour of `ingest_data` with an empty list was not checked (`keckdrpframework` is not installed here) |
| S02 | scripts/reduce_kcwi.py:443 | `-f` branch: `args.red and 'kb' in frame.lower()` | same | same | proceed | Yes | B8, B5 |
| S03 | scripts/reduce_kcwi.py:462 | `-l/--list` branch: `"#" not in frame`, `args.blue and 'kr' in frame.lower()` | same | same (`break` leaves the file loop only) | proceed | Yes | B8, B5 |
| S04 | scripts/reduce_kcwi.py:468 | `-l` branch: `args.red and 'kb' in frame.lower()` | same | same | proceed | Yes | B8, B5 |
| S05 | scripts/kcwi_rti.py:412 | `-f` branch: `args.blue and ('kr' in frame or 'KR' in frame)` | same | text containing q/Q: `frames = []`, `break` | proceed | **No effect even for q**: L423 calls `framework.ingest_data(None, args.frames, False)` with the original list, not `frames`, so answering Q does not stop the run | B8, B5; quit is ineffective in this branch (observed, not in section 13) |
| S06 | scripts/kcwi_rti.py:418 | `-f` branch: `args.red and ('kb' in frame or 'KB' in frame)` | same | same | proceed | as S05 | B8, B5, as S05 |
| S07 | scripts/kcwi_rti.py:434 | `-l` branch, `args.blue and ('kr' ...)` | same | q/Q: `frames = []`, `break` | proceed | Yes (`ingest_data(None, frames, False)` at L445 receives the emptied list) | B8, B5 |
| S08 | scripts/kcwi_rti.py:440 | `-l` branch, `args.red and ('kb' ...)` | same | same | proceed | Yes | B8, B5 |

Other notes on input behaviour that apply to every row: the `isnumeric()` loop in P25 is the only place that re-asks on an empty reply by design; all "any = next" prompts accept text and ignore it; `EOFError` from a closed standard input is never caught, so a run started with piped input that runs out of lines ends the primitive with an exception (the framework logs it and continues, section 6.1 of the upgrade plan).

## 4. Sleep list (36 sites)

Every `time.sleep(...)` that uses `plot_pause`, in the order of `grep -rn "plot_pause" kcwidrp/primitives | grep sleep` grouped by file. The paired prompt is the `input()` in the sibling branch of the same `if plot_level >= N / else` statement.

| # | Sleep site | Paired prompt | Note |
|---|---|---|---|
| 1 | SubtractScatteredLight.py:109 | P50 (L107) | |
| 2 | MakeMasterSky.py:349 | P34 (L347) | |
| 3 | MakeMasterSky.py:437 | P35 (L435) | |
| 4 | MakeMasterSky.py:609 | P36 (L607) | |
| 5 | TraceBars.py:146 | P51 (L144) | |
| 6 | FindBars.py:140 | P06 (L138) | |
| 7 | FindBars.py:169 | P07 (L165) | dead code (nested `if do_inter` / `else`) |
| 8 | ReadAtlas.py:195 | P39 (L193) | |
| 9 | ReadAtlas.py:233 | P40 (L224) | also sets `q = None` to leave the loop |
| 10 | SubtractOverscan.py:118 | P48 (L116) | |
| 11 | SubtractOverscan.py:142 | P49 (L140) | dead code (`pl>=2` is always true inside the `pl>=3` block) |
| 12 | MakeMasterFlat.py:272 | P21 (L270) | |
| 13 | MakeMasterFlat.py:333 | P22 (L331) | |
| 14 | MakeMasterFlat.py:377 | P23 (L375) | |
| 15 | MakeMasterFlat.py:454 | P24 (L452) | runs at pl 0 as well |
| 16 | MakeMasterFlat.py:496 | P26 (L494) | dead code |
| 17 | MakeMasterFlat.py:553 | P27 (L551) | |
| 18 | MakeMasterFlat.py:660 | P28 (L658) | |
| 19 | MakeMasterFlat.py:729 | P29 (L727) | |
| 20 | MakeMasterFlat.py:749 | P30 (L747) | |
| 21 | MakeMasterFlat.py:768 | P31 (L766) | |
| 22 | MakeMasterFlat.py:786 | P32 (L784) | |
| 23 | MakeMasterFlat.py:890 | P33 (L888) | |
| 24 | MakeInvsens.py:587 | P16 (L585) | `2. * plot_pause` |
| 25 | MakeInvsens.py:608 | P17 (L606) | `2. * plot_pause` |
| 26 | MakeInvsens.py:627 | P18 (L625) | `2. * plot_pause` |
| 27 | MakeInvsens.py:686 | P19 (L655) | `2. * plot_pause`; follows `done = True` |
| 28 | MakeCube.py:284 | P12 (L280) | live branch after the user answers q (not dead) |
| 29 | SolveArcs.py:495 | P45 (L493) | |
| 30 | SolveArcs.py:535 | P46 (L533) | |
| 31 | SolveArcs.py:575 | P47 (L573) | |
| 32 | ArcOffsets.py:123 | P02 (L121) | |
| 33 | GetAtlasLines.py:564 | P11 (L562) | |
| 34 | FitCenter.py:372 | P09 (L370) | |
| 35 | FitCenter.py:400 | P10 (L398) | |
| 36 | MakeMasterBias.py:156 | P20 (L154) | |

Prompts that have no `time.sleep` sibling (they sit in blocks that only run at the level that prompts, or are loops that ask every time): P01, P03, P04, P05, P08, P13, P14, P15, P25, P37, P38, P41, P42, P43, P44. Two of these matter for the replacement viewer because a figure is pushed at a lower level with no pause and no prompt: F39 (SolveAIT.py, pl 1: pushed and saved, nothing waits) and F25 (MakeMasterFlat.py "Deriv", always pushed, and followed at pl<2 by the automatic peak choice with no pause). Dead sleeps are rows 7, 11 and 16.

## 5. Counts and verification

Commands were run from the repository root at `a360685`.

| Quantity | Expected | Actual | Match |
|---|---|---|---|
| `input(` lines in `kcwidrp/primitives` | 52 | 52 | yes |
| Primitives containing `input(` | 17 | 17 | yes |
| `input(` lines in `scripts/reduce_kcwi.py` + `scripts/kcwi_rti.py` | 8 | 8 (4 + 4) | yes |
| `bokeh_plot(` calls in `kcwidrp/primitives` | 49 | 49 | yes |
| `save_plot(` calls in `kcwidrp/primitives` | 20 | 20 | yes |
| `time.sleep(...plot_pause...)` sites | 36 | 36 | yes |
| `figure(` calls excluding `StartBokeh.py` (figure rows) | 48 | 48 | yes |
| Primitives with a figure row | 16 | 16 | yes |
| `figure(` calls including `StartBokeh.py:41` | 49 | 49 | yes |

No count differs from its expected value. How the counts reconcile:

- 52 `input(` lines = 51 logical prompts (P01-P51) + 1, because `MakeMasterFlat.py:471` and `:473` are one logical prompt (P25). The scripts add 8 more (S01-S08), for 59 `input(` lines and 59 logical prompts in total.
- 49 `bokeh_plot(` calls = 48 figure rows + 1: `SubtractOverscan.py:131` pushes the `gridplot` of the F46 figures in addition to the push of each F46 figure at `:113`. Every other figure is pushed exactly once.
- 20 `save_plot(` calls = 19 directly saved figures + the one `gridplot` save at `SubtractOverscan.py:135`. The plan's "21" for `save_plot` counts the `def save_plot` line in `kcwi_plotting.py`, which is outside `kcwidrp/primitives`.
- 36 sleep sites include 3 dead ones (rows 7, 11, 16 of section 4) and 1 that is live only after the user answers q (row 28).
- 48 + 1 = 49 `figure(` calls: the extra one is `StartBokeh.py:41`, the placeholder for the Bokeh session.

Per-file `figure(` counts (from `grep -rn "figure(" kcwidrp/primitives | awk -F: '{print $1}' | uniq -c`): ArcOffsets 2, ExtractArcs 2, FindBars 3, FitCenter 3, GetAtlasLines 1, MakeCube 1, MakeInvsens 7, MakeMasterBias 1, MakeMasterFlat 13, MakeMasterSky 3, ReadAtlas 2, SolveAIT 1, SolveArcs 6, StartBokeh 1 (not a row), SubtractOverscan 1, SubtractScatteredLight 1, TraceBars 1.

Commands and output:

```
$ git rev-parse --short HEAD
a360685
$ grep -rn "input(" kcwidrp/primitives | wc -l
      52
$ grep -rln "input(" kcwidrp/primitives | wc -l
      17
$ grep -n "input(" kcwidrp/scripts/reduce_kcwi.py kcwidrp/scripts/kcwi_rti.py
kcwidrp/scripts/reduce_kcwi.py:437:                qstr = input('Proceed? <cr>=yes or Q=quit: ')
kcwidrp/scripts/reduce_kcwi.py:443:                qstr = input('Proceed? <cr>=yes or Q=quit: ')
kcwidrp/scripts/reduce_kcwi.py:462:                        qstr = input('Proceed? <cr>=yes or Q=quit: ')
kcwidrp/scripts/reduce_kcwi.py:468:                        qstr = input('Proceed? <cr>=yes or Q=quit: ')
kcwidrp/scripts/kcwi_rti.py:412:                qstr = input('Proceed? <cr>=yes or Q=quit: ')
kcwidrp/scripts/kcwi_rti.py:418:                qstr = input('Proceed? <cr>=yes or Q=quit: ')
kcwidrp/scripts/kcwi_rti.py:434:                        qstr = input('Proceed? <cr>=yes or Q=quit: ')
kcwidrp/scripts/kcwi_rti.py:440:                        qstr = input('Proceed? <cr>=yes or Q=quit: ')
$ grep -n "input(" kcwidrp/scripts/reduce_kcwi.py kcwidrp/scripts/kcwi_rti.py | wc -l
       8
$ grep -rn "bokeh_plot(" kcwidrp/primitives | wc -l
      49
$ grep -rn "save_plot(" kcwidrp/primitives | wc -l
      20
$ grep -rn "plot_pause" kcwidrp/primitives | grep -c "time.sleep"
36
$ grep -rn "figure(" kcwidrp/primitives | wc -l
      49
$ grep -rn "figure(" kcwidrp/primitives | grep -v StartBokeh | wc -l
      48
$ grep -rln "figure(" kcwidrp/primitives | grep -v StartBokeh | wc -l
      16
$ grep -rln "input(" kcwidrp/primitives | wc -l
      17
$ grep -rn "gridplot(" kcwidrp/primitives | wc -l
       2
```

## 6. Easy-to-miss items from bokeh plan 0.2

| Item | Present in the feature summary | Source |
|---|---|---|
| `figure(y_range=...)` twice in MakeMasterFlat | Yes, 2.1 | MakeMasterFlat.py:647 (F28), :875 (F33) |
| `figure(y_axis_type='log')` in MakeInvsens | Yes, 2.1 and 2.6 | MakeInvsens.py:594 (F17) |
| `line_alpha` on `circle` (MakeMasterSky, MakeMasterFlat) | Yes, 2.2 | MakeMasterSky.py:339, :421, :598; MakeMasterFlat.py:651, :653, :655, :878 |
| `scatter` with no `marker` (FitCenter) | Yes, 2.2 and 2.6 | FitCenter.py:331 (also TraceBars.py:140, a second case not listed in the plan) |
| `Range1d(start=, end=)` keyword form | Yes, 2.4 | SolveArcs.py:239, :240 (the other `Range1d` uses are positional) |
| `y_range.end` via `set_plot_lims` | Yes, 2.3 and 2.5 | `kcwi_plotting.py:115`; also direct at MakeMasterSky.py:604 |
| `image(..., palette=, level=)` | Yes, 2.2 and 2.6 | MakeCube.py:276-277 |
| `map` / `close` not listed as figure methods | Yes, note in 2.3 | MakeCube.py:244-246, FitCenter.py:302-304 |
