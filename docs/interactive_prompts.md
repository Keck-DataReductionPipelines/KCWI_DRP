# Interactive prompt catalogue (task P.0)

Purpose: the text that every interactive prompt of the DRP prints in PR3a
(upgrade plan §8). For each logical prompt of `docs/plot_inventory.md` §3
(P01–P51 in the primitives, S01–S08 in the scripts) this file gives the
one-line `what` printed before the prompt and the long `help` printed when
the user types `?` or `help`. P.3 and P.4 paste them into the call sites as
the `what=` and `help=` keyword arguments of the `kcwidrp.core.kcwi_interactive`
helpers (§8.2, "Interface fixed for parallel work"). No text registry is
created: the strings live at the call site.

Source read for this file: branch `upgrade/pr3a-prompts` at `cb7444c`
(PR1, the Matplotlib port, is merged). Line numbers below are from that
commit; where they differ from the inventory (taken at `a360685`) both are
given. Colours were checked against the figure code at each site; a glyph
call with no colour is drawn in the wrapper's default `#1f77b4` (blue,
`kcwidrp/core/_plot_render.py:35`), and a marker given only `fill_color`
keeps a default-blue outline (`_plot_render.py:238`). `scatter` without a
`marker` draws circles (`_plot_render.py:236`).
Since 2026-11-09 the plots and their help texts follow the
[plot colour scheme](#plot-colour-scheme).

Order: sections and tables follow the order in which a reduction reaches
the prompts: the scripts first (before the framework starts), then the
primitives in the order of the event table in
`kcwidrp/pipelines/kcwi_pipeline.py` (bias, continuum bars, arcs, flats,
objects, nod-and-shuffle). Within a primitive, prompts are in line order.
SolveAIT is not in the event table (AIT data only); it sits after
SolveArcs, its counterpart for standard data. The IDs keep the inventory's numbering
(alphabetical by file) because tests and other documents refer to them.

## How to read a section

Each section has:

- **Site**: `file:line` of today's `input(` call, and the inventory's line
  if different.
- **Helper**: which function the site will use, with the arguments that
  matter. `ask_value` and `ask_yes_quit` take a `prompt` and no `what`; at
  those sites the call site prints `what` itself just before the call
  (or P.1 adds a `what=` argument; P.4 decides).
- **Figure**: the inventory figure ID on screen when the prompt is asked.
- **Changes reduction**: yes/no, and for yes, which outputs.
- **what**: one line, at most 100 characters. Where it contains `%d`/`%s`,
  the site fills them in with `%` formatting before the call.
- **help**: the long text, in the template of §8.2:
  *What you see* / *What to look for* / *Your options* / *Default*. It is
  written for an observer at the telescope or an astronomer reducing data
  later. Lines are at most 76 characters so the text fits in a terminal.
- For value prompts: **Accepted input** (exactly what today's code accepts,
  D12; stricter checks that change results for valid input wait for PR3b),
  **Prompt** (proposed prompt string) and **Effect** (the format string for
  the INFO log line; `ask_value` formats it with `{value}` and prefixes
  `"<who>: "`).
- **Reviewer note** (where needed): for the maintainer, not part of the user
  text. Notes say where a prompt is reached only through a known bug's path,
  or where the code does something the text has to work around.
- **Status**: filled in by the maintainer during the P.0 human check:
  `approved` (text correct as written) or `corrected` (text changed in
  review; the correction is made in this file before P.3/P.4 paste it).

The help texts describe the behaviour **after PR3a**: `?`/`help` reprints
the help, `q`/`quit` (exact, any case) quits a loop, a bar number outside
the valid range is rejected with the range shown (B4), and invalid input to
a value prompt re-prompts keeping the last accepted value (B1c, B2, B6). At
end of input (Ctrl-D, or a piped answers file that runs out) every prompt
takes its default (D4). "Default" in each help text is what `<cr>` does,
which is also what the reduction does when it is not interactive
(`plot_level` below the prompt's level).

Common words used in the texts: *bar* = one of the 120 continuum-bar
traces (5 per slice, numbered 0–119, slice = bar // 5); *REFBAR* = the
reference bar (config, default 57); *slice boundaries* = the black
dash-dot vertical lines that `oplot_slices` draws every 5 bars
(`kcwi_plotting.py:807`).

## Plot colour scheme

Every plot uses one colour (and line style) per role, so the same thing
looks the same from plot to plot (maintainer, 2026-11-09). The help texts
describe the plots in this scheme; the per-site list is in
[Colours named per site](#colours-named-per-site).

| Role | Colour / style |
|---|---|
| Data being inspected (points, spectra, profiles, histogram bars) | blue |
| Second data series of the same kind (raw background, corrected data, continuum-bar spectrum, red-end slice) | cyan |
| Reference or comparison data (atlas, REFBAR arc, the star's reference spectrum, reference slice) | black |
| Fit, model, smoothing | red solid |
| Mean / ±σ, ±RMS of a fit or of the data | red solid / red dashed |
| Mean of a subset (per-slice mean, P46) | black, thick, over the subset only |
| CWAVE | orange, dash-dot-dot (used for nothing else) |
| WAVGOOD | green solid verticals |
| WAVALL | grey solid verticals |
| Found or adopted value (peak, centroid, adopted dispersion, kept lines) | magenta; P11 kept lines green (maintainer, 2026-11-09; P11 has no WAVGOOD) |
| Predicted or calculated value (grating equation, predicted ledge, ATOFF) | black dashed |
| Limits and windows (fit limits, fit windows, buffer zones, WavLim, Cor lim, source masks) | purple; a second set dashed or dotted |
| Line masks | purple dashed pairs |
| Thresholds / zero lines | grey dashed / grey dotted |
| Masked points, lines rejected before a fit | grey x |
| Lines rejected as outliers in residual clipping | gold |
| GetAtlasLines rejection classes (P11) | cyan WinRej, pink FitRej, olive ParRej, brown DupRej, grey FntRej |
| Slice boundaries | black dash-dot |
| Images (P12) | Spectral11 with a colour bar; NaN grey; pixel boundaries as black ticks at the top and bottom |

Rules:

- Where a plot shows several series of the same kind, each with its own
  fit (P28–P32), each fit is drawn in its series' colour; there the
  reference slice is black, the blue-end slice blue and the red-end slice
  cyan.
- Every plot whose x axis covers the whole wavelength range shows
  WAVGOOD in green, also in pixel space where a wavelength solution
  exists when the plot is drawn. WAVGOOD is first computed in SolveGeom,
  so plots made before it (P01, P03, P04, P11, P39, P40, P43, P44) have
  none. P50 (CCD rows) has none: a row is a different wavelength in each
  slice. P12 shows WAVGOOD as green horizontals (y is wavelength).
- Every line or marker a help text names has a legend entry.
- Markers have an outline of their fill colour.
- In the help text each colour item is on its own line; continuation
  lines are indented by two more spaces.

## Check these first

The P.0 author was least sure of the science in these prompts. Each
is marked ⚠ in the index below and has a **⚠ Check first** line under
its heading.

| Prompts | What to check |
|---|---|
| [P37](#p37--sky-scale-blue-camera), [P38](#p38--sky-scale-red-camera) | Why a sky scale is needed in aborted nod-and-shuffle (the agent's own reading) |
| [P24](#p24--bm-ledge), [P25](#p25--bm-ledge-peak-choice) | The description of the BM ledge and of a good peak |
| [P02](#p02--bar-offsets-summary), [P09](#p09--central-wavelengths), [P45](#p45--wavelength-coefficients) | The claim that a good result repeats smoothly from slice to slice |
| [P14](#p14--fit-wavelength-limits), [P15](#p15--line-masks) | The features suggested for masking (Balmer lines, O2 B/A bands at 6870/7600 Å) |
| [P20](#p20--bias-noise-histogram), [P43](#p43--residuals-of-one-bar), [P47](#p47--fit-rms-per-bar) | No numbers for read noise or RMS, only "close to the amp's usual value" |
| [P50](#p50--scattered-light) | What the ±180 px columns used for scattered light sample on the detector |

## Index

| ID | Site | Helper | Changes reduction | Status |
|---|---|---|---|---|
| S01 | reduce_kcwi.py:422 | ask_yes_quit | **Yes** (frames ingested) | |
| S02 | reduce_kcwi.py:428 | ask_yes_quit | **Yes** | |
| S03 | reduce_kcwi.py:447 | ask_yes_quit | **Yes** | |
| S04 | reduce_kcwi.py:453 | ask_yes_quit | **Yes** | |
| S05 | kcwi_rti.py:400 | ask_yes_quit | No (see note, B11) | |
| S06 | kcwi_rti.py:406 | ask_yes_quit | No (see note, B11) | |
| S07 | kcwi_rti.py:422 | ask_yes_quit | **Yes** | |
| S08 | kcwi_rti.py:428 | ask_yes_quit | **Yes** | |
| P48 | SubtractOverscan.py:116 | pause | No | |
| P49 | SubtractOverscan.py:140 | pause | No | |
| ⚠ P20 | MakeMasterBias.py:154 | pause | No | |
| P05 | FindBars.py:97 | pause | No | |
| P06 | FindBars.py:136 | pause | No | |
| P07 | FindBars.py:163 | ask_next | No | |
| P51 | TraceBars.py:144 | pause | No | |
| P03 | ExtractArcs.py:149 | ask_next | No | |
| P04 | ExtractArcs.py:215 | ask_next | No | |
| P01 | ArcOffsets.py:92 | ask_next (n_items) | No | |
| ⚠ P02 | ArcOffsets.py:121 | pause | No | |
| P39 | ReadAtlas.py:193 | pause | No | |
| P40 | ReadAtlas.py:224 | ask_value | **Yes** | |
| P08 | FitCenter.py:339 | ask_next (n_items) | No | |
| ⚠ P09 | FitCenter.py:370 | pause | No | |
| P10 | FitCenter.py:398 | pause | No | |
| P11 | GetAtlasLines.py:562 | pause | No | |
| P42 | SolveArcs.py:262 | ask_next | No | |
| ⚠ P43 | SolveArcs.py:412 | pause | No | |
| P44 | SolveArcs.py:438 | ask_next (n_items) | No | |
| ⚠ P45 | SolveArcs.py:493 | pause | No | |
| P46 | SolveArcs.py:533 | pause | No | |
| ⚠ P47 | SolveArcs.py:573 | pause | No | |
| P41 | SolveAIT.py:172 | pause | No | |
| P12 | MakeCube.py:280 | ask_next | No | |
| ⚠ P50 | SubtractScatteredLight.py:107 | pause | No | |
| P21 | MakeMasterFlat.py:270 | pause | No | |
| P22 | MakeMasterFlat.py:331 | pause | No | |
| P23 | MakeMasterFlat.py:375 | pause | No | |
| ⚠ P24 | MakeMasterFlat.py:452 | pause | No | |
| ⚠ P25 | MakeMasterFlat.py:476/478 | ask_value | **Yes** | |
| P26 | MakeMasterFlat.py:500 | pause | No | |
| P27 | MakeMasterFlat.py:557 | pause | No | |
| P28 | MakeMasterFlat.py:664 | pause | No | |
| P29 | MakeMasterFlat.py:733 | pause | No | |
| P30 | MakeMasterFlat.py:753 | pause | No | |
| P31 | MakeMasterFlat.py:772 | pause | No | |
| P32 | MakeMasterFlat.py:790 | pause | No | |
| P33 | MakeMasterFlat.py:894 | pause | No | |
| P34 | MakeMasterSky.py:347 | pause | No | |
| P35 | MakeMasterSky.py:435 | pause | No | |
| P36 | MakeMasterSky.py:607 | pause | No | |
| P13 | MakeInvsens.py:196 | ask_next | No | |
| ⚠ P14 | MakeInvsens.py:364 | ask_value | **Yes** | |
| ⚠ P15 | MakeInvsens.py:468 | ask_value | **Yes** | |
| P16 | MakeInvsens.py:584 | pause | No | |
| P17 | MakeInvsens.py:605 | pause | No | |
| P18 | MakeInvsens.py:624 | pause | No | |
| P19 | MakeInvsens.py:654 | ask_value | **Yes** | |
| ⚠ P37 | NandshuffSubtractSky.py:185 | ask_value | **Yes** | |
| ⚠ P38 | NandshuffSubtractSky.py:269 | ask_value | **Yes** (header only, see note) | |

---

## Scripts

All eight prompts run before the framework starts, with no figure and no
`plot_level` guard. The helper is `ask_yes_quit(logger, prompt, help,
who=...)`: `<cr>` proceeds, `q`/`quit` quits, other text asks again. The
line printed before the prompt today (`Blue channel requested, but red
files in list` or the red variant) becomes `what`.

Common reviewer notes: B8 (PR3b): the channel test is a substring match on
the whole path (`'kr' in frame.lower()`; in `kcwi_rti` `'kr'`/`'KR'`), so a
directory name containing `kr`/`kb` triggers the prompt, and it repeats for
every matching file. Today any answer without a `q` proceeds (e.g. `y`);
after PR3a such text re-prompts (B5/§8.2 interface).

### S01 — reduce_kcwi -f, red file in a blue reduction

- Site: `kcwidrp/scripts/reduce_kcwi.py:422` (inventory :437)
- Helper: `ask_yes_quit(logger, prompt="Proceed? <cr>=yes or q=quit", help, who="reduce_kcwi")`
- Changes reduction: **Yes**: quit empties the frame list; `ingest_data` is called with no frames, so nothing is reduced.
- Log on quit (suggested): `reduce_kcwi: user quit at the channel check; no frames ingested`
- Status:

**what**
```text
Blue channel requested (-b), but this file name looks like a red-channel (kr) file.
```

**help**
```text
What you see:
  No plot. The reduction was started for the blue channel, and the
  path of one of the frames given with -f contains "kr" (red files
  are named kr<date>_<n>.fits).
What to look for:
  Check the file: a red frame in a blue reduction will not combine
  with the blue calibrations. A directory name containing "kr" also
  triggers this check (known issue); then it is safe to proceed.
Your options:
  <cr>      keep this frame and go on (the next file is checked too).
  q / quit  stop: no frames are reduced.
Default:
  <cr> proceeds with this frame.
```

### S02 — reduce_kcwi -f, blue file in a red reduction

- Site: `reduce_kcwi.py:428` (inventory :443)
- Helper: as S01
- Changes reduction: **Yes** (as S01)
- Status:

**what**
```text
Red channel requested (-r), but this file name looks like a blue-channel (kb) file.
```

**help**
```text
What you see:
  No plot. The reduction was started for the red channel, and the path
  of one of the frames given with -f contains "kb" (blue files are
  named kb<date>_<n>.fits).
What to look for:
  Check the file: a blue frame in a red reduction will not combine
  with the red calibrations. A directory name containing "kb" also
  triggers this check (known issue); then it is safe to proceed.
Your options:
  <cr>      keep this frame and go on (the next file is checked too).
  q / quit  stop: no frames are reduced.
Default:
  <cr> proceeds with this frame.
```

### S03 — reduce_kcwi -l, red file in a blue reduction

- Site: `reduce_kcwi.py:447` (inventory :462)
- Helper: as S01
- Changes reduction: **Yes**: quit empties the list read so far and stops reading the file; `ingest_data` gets no frames (the `<list>_ingest` stamp file is still written).
- Status:

**what**
```text
Blue channel requested (-b), but this file in the list looks like a red-channel (kr) file.
```

**help**
```text
What you see:
  No plot. The reduction was started for the blue channel, and a line
  of the file list given with -l contains "kr" (red files are named
  kr<date>_<n>.fits).
What to look for:
  Check the list: a red frame in a blue reduction will not combine
  with the blue calibrations. A directory name containing "kr" also
  triggers this check (known issue); then it is safe to proceed.
Your options:
  <cr>      keep this frame and go on (the next line is checked too).
  q / quit  stop: no frames from the list are reduced.
Default:
  <cr> proceeds with this frame.
```

### S04 — reduce_kcwi -l, blue file in a red reduction

- Site: `reduce_kcwi.py:453` (inventory :468)
- Helper: as S01
- Changes reduction: **Yes** (as S03)
- Status:

**what**
```text
Red channel requested (-r), but this file in the list looks like a blue-channel (kb) file.
```

**help**
```text
What you see:
  No plot. The reduction was started for the red channel, and a line
  of the file list given with -l contains "kb" (blue files are named
  kb<date>_<n>.fits).
What to look for:
  Check the list: a blue frame in a red reduction will not combine
  with the red calibrations. A directory name containing "kb" also
  triggers this check (known issue); then it is safe to proceed.
Your options:
  <cr>      keep this frame and go on (the next line is checked too).
  q / quit  stop: no frames from the list are reduced.
Default:
  <cr> proceeds with this frame.
```

### S05 — kcwi_rti -f, red file in a blue reduction

- Site: `kcwidrp/scripts/kcwi_rti.py:400` (inventory :412)
- Helper: `ask_yes_quit(..., who="kcwi_rti")`
- Changes reduction: **No, today**: B11 (PR3b): the `-f` branch ingests `args.frames`, not the emptied list, so `q` stops only the checking; all frames are still reduced.
- Status:
- Reviewer note: the help text states the B11 behaviour honestly; when B11 is fixed in PR3b, replace the `q / quit` line with the S01 wording.

**what**
```text
Blue channel requested (-b), but this file name looks like a red-channel (kr) file.
```

**help**
```text
What you see:
  No plot. The real-time reduction was started for the blue channel,
  and the path of one of the frames given with -f contains "kr" or
  "KR" (red files are named kr<date>_<n>.fits).
What to look for:
  Check the file: a red frame in a blue reduction will not combine
  with the blue calibrations. A directory name containing "kr" also
  triggers this check (known issue); then it is safe to proceed.
Your options:
  <cr>      keep this frame and go on (the next file is checked too).
  q / quit  stop checking. Note: in this version all frames given
            with -f are still reduced (known bug); press Ctrl-C to
            stop the run.
Default:
  <cr> proceeds with this frame.
```

### S06 — kcwi_rti -f, blue file in a red reduction

- Site: `kcwi_rti.py:406` (inventory :418)
- Helper: as S05
- Changes reduction: **No, today** (B11, as S05)
- Status:

**what**
```text
Red channel requested (-r), but this file name looks like a blue-channel (kb) file.
```

**help**
```text
What you see:
  No plot. The real-time reduction was started for the red channel,
  and the path of one of the frames given with -f contains "kb" or
  "KB" (blue files are named kb<date>_<n>.fits).
What to look for:
  Check the file: a blue frame in a red reduction will not combine
  with the red calibrations. A directory name containing "kb" also
  triggers this check (known issue); then it is safe to proceed.
Your options:
  <cr>      keep this frame and go on (the next file is checked too).
  q / quit  stop checking. Note: in this version all frames given
            with -f are still reduced (known bug); press Ctrl-C to
            stop the run.
Default:
  <cr> proceeds with this frame.
```

### S07 — kcwi_rti -l, red file in a blue reduction

- Site: `kcwi_rti.py:422` (inventory :434)
- Helper: as S05
- Changes reduction: **Yes**: quit empties the list; `ingest_data` gets no frames (the `<list>_ingest` stamp file is still written).
- Status:

**what**
```text
Blue channel requested (-b), but this file in the list looks like a red-channel (kr) file.
```

**help**
```text
What you see:
  No plot. The real-time reduction was started for the blue channel,
  and a line of the file list given with -l contains "kr" or "KR"
  (red files are named kr<date>_<n>.fits).
What to look for:
  Check the list: a red frame in a blue reduction will not combine
  with the blue calibrations. A directory name containing "kr" also
  triggers this check (known issue); then it is safe to proceed.
Your options:
  <cr>      keep this frame and go on (the next line is checked too).
  q / quit  stop: no frames from the list are reduced.
Default:
  <cr> proceeds with this frame.
```

### S08 — kcwi_rti -l, blue file in a red reduction

- Site: `kcwi_rti.py:428` (inventory :440)
- Helper: as S05
- Changes reduction: **Yes** (as S07)
- Status:

**what**
```text
Red channel requested (-r), but this file in the list looks like a blue-channel (kb) file.
```

**help**
```text
What you see:
  No plot. The real-time reduction was started for the red channel,
  and a line of the file list given with -l contains "kb" or "KB"
  (blue files are named kb<date>_<n>.fits).
What to look for:
  Check the list: a blue frame in a red reduction will not combine
  with the red calibrations. A directory name containing "kb" also
  triggers this check (known issue); then it is safe to proceed.
Your options:
  <cr>      keep this frame and go on (the next line is checked too).
  q / quit  stop: no frames from the list are reduced.
Default:
  <cr> proceeds with this frame.
```

---

## SubtractOverscan

### P48 — overscan, per amplifier

- Site: `SubtractOverscan.py:116`; prompt at pl >= 3, sleep at pl 2 (:118)
- Helper: `pause(...)` inside a kept `plot_level >= 3` test; the pl 2 sleep must stay (P.3: `pause` alone would prompt at pl 2).
- Figure: F46 (pl >= 2; the amplifier just drawn)
- Changes reduction: No
- Status:

**what**
```text
Overscan of one amplifier along the readout, with the fit used for the noise estimate.
```

**help**
```text
What you see:
  Title "Img # <n> OSCAN [x0:x1, y0:y1] amp <a>, noise: <s> e-/px".
  Blue circles: median of the overscan columns for each row (counts).
  Red thick line: polynomial fit; the residuals give the read noise
    in the title (OSCNRN<a>).
  The value subtracted from the amplifier is a single number, the
  median of the overscan region (OSCNVAL<a>).
What to look for:
  Good: a flat, smooth level the fit follows; a ramp only in the first
  or last 50 rows (excluded from the fit).
  Bad: jumps, strong gradients or a noise value much higher than usual
  for this amplifier: an electronics problem in this frame.
Your options:
  <cr>  continue (next amplifier). Display only.
Default:
  <cr> continues. At plot_level 2 each amplifier is shown for
  plot_pause seconds; below plot_level 2 nothing is drawn.
```

### P49 — overscan, all amplifiers

- Site: `SubtractOverscan.py:140` (sleep :142 is dead code)
- Helper: `pause(...)` inside the existing `plot_level >= 3` block.
- Figure: the grid of the F46 figures (saved as `oscan_<frame>.png` before the prompt)
- Changes reduction: No
- Status:
- Reviewer note: `bokeh_clear` (:143) now returns when the session is `None` (B10 fixed in T.11, `bokeh_plotting.py:38`).

**what**
```text
Overscan of all amplifiers of this frame side by side.
```

**help**
```text
What you see:
  The overscan plots of all amplifiers of this frame in one grid, each
  titled with its noise.
  Blue circles: row medians.
  Red line: fit.
What to look for:
  Good: all amplifiers flat and smooth, with similar noise.
  Bad: one amplifier with jumps, gradients or much higher noise than
  the others.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. Below plot_level 3 this grid is not shown.
```

## MakeMasterBias

### P20 — bias noise histogram

> ⚠ **Check first:** No numbers for read noise or RMS, only "close to the amp's usual value".

- Site: `MakeMasterBias.py:154` (sleep :156)
- Helper: `pause(...)`
- Figure: F20 (pl >= 1; one per amplifier; saved as `bias_..._amp<n>_rdnoise.png` after the prompt)
- Changes reduction: No
- Status:

**what**
```text
Read-noise histogram of one amplifier, from the difference of two bias frames.
```

**help**
```text
What you see:
  Title "BIAS NOISE amp <n> = <sigma>" (e-).
  Blue bars: histogram of (bias 2 - bias 3) x gain / sqrt(2) over
    amplifier n, within the 3.5 sigma clipping range.
  Red thick line: Gaussian of the measured sigma (centred on 0, scaled
    to the histogram peak).
  Red dashed verticals: +- sigma.
  Sigma is stored as BIASRN<n> in the master bias header.
What to look for:
  Good: a symmetric, Gaussian histogram centred on zero that matches
  the red curve; sigma close to the usual read noise of this amp in
  this readout mode (slow/fast, gain).
  Bad: a broad, skewed or spiky histogram, or an offset from zero:
  pattern noise, or a bad bias frame in the stack.
Your options:
  <cr>  continue (next amplifier). Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

## FindBars

### P05 — bar peak search, each try

- Site: `FindBars.py:97` (inventory :99); no sleep branch
- Helper: `pause(...)` inside the existing `plot_level >= 3` block (always prompts there).
- Figure: F05 (pl >= 3)
- Changes reduction: No
- Status:

**what**
```text
Continuum-bar profile across one CCD row, with the peaks found above the threshold.
```

**help**
```text
What you see:
  Title "BARS TRACE at = <row>, TRY <n>": the continuum-bars image
  across CCD x, median of a few rows around <row>.
  Blue line: the profile.
  Magenta x: peaks found above the threshold.
  Grey dashed: threshold (mean + 0.5 x standard deviation).
What to look for:
  Good: exactly one magenta x on each of the 120 bars (24 groups of 5),
  none on noise or between bars.
  Bad: missing or extra peaks. Then the DRP moves to another row and
  tries again (up to 5 tries; one plot per try). If no try finds 120
  bars, FindBars logs an error and the geometry cannot be solved.
Your options:
  <cr>  continue. Display only: the reduction is not changed.
Default:
  <cr> continues. Below plot_level 3 this plot is not shown.
```

### P06 — bars found at the middle row

- Site: `FindBars.py:136` (inventory :138; sleep at :138)
- Helper: `pause(...)`
- Figure: F06 (pl >= 1, only when 120 bars were found)
- Changes reduction: No
- Status:

**what**
```text
Continuum-bar profile at the row used for bar finding: all bars found.
```

**help**
```text
What you see:
  Title "BARS MID TRACE Thresh = <value>": the continuum-bars profile
  across CCD x at the row where all bars were found.
  Blue line: profile.
  Magenta x: the bar peaks.
  Grey dashed: peak threshold.
What to look for:
  Good: one magenta x on top of every bar, 24 groups of 5, regular
  spacing.
  Bad: a magenta x on a noise spike or between two bars; this shifts
  the bar numbering and every later step that uses it.
Your options:
  <cr>  continue. Display only: the reduction is not changed.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P07 — per-bar centroid

- Site: `FindBars.py:163` (inventory :165)
- Helper: `ask_next(self.logger, what, help, who="FindBars")`; `QUIT` sets `do_inter = False`.
- Figure: F07 (pl >= 3)
- Changes reduction: No
- Status:
- Reviewer note: the sleep at :167 is dead code (inventory §4 row 7); P.3 can drop it.

**what**
```text
Profile of this bar across CCD x and its centroid.
```

**help**
```text
What you see:
  Title "BAR <n> CENTRD = <x>": bar n across CCD x, minimum subtracted.
  Blue line and circles: the profile ("Bar Trace").
  Magenta vertical: the centroid, the x position used for this bar
    (its height spans threshold to peak of the raw profile, so it need
    not touch the blue curve).
What to look for:
  Good: one clean, roughly symmetric peak with the magenta line at its
  centre.
  Bad: a double or truncated peak, or a flat profile: the centroid, and
  the trace that starts from it, will be off.
Your options:
  <cr>      show the next bar.
  q / quit  stop showing centroids (all are still computed).
  Display only: no answer changes the reduction.
Default:
  <cr> shows the next bar. Below plot_level 3 these plots are not shown.
```

## TraceBars

### P51 — bar traces

- Site: `TraceBars.py:144` (sleep :146)
- Helper: `pause(...)`
- Figure: F48 (pl >= 1; saved as `bars_<frame>_....png` after the prompt)
- Changes reduction: No
- Status:

**what**
```text
Continuum-bar traces: control points used for the geometric solution.
```

**help**
```text
What you see:
  Title "SPATIAL CONTROL POINTS" on CCD x, y.
  Blue x: centroids of each bar traced up and down the CCD from the
    middle row.
  Magenta circles: the bar positions found on the middle row.
What to look for:
  Good: 120 continuous traces (24 groups of 5) running the length of
  the image, each starting from a magenta circle, without gaps or
  jumps.
  Bad: a trace that stops early, jumps to a neighbouring bar, or is
  missing: the geometry of that slice will be wrong.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

---

## ExtractArcs

### P03 — extracted arc, AIT branch

- Site: `ExtractArcs.py:149`
- Helper: `ask_next(self.logger, what, help, who="ExtractArcs")`; `QUIT` sets `do_plot = False`.
- Figure: F03 (pl >= 3, `NBARS < 60`)
- Changes reduction: No
- Status:
- Reviewer note: reached only with `NBARS < 60`, i.e. the KCRM AIT configuration ("This is strictly for AIT data!"). Normal KCWI reductions (NBARS = 120) see P04.

**what**
```text
Arc spectrum extracted along this bar's trace, with its background estimate.
```

**help**
```text
What you see:
  Title "ARC # <n>": the arc spectrum extracted along bar n's trace,
  against CCD y pixel.
  Blue: arc flux summed across the extraction window.
  Cyan: background, the window's minimum times its width, row by row.
  Red: smoothed background (Savitzky-Golay); this is subtracted.
What to look for:
  Good: sharp arc lines on a low background; the red curve runs
  smoothly under the blue continuum between lines.
  Bad: red curve following line peaks or lying well above the
  inter-line level (the window is off the bar: check the trace).
Your options:
  <cr>      show the next bar.
  q / quit  stop showing extracted arcs (all bars are still extracted).
  Display only: no answer changes the reduction.
Default:
  <cr> shows the next bar. Below plot_level 3 these plots are not shown.
```

### P04 — extracted arc, standard branch

- Site: `ExtractArcs.py:215`
- Helper: `ask_next(self.logger, what, help, who="ExtractArcs")`; `QUIT` sets `do_plot = False`.
- Figure: F04 (pl >= 3, `NBARS >= 60`)
- Changes reduction: No
- Status:

**what**
```text
Arc spectrum extracted at this bar (warped arc image), with its fitted background.
```

**help**
```text
What you see:
  Title "ARC # <n>": the arc spectrum of bar n, the median across the
  bar in the geometrically warped arc image, against CCD y pixel.
  Blue: arc spectrum.
  Red: background model, a cubic fitted to the minima of 16 sections
    of the spectrum; it is subtracted from the arc.
What to look for:
  Good: sharp arc lines; the red curve follows the low level between
  lines and stays below the line-free continuum.
  Bad: red curve cutting through the spectrum or bending up at the ends
  (e.g. a section without a line-free minimum); a featureless or very
  faint blue spectrum (bar missed by the trace).
Your options:
  <cr>      show the next bar.
  q / quit  stop showing extracted arcs (all bars are still extracted).
  Display only: no answer changes the reduction.
Default:
  <cr> shows the next bar. Below plot_level 3 these plots are not shown.
```

## ArcOffsets

### P01 — per-bar arc cross-correlation

- Site: `kcwidrp/primitives/ArcOffsets.py:92`
- Helper: `ask_next(self.logger, what, help, n_items=nbars, current=arc_number, who="ArcOffsets")`; `QUIT` sets `do_plot = False`, an int sets `next_bar_to_plot`, `NEXT` sets it to `arc_number + 1`.
- Figure: F01 (pl >= 2)
- Changes reduction: No
- Status:

**what**
```text
Arc spectrum of this bar, shifted by its cross-correlation offset, over the reference bar.
```

**help**
```text
What you see:
  Title "BAR OFFSET for Arc: <bar> Slice: <slice> = <offset>".
  Black: arc spectrum of the reference bar (REFBAR) along CCD y.
  Blue: arc spectrum of this bar, shifted by the offset (px) that
    maximises its cross-correlation with the reference bar.
  Both are tapered at the ends (Tukey window, TUKEYALPHA).
What to look for:
  Good: blue and black arc lines coincide line for line over the whole
  range; different relative line strengths are normal.
  Bad: blue lines displaced from black everywhere, or aligned in one
  part only. A bright line near an end can pull the correlation; a
  larger TUKEYALPHA reduces its weight.
Your options:
  <cr>      show the next bar.
  <number>  jump to that bar (after the current one, up to 119);
            bars in between are not shown. Other numbers are refused.
  q / quit  stop showing bars; offsets are still found for all bars.
  Display only: no answer changes the reduction.
Default:
  <cr> shows the next bar. Below plot_level 2 these plots are not shown.
```

### P02 — bar offsets summary

> ⚠ **Check first:** The claim that a good result repeats smoothly from slice to slice.

- Site: `ArcOffsets.py:121` (sleep at :123)
- Helper: `pause(self.config.instrument, self.logger, what=..., help=..., who="ArcOffsets")`
- Figure: F02 (pl >= 1; saved as `arc_<frame>_..._ta<alpha>_baroffs.png` after the prompt)
- Changes reduction: No
- Status:

**what**
```text
Offset of every bar's arc spectrum from the reference bar (px), by bar number.
```

**help**
```text
What you see:
  Title "BAR OFFSETS": the offset (px along CCD y) of each bar's arc
  spectrum from the reference bar, against bar number.
  Blue diamonds: the offsets.
  Black dash-dot lines: slice boundaries (every 5 bars).
What to look for:
  Good: the reference bar is at 0; the five bars of a slice have similar
  offsets and the pattern changes smoothly from slice to slice.
  Bad: an isolated point far from its neighbours is a failed
  cross-correlation for that bar (check that bar in the next steps).
  These offsets are the starting points of the per-bar wavelength fits.
Your options:
  <cr>  continue. Display only: the reduction is not changed.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn. The PNG is saved either way
  at plot_level >= 1.
```

## ReadAtlas

### P39 — arc–atlas cross-correlation

- Site: `ReadAtlas.py:193` (sleep :195)
- Helper: `pause(...)`
- Figure: F37 (pl >= 1)
- Changes reduction: No
- Status:

**what**
```text
Cross-correlation of the reference-bar arc with the lamp atlas: offset found.
```

**help**
```text
What you see:
  Title "ATLAS OFFSET = <n> px": cross-correlation of the reference-bar
  arc spectrum (on its calculated wavelengths) with the atlas spectrum,
  against offset in atlas pixels (central third of the range).
  Blue: cross-correlation.
  Magenta vertical ("Peak"): offset found.
  Black dashed vertical: offset requested with ATOFF (only if set).
What to look for:
  Good: one sharp peak clearly above the rest, with the magenta line on
  it.
  Bad: several similar peaks or a broad hump: the offset is ambiguous.
  Check the overlay in the next plot and set the offset there.
Your options:
  <cr>  continue to the overlay. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P40 — atlas offset entry

- Site: `ReadAtlas.py:224` (pl 1: sleep and leave the loop, :233)
- Helper: `ask_value(self.logger, prompt, help, parse=int, default=None, effect, who="ReadAtlas")` in the existing `while q:` loop; `None` leaves it; an int sets `offset_pix`, `offset_wav = offset_pix * refdisp` and redraws.
- Figure: F38 (pl >= 1, redrawn each pass; the last pass is saved as `arc_..._atoff.png`)
- Changes reduction: **Yes**: stored as `offset_pix`/`offset_wave`; FitCenter subtracts `offset_wave` from every bar's starting wavelength, so it changes the central fits, the atlas line list (`arc_..._atlines.txt`), the per-bar solutions, the geometry `*_geom.pkl` from SolveGeom, and the wavelength scale of every cube reduced with this arc. No header keyword records it; the INFO line "Final arc-atlas offset" does.
- Prompt: `New offset <int px> (<cr> = accept)`
- Accepted input (D12): anything `int()` reads (`12`, `-3`, `+4`); no range check. Decimals and text (including `q`): "Try again: integer pixel values accepted", asked again.
- Effect: `atlas offset set by user to {value} px (calculated %d px)`
- Status:

**what**
```text
Reference-bar arc over the lamp atlas: accept the offset or enter a new one (atlas px).
```

**help**
```text
What you see:
  Title "ATLAS OFFSET = <A> Ang (<n> px)".
  Blue: reference-bar arc spectrum (bar REFBAR), normalised, on its
    calculated wavelengths minus the current offset.
  Black: atlas spectrum of the lamp, normalised.
  Orange dash-dot-dot: CWAVE.
What to look for:
  Good: blue and black lines coincide line for line over the range
  (relative strengths differ; small misalignments at the ends are
  fixed later, since the dispersion is not yet refined).
  Bad: a constant shift between the two patterns. Enter a new offset:
  a larger offset moves the blue spectrum to shorter wavelengths; the
  title gives the offset in A and px, i.e. the atlas A per px.
Your options:
  <cr>      accept the offset shown and continue.
  <number>  try that offset (integer atlas px); the plot is redrawn and
            you are asked again.
  The offset sets the starting wavelengths of every bar's fit, so it
  changes the whole wavelength solution and all cubes using this arc.
Default:
  <cr> keeps the offset shown: the calculated one (or ATOFF if set)
  unless you changed it. Below plot_level 2 that offset is used.
```

## FitCenter

### P08 — central dispersion scan, per bar

- Site: `FitCenter.py:339`
- Helper: `ask_next(self.logger, what, help, n_items=len(results), current=ir, who="FitCenter")`; `QUIT` sets `do_inter = False`.
- Figure: F08 (pl >= 2)
- Changes reduction: No
- Status:

**what**
```text
Atlas cross-correlation peak against trial central dispersion for this bar.
```

**help**
```text
What you see:
  Title "CENTRAL DISPERSION FIT for Bar: <n> Slice: <s>".
  x: trial central dispersion (A/px), within 5 % of the calculated one.
  y: peak of the cross-correlation of this bar's central arc spectrum
  with the atlas at that dispersion.
  Blue circles and line: the scan ("Data").
  Magenta vertical: dispersion adopted (maximum of the interpolated
    scan, searched in the middle third of the range).
  Black dashed vertical: dispersion calculated from the grating
    equation.
What to look for:
  Good: one clear maximum, magenta line on it, close to the black line.
  Bad: flat, noisy or multi-peaked scan, or the maximum at the edge of
  the middle third: the atlas offset (ReadAtlas) or the arc is suspect.
Your options:
  <cr>      show the next bar.
  <number>  jump to that bar (after the current one, up to 119).
  q / quit  stop showing bars. Display only: nothing is changed.
Default:
  <cr> shows the next bar. Below plot_level 2 these plots are not shown.
```

### P09 — central wavelengths

> ⚠ **Check first:** The claim that a good result repeats smoothly from slice to slice.

- Site: `FitCenter.py:370` (sleep :372)
- Helper: `pause(...)`
- Figure: F09 (pl >= 1; saved as `arc_..._cwaves.png` after the prompt)
- Changes reduction: No
- Status:

**what**
```text
Central wavelength found for every bar, against bar number.
```

**help**
```text
What you see:
  Title "CENTRAL VALUES", y "Central Wavelength (A)".
  Blue x: the central wavelength fitted for each bar.
  Orange dash-dot-dot: CWAVE, the grating central wavelength from the
    header.
  Black dash-dot lines: slice boundaries.
What to look for:
  Good: values close to CWAVE, a smooth repeating pattern from slice to
  slice, no isolated points.
  Bad: a bar far from its neighbours has a wrong central fit; its final
  wavelength solution starts from a wrong guess.
Your options:
  <cr>  continue. Display only: the reduction is not changed.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P10 — central dispersions

- Site: `FitCenter.py:398` (sleep :400)
- Helper: `pause(...)`
- Figure: F10 (pl >= 1; saved as `arc_..._cdisp.png`)
- Changes reduction: No
- Status:

**what**
```text
Central dispersion found for every bar, against bar number.
```

**help**
```text
What you see:
  Title "CENTRAL VALUES", y "Central Dispersion (A)".
  Blue x: the central dispersion (A/px) fitted for each bar.
  Black dashed line: dispersion calculated from the grating equation.
  Black dash-dot lines: slice boundaries.
What to look for:
  Good: all bars close to each other and to the black line, varying
  smoothly across the field.
  Bad: isolated bars far off, or values at the edge of the scanned range
  (5 % either side of the black line): a failed central fit.
Your options:
  <cr>  continue. Display only: the reduction is not changed.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

## GetAtlasLines

### P11 — atlas line list

- Site: `GetAtlasLines.py:562` (sleep :564)
- Helper: `pause(...)`
- Figure: F11 (pl >= 1, only when no LINELIST is configured; saved as `arc_..._atlines.png`)
- Changes reduction: No
- Status:
- Reviewer note (not a prompt bug, not in §13): when more than 400 lines survive, `rej_fnt_w = refws[:hlim]` (L508) is taken after `refws` was already trimmed (L506), so the magenta "FntRej" diamonds mark the wrong lines (a subset of the kept ones). The text below says only what the colour means.

**what**
```text
Atlas lines selected for the wavelength fit, and those rejected, over the arc.
```

**help**
```text
What you see:
  Title "ATLAS LINES Ngood = <n>, Nrej = <m>", normalised flux.
  Blue: reference-bar arc spectrum on its central-fit wavelengths.
  Black: atlas spectrum of the arc lamp.
  Large green diamonds: atlas lines kept for the fit ("Kept").
  Small diamonds, rejected lines (plotted at the arc line's height):
    cyan: WinRej (no clean window).
    pink: FitRej (Gaussian fit failed).
    olive: ParRej (position/width out of bounds).
    brown: DupRej (duplicate).
    grey: FntRej (faintest dropped when > 400 lines).
  Purple verticals: wavelength limits used ("WavLim").
  x range: the WavLim range plus 3 % on each side.
What to look for:
  Good: blue and black spectra line up; kept lines sit on clean atlas
  peaks spread over the whole range between the purple lines.
  Bad: blue and black offset (atlas offset wrong); kept lines missing
  from part of the range, which weakens the solution there.
Your options:
  <cr>  continue. Display only (the list is written to
        arc_..._atlines.txt either way).
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

## SolveArcs

### P42 — individual arc line fits

- Site: `SolveArcs.py:262`; today's prompt string is prefixed with `ptitle`
- Helper: `ask_next(self.logger, what, help, who="SolveArcs")` with `what` built from `ptitle`; `QUIT` sets `do_inter = False`.
- Figure: F40 (pl >= 3, bar `next_bar_to_plot` only)
- Changes reduction: No
- Status:

**what** (filled in at the site with `ib, iw + 1, len(at_wave), peak, aw`)
```text
Line fit: bar %d, line %d/%d, arc peak at %.1f px, atlas %.2f A.
```

**help**
```text
What you see:
  Title "ATLAS/ARC LINE FITS Bar# <n> - line <i>/<N> ...". Two x axes:
  wavelength (bottom, atlas) and CCD y pixel (top, arc), matched over
  this window by the central wavelength fit.
  Atlas:
    Black line and circles: atlas spectrum.
    Black dashed vertical: the atlas line.
  Arc:
    Blue circles: smoothed data in the fit window.
    Red line: interpolation.
    Magenta vertical: interpolated peak, the position used.
    Magenta dashed: centroid.
    Magenta dotted: Gaussian peak.
What to look for:
  Good: one clean arc line whose three magenta verticals coincide, at
  the same place as the black atlas line.
  Bad: blended or distorted lines, or an arc line far from the atlas
  line (misidentification); the residual clipping may then reject it.
Your options:
  <cr>      show the next line of this bar.
  q / quit  stop the line plots for all bars (bar plots continue).
  Display only: no answer changes the reduction.
Default:
  <cr> shows the next line. Below plot_level 3 these are not shown.
```

### P43 — residuals of one bar

> ⚠ **Check first:** No numbers for read noise or RMS, only "close to the amp's usual value".

- Site: `SolveArcs.py:412`; no sleep branch (inside `master_inter`)
- Helper: `pause(...)` inside the existing `master_inter` block (pl >= 2, so it always prompts).
- Figure: F41 (pl >= 2, bar `next_bar_to_plot` only)
- Changes reduction: No
- Status:

**what**
```text
Wavelength-fit residuals of this bar, with kept and rejected lines.
```

**help**
```text
What you see:
  Title "RESIDUALS for Bar <n>, Slice <s>, RMS = <r> A (<p> px),
  N = <k>": for each line, the wavelength the fit gives at the arc
  line's pixel minus the atlas wavelength (A). RMS is that of the kept
  lines; the px value divides it by this bar's dispersion.
  Blue diamonds: lines kept.
  Gold diamonds: lines rejected as outliers in the iterative clipping.
  Grey x (at 0): lines rejected before the fit (no clean window,
    failed Gaussian, centroid off).
  Grey dotted: zero.
  Red dashed: +- RMS.
  Purple verticals ("WavLim"): the wavelength limits of the atlas line
    list (as in P11); lines are only used between them.
  Orange dash-dot-dot: CWAVE.
  No WAVGOOD lines: WAVGOOD is only known after SolveGeom.
What to look for:
  Good: blue points scatter randomly about zero over the whole range;
  RMS (A) well below this bar's dispersion, so the px value in the
  title is well below 1.
  Bad: a trend or curve with wavelength (misidentified lines or too
  low a polynomial order), large gaps, many gold or grey lines (often
  a wrong atlas offset).
Your options:
  <cr>  show the fit overlay for this bar. Display only.
Default:
  <cr> continues. Below plot_level 2 this plot is not shown.
```

### P44 — fit overlay of one bar, bar selection

- Site: `SolveArcs.py:438`
- Helper: `ask_next(self.logger, what, help, n_items=len(self.context.arcs), current=ib, who="SolveArcs")`; `QUIT` sets `master_inter = False`; an int sets `next_bar_to_plot`; `NEXT` sets it to `ib + 1`.
- Figure: F42 (pl >= 2)
- Changes reduction: No
- Status:

**what**
```text
Arc spectrum of this bar on its final wavelength solution, over the atlas.
```

**help**
```text
What you see:
  Title "ATLAS/ARC FIT for Bar <n> ...".
  Blue line: this bar's arc on its final wavelength solution.
  Blue diamonds: lines used in the fit, at their fitted wavelengths.
  Black: atlas spectrum (scaled).
  Magenta diamonds ("Kept"): every line of the atlas line list.
  Gold: lines rejected as residual outliers.
  Grey x: lines rejected before the fit.
  Orange dash-dot-dot: CWAVE.
What to look for:
  Good: blue arc lines fall on the black atlas lines over the whole
  range, with blue diamonds on most magenta ones.
  Bad: blue lines shifted from black in part of the range: the
  solution for this bar is wrong there.
Your options:
  <cr>      show the next bar.
  <number>  jump to that bar (after the current one, up to 119); at
            plot_level 3 the line plots follow the jump.
  q / quit  stop the bar review (line plots stop too). The summary
            plots (coefficients, lines per bar, RMS) still follow.
  Display only: no answer changes the reduction.
Default:
  <cr> shows the next bar. Below plot_level 2 these are not shown.
```

### P45 — wavelength coefficients

> ⚠ **Check first:** The claim that a good result repeats smoothly from slice to slice.

- Site: `SolveArcs.py:493` (sleep :495)
- Helper: `pause(...)`
- Figure: F43 (pl >= 1; one per coefficient, C0 first; saved as `arc_..._coef<k>.png` after the prompt)
- Changes reduction: No
- Status:

**what**
```text
One coefficient of the wavelength solution for every bar.
```

**help**
```text
What you see:
  Title "COEF <k> OF C0-C<n> VALUES <Ck> = <mean> +- <std>": coefficient
  k of each bar's wavelength polynomial in CCD pixel, against bar
  number. The polynomial has n + 1 coefficients, C0 to Cn: C0 is the
  wavelength (A) at the first pixel, C1 the dispersion (A/px), higher
  terms the curvature. One plot per coefficient.
  Blue diamonds: values.
  Black dash-dot lines: slice boundaries.
What to look for:
  Good: a smooth pattern within and between slices, no isolated
  points.
  Bad: a bar far off its neighbours has a poor solution (compare its
  N lines and RMS in the next plots).
Your options:
  <cr>  next coefficient (or the next plot). Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P46 — number of lines per bar

- Site: `SolveArcs.py:533` (sleep :535)
- Helper: `pause(...)`
- Figure: F44 (pushed at pl >= 1; saved as `arc_..._nlines.png` at every plot_level)
- Changes reduction: No
- Status:

**what**
```text
Number of arc lines used in the wavelength fit of every bar.
```

**help**
```text
What you see:
  Title "FIT STATS <Nlns> = <mean> +- <std>": number of lines in each
  bar's final fit, against bar number.
  Blue diamonds ("Bars"): values.
  Red ("Mean"): mean of all bars.
  Red dashed: mean +- 1 standard deviation.
  Black thick, one per slice ("Slice mean"): mean of the slice's 5 bars.
  Black dash-dot: slice boundaries.
What to look for:
  Good: similar numbers for all bars and slices.
  Bad: bars or slices with far fewer lines than the rest; their
  solutions are less constrained, especially at the ends of the range.
  A red mean far below what this grating usually gives means the line
  identification failed: check the atlas offset (P39, P40) and the
  line list (P11).
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn (the PNG is still saved).
```

### P47 — fit RMS per bar

> ⚠ **Check first:** No numbers for read noise or RMS, only "close to the amp's usual value".

- Site: `SolveArcs.py:573` (sleep :575)
- Helper: `pause(...)`
- Figure: F45 (pushed at pl >= 1; saved as `arc_..._resid.png` at every plot_level)
- Changes reduction: No
- Status:

**what**
```text
RMS of the wavelength fit of every bar (Angstrom).
```

**help**
```text
What you see:
  Title "FIT STATS <RMS> = <mean> +- <std> A (<p> px)": RMS residual
  (A) of each bar's final fit, against bar number. A bar's RMS is that
  of its kept lines' residuals: the wavelength the fit gives at the arc
  line's pixel minus the atlas wavelength (as in P43). The px value
  divides the mean by the mean dispersion.
  Blue diamonds ("Bars"): values.
  Red ("Mean"): mean.
  Red dashed: mean +- 1 standard deviation.
  Black dash-dot: slice boundaries.
What to look for:
  Good: uniformly low RMS, well below the dispersion, so the px value
  in the title is well below 1.
  Bad: isolated bars or whole slices with high RMS: check them with
  the bar review (plot_level 2) before trusting those slices.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn (the PNG is still saved).
```

## SolveAIT

### P41 — AIT bar spectra

- Site: `SolveAIT.py:172`; no sleep branch (at pl 1 the figure is pushed and saved with no pause)
- Helper: `pause(...)` inside a kept `plot_level >= 2` guard, so that plot_level 1 does not gain a sleep.
- Figure: F39 (pl >= 1; saved as `arc_..._bar<n>.png` **before** the prompt)
- Changes reduction: No
- Status:
- Reviewer note: SolveAIT is used only for KCRM AIT data ("only in case there is a need to reprocess the AIT data"); normal reductions never reach this prompt.

**what**
```text
AIT data: arc and bar spectrum of this trace on the output wavelength grid.
```

**help**
```text
What you see:
  Title "ARC # <n>": spectra of AIT trace n, resampled onto the output
  wavelength grid.
  Blue: arc spectrum.
  Cyan: continuum-bar spectrum.
  Green verticals: good wavelength range (WAVGOOD).
  Grey verticals: full wavelength range (WAVALL).
What to look for:
  Good: sharp arc lines at the right wavelengths, smooth bar spectrum,
  both covering the green range.
  Bad: smeared or doubled lines (poor solution for this trace), or the
  data not reaching the green limits.
Your options:
  <cr>  continue (next trace). Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown and saved with no
  pause; at plot_level 0 it is not drawn.
```

## MakeCube

### P12 — warped slice images

- Site: `MakeCube.py:280`
- Helper: `ask_next(self.logger, what, help, who="MakeCube")` in the existing `if do_inter:` branch; `QUIT` sets `do_inter = False`. The `else` sleep (:284) stays live after a quit (P.3 keeps it, via `pause` or as today).
- Figure: F12 (pl >= 3; 24 slices)
- Changes reduction: No
- Status:

**what**
```text
One slice of the cube after the geometric and wavelength transformation.
```

**help**
```text
What you see:
  Title "WARPED Slice <n>": slice n resampled onto the cube grid.
  x: position along the slice (cube px); y: wavelength (A), shorter
  wavelengths at the bottom.
  Colour: signal (e-) in each cube pixel; the colour bar on the right
    gives the e- for each colour, from the slice's lowest value
    (purple) to its highest (dark red).
  Black ticks at the top and bottom ("Pixel"): one at each boundary
    between pixels along the slice.
  Green horizontals ("WavGood"): WAVGOOD, the wavelength range with
    good data in every slice.
  Grey blocks: no data (NaN).
What to look for:
  Good: sky and arc lines run straight across the slice (constant
  wavelength); a point source is a straight vertical streak; the slice
  edges are clean.
  Bad: tilted or curved lines, or steps between parts of the slice: the
  geometry (arc/bars solution) is wrong for this slice.
Your options:
  <cr>      show the next slice.
  q / quit  stop asking; the remaining slices are still drawn, each for
            plot_pause seconds.
  Display only: no answer changes the reduction.
Default:
  <cr> shows the next slice. Below plot_level 3 these are not shown.
```

## SubtractScatteredLight

### P50 — scattered light

> ⚠ **Check first:** What the ±180 px columns used for scattered light sample on the detector.

- Site: `SubtractScatteredLight.py:107` (sleep :109)
- Helper: `pause(...)`
- Figure: F47 (pl >= 1, not for nod-and-shuffle frames or with `skipscat`; saved as `scat_<frame>_....png` after the prompt)
- Changes reduction: No
- Status:
- Reviewer note: I describe the measured region only as "the central columns"; please add the physical reason (which part of the detector these columns sample) if wanted.

**what**
```text
Scattered-light profile along CCD y and the smooth model subtracted from every column.
```

**help**
```text
What you see:
  Title "SCATTERED LIGHT FIT": for every CCD row, the median of a band
  of columns at the centre of the CCD (+- 180 unbinned px), which runs
  along the spectral (y) direction.
  y axis in e-: CorrectGain has already converted counts to electrons
  (it runs just before this step).
  Blue circles ("Scat"): the medians.
  Red thick line ("Smooth"): Savitzky-Golay smoothing (window 151, 303
    or 501 rows, wider for noisier data); this curve is subtracted
    from every column.
What to look for:
  Good: the red curve follows the slow trend of the points without
  following individual bumps.
  Bad: the curve chasing features, or a level that looks like real
  signal rather than a smooth background: it would be subtracted from
  the science data.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

## MakeMasterFlat

### P21 — wavelength slope fit (internal flats)

- Site: `MakeMasterFlat.py:270` (sleep :272)
- Helper: `pause(...)`
- Figure: F21 (pl >= 1, internal flats (SFLAT) only)
- Changes reduction: No
- Status:

**what**
```text
Internal flat: lamp spectrum slope at the centre of the wavelength range, reference slice.
```

**help**
```text
What you see:
  Title "WAVE SLOPE FIT": the un-vignetted part of the reference slice
  in a narrow wavelength window around the centre of the flat's
  wavelength range (cut at 5620 A on blue, 5580 A on red); the window
  is 1/15 of that range wide. x: wavelength (A).
  The reference slice is slice 9 in all the flat plots; it is not the
  reference bar of the arc plots.
  Blue circles: counts.
  Red: quadratic fit, divided out before the vignetting is measured.
What to look for:
  Good: the red line follows the trend of the points.
  Bad: a fit pulled away by outliers, which then tilts the vignetting
  correction in the next plots.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P22 — vignetting

- Site: `MakeMasterFlat.py:331` (sleep :333)
- Helper: `pause(...)`
- Figure: F22 (pl >= 1, internal flats only)
- Changes reduction: No
- Status:

**what**
```text
Internal flat: vignetting along the slice edge and its linear correction.
```

**help**
```text
What you see:
  Title "Vignetting": slope-corrected flat ratio against position along
  the reference slice (px); y range fixed at 0.5 - 1.25.
  Internal flats only: their lamp light falls off near one edge of
  every slice (left on blue, right on red). The DRP fits a line to the
  fall-off and one to the flat part; where they cross is the knee.
  Blue circles ("Data"): data.
  Red dashed ("Vign."): linear fit to the vignetted part.
  Purple dashed verticals ("Vign. window"): the window of that fit.
  Red line ("UnVign."): linear fit to the un-vignetted part, drawn
    over its window and out to the vignetted edge.
  Purple verticals ("UnVign. window"): the window of that fit.
  Magenta vertical ("Crossing"): where the two fits cross (the knee).
  Purple dotted verticals ("Buffer"): buffer zone either side of the
    crossing (next plot).
What to look for:
  Good: each fit follows the data in its own region and they cross at
  the knee where the fall-off begins.
  Bad: fits that miss the data, or a crossing far from the knee: the
  internal flat will over- or under-correct that slice edge.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P23 — vignetting buffer region

- Site: `MakeMasterFlat.py:375` (sleep :377)
- Helper: `pause(...)`
- Figure: F23 (pl >= 1, internal flats only)
- Changes reduction: No
- Status:

**what**
```text
Internal flat: smooth fit across the buffer zone between corrected and uncorrected parts.
```

**help**
```text
What you see:
  Title "Buffer Region": the reference-slice ratio data within the
  buffer zone (+- 6 unbinned px) around the knee of P22, against
  slice position.
  Why: on the vignetted side of the buffer, the flat is scaled up to
  the red un-vignetted line of P22, which removes the lamp's fall-off.
  A sharp switch at the knee would leave a step, so inside the buffer
  the flat is instead scaled by that line divided by this cubic fit.
  Blue circles ("Data"): data.
  Red ("Fit"): cubic fit used to join the corrected and uncorrected
    parts smoothly.
What to look for:
  Good: a smooth red curve through the points.
  Bad: a curve that swings away from the data at either end, which
  leaves a step in the flat at the edge of the buffer.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P24 — BM ledge

> ⚠ **Check first:** The description of the BM ledge and of a good peak.

- Site: `MakeMasterFlat.py:452` (sleep :454)
- Helper: `pause(...)`
- Figure: F24 (BM grating, not a twilight flat, ledge inside WAVGOOD; **no `plot_level` guard**: drawn and paused at every level)
- Changes reduction: No
- Status:
- Reviewer note: the figure and its pause run at plot_level 0 too (inventory §1 note; with no viewer this is the B9 path fixed in T.11). `pause` keeps that: it sleeps below plot_level 2.

**what**
```text
BM grating flat: the 'ledge' step near its predicted wavelength, reference slice.
```

**help**
```text
What you see:
  Title "Ledge": reference-slice flat within +- 25 A of the wavelength
  where the BM grating's ledge (a step in the flat-field spectrum) is
  predicted from CWAVE.
  Blue circles: smoothed data.
  Red line: spline fit.
  Black dashed vertical: predicted ledge wavelength.
What to look for:
  Good: a clear step near the black line. Its position is found from the
  derivative in the next plot.
  Bad: no step, or several similar steps; then the peak choice in the
  next step needs care.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. Below plot_level 2 the plot is shown for plot_pause
  seconds (also at plot_level 0).
```

### P25 — BM ledge peak choice

> ⚠ **Check first:** The description of the BM ledge and of a good peak.

- Site: `MakeMasterFlat.py:476` and `:478` (one logical prompt; inventory :471/:473); the line `Please indicate the integer pixel value of the peak` is printed first (:475)
- Helper: `ask_value(self.logger, prompt, help, parse, default=<automatic index>, effect, validate, who="MakeMasterFlat")`, replacing the `isnumeric()` loop; the automatic index is `peaks[properties['peak_heights'].argmax()]`.
- Figure: F25 "Deriv" (drawn at every plot_level)
- Changes reduction: **Yes**: the ledge position `apk`, hence the step ratio and the correction of the reference-slice data at and redward of it, the master illumination fit, and the master flat (`*_mflat.fits` for internal, `*_mdome.fits` for dome flats; twilight flats skip this block). The ratio is logged ("BM ledge ratio"); no header keyword records the choice.
- Prompt: `Peak? <int> (<cr> = highest peak, index %d)`
- Accepted input: a non-negative integer, digits only (`str.isnumeric()` and `int()`), from 0 to 59 (`len(xvals) - 1`; the derivative has 60 points). Out of range (today `IndexError`, B1c) and characters such as `²` (today `ValueError`) re-prompt in PR3a. Empty: today re-asks; in PR3a it takes the default (allowed by D12, since it is invalid today).
- Effect: `BM ledge peak set by user to index {value} (automatic choice: index %d)`
- Status:
- Reviewer notes: (1) Since the B1/B1b hotfix (D23, PR #223) the prompt is reached only with **two or more** peaks (zero peaks: warning, no ledge correction; one peak: used). (2) The automatic choice now indexes correctly (`peaks[argmax]`), but it is still just the highest derivative peak; the code itself warns "This may not be the correct peak!". The help says so. (3) The inventory's "B1/B1b apply to the non-interactive branch" no longer holds at this commit.

**what**
```text
Several BM ledge peaks found: enter the index (x value) of the peak that is the ledge.
```

**help**
```text
What you see:
  Title "Deriv": derivative of the ledge fit, in 0.5 A steps, against
  index 0-59 (axis "px (Wavelength)"; index i is at predicted ledge
  wavelength - 15 A + 0.5 i A).
  Blue circles: derivative.
  Black dashed vertical (index 30): the predicted ledge wavelength.
  Magenta verticals: peaks found (height >= 20); there is more than
    one.
What to look for:
  The ledge is a single sharp peak, usually the strongest one close to
  the black line. Smaller or broader peaks come from noise or spectral
  structure in the flat.
Your options:
  <index>  put the ledge at that index (read it off the x axis). The
           step size is measured from linear fits 3-8 A either side,
           and the reference-slice flat redward of it is divided by it.
           This changes the master flat (*_mflat or *_mdome) and every
           frame flat-fielded with it.
  <cr>     use the automatic choice, the highest peak.
Default:
  <cr> takes the highest peak, as the reduction does below plot_level
  2. The highest peak is not always the ledge: check the next plots
  (BM Ledge Region) and rerun at plot_level 2 if the step remains.
```

### P26 — peak of ledge

- Site: `MakeMasterFlat.py:500` (inventory :494)
- Helper: `pause(...)` inside the existing `plot_level >= 3` block (always prompts there).
- Figure: F26 (pl >= 3)
- Changes reduction: No
- Status:
- Reviewer note: the sleep at :502 is dead code (inventory §4 row 16).

**what**
```text
BM ledge: derivative against wavelength with the chosen ledge position.
```

**help**
```text
What you see:
  Title "Peak of ledge": the derivative of the ledge fit against
  wavelength (A).
  Blue circles: data.
  Magenta vertical ("Pk"): the ledge position used.
What to look for:
  Good: the magenta line sits on the sharp peak you intended.
  Bad: it sits on a noise peak; the correction will put a step into
  the flat instead of removing one.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. Below plot_level 3 this plot is not shown.
```

### P27 — BM ledge region

- Site: `MakeMasterFlat.py:557` (inventory :551; sleep :559)
- Helper: `pause(...)`
- Figure: F27 (pl >= 1)
- Changes reduction: No
- Status:

**what**
```text
BM ledge: flat data before and after the step correction.
```

**help**
```text
What you see:
  Title "BM Ledge Region", within +- 25 A of the predicted ledge.
  Blue circles: reference-slice flat data.
  Cyan circles: the same data redward of the ledge after dividing by
    the step ratio ("Corrected").
  Red: spline fit to the smoothed data.
  Red dashed: blueward linear fit.
  Purple dashed verticals: the window of the blueward fit.
  Red dotted: redward linear fit.
  Purple dotted verticals: the window of the redward fit.
  Magenta vertical ("Pk"): ledge position.
What to look for:
  Good: the cyan points continue the level of the blue points below
  the ledge, with no step left.
  Bad: a step remains or is overcorrected: the wrong peak was taken.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P28 — blue/red slice fits

- Site: `MakeMasterFlat.py:664` (inventory :658; sleep :666)
- Helper: `pause(...)`
- Figure: F28 (pl >= 1; saved as `redblue_<frame>_....png` after the prompt)
- Changes reduction: No
- Status:

**what**
```text
Flat spectra of the reference slice and of the slices reaching furthest blue and red.
```

**help**
```text
What you see:
  Title "Blue/Red fits", flux (e-) against wavelength. Each slice is
  drawn as data (dots) and spline fit (line) in one colour.
  Black: reference slice (9).
  Blue: the slice reaching furthest to the blue (12 for the blue
    camera, 11 for red).
  Cyan: the slice reaching furthest to the red (23 for the blue
    camera, 0 for red).
  Green verticals: WAVGOOD.
  The end slices extend the illumination spectrum beyond the
  reference slice's wavelength range.
What to look for:
  Good: each line follows its dots; the curves have the same shape
  where they overlap.
  Bad: a fit that departs from its data, or very different shapes
  in the overlap.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P29 — blue end fits

- Site: `MakeMasterFlat.py:733` (inventory :727; sleep :735)
- Helper: `pause(...)`
- Figure: F29 (pl >= 1 and more than one point in the blue overlap)
- Changes reduction: No
- Status:

**what**
```text
Blue overlap: reference-slice and blue-slice flat fits.
```

**help**
```text
What you see:
  Title "Blue fits": the spline fits in the overlap near the blue end
  of the reference slice (5-10 % into its range).
  Black line and circles: reference slice.
  Blue line and circles: the blue-end slice.
What to look for:
  Good: both curves smooth, with similar shapes; their ratio (next
  plot) scales the blue slice to join the reference.
  Bad: noisy or diverging curves: the blue extension of the flat will
  have a step or slope error.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P30 — blue end ratio

- Site: `MakeMasterFlat.py:753` (inventory :747; sleep :755)
- Helper: `pause(...)`
- Figure: F30
- Changes reduction: No
- Status:

**what**
```text
Blue overlap: ratio of reference to blue-slice fit, and its linear fit.
```

**help**
```text
What you see:
  Title "Blue ratios".
  Black line and circles: reference-slice fit divided by blue-slice
    fit.
  Red line and circles: linear fit to that ratio, used to scale the
    blue-slice data blueward of the reference slice.
What to look for:
  Good: a smooth, nearly constant ratio, well followed by the red
  line.
  Bad: a strongly curved or noisy ratio; the linear scaling will not
  match it. (If either fit crosses zero only an offset is applied.)
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P31 — red end fits

- Site: `MakeMasterFlat.py:772` (inventory :766; sleep :774)
- Helper: `pause(...)`
- Figure: F31 (pl >= 1 and more than one point in the red overlap)
- Changes reduction: No
- Status:

**what**
```text
Red overlap: reference-slice and red-slice flat fits.
```

**help**
```text
What you see:
  Title "Red fits": the spline fits in the overlap near the red end of
  the reference slice (90-95 % into its range).
  Black line and circles: reference slice.
  Cyan line and circles: the red-end slice.
What to look for:
  Good: both curves smooth, with similar shapes; their ratio (next
  plot) scales the red slice to join the reference.
  Bad: noisy or diverging curves: the red extension of the flat will
  have a step or slope error.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P32 — red end ratio

- Site: `MakeMasterFlat.py:790` (inventory :784; sleep :792)
- Helper: `pause(...)`
- Figure: F32
- Changes reduction: No
- Status:

**what**
```text
Red overlap: ratio of reference to red-slice fit, and its linear fit.
```

**help**
```text
What you see:
  Title "Red ratios".
  Black line and circles: reference-slice fit divided by red-slice
    fit.
  Red line and circles: linear fit to that ratio, used to scale the
    red-slice data redward of the reference slice.
What to look for:
  Good: a smooth, nearly constant ratio, well followed by the red
  line.
  Bad: a strongly curved or noisy ratio; the linear scaling will not
  match it. (If either fit crosses zero only an offset is applied.)
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P33 — master illumination

- Site: `MakeMasterFlat.py:894` (inventory :888; sleep :896)
- Helper: `pause(...)`
- Figure: F33 (pl >= 1; saved as `flat_<frame>_....png` after the prompt)
- Changes reduction: No
- Status:

**what**
```text
Master illumination spectrum: combined flat data over the full range, with fit.
```

**help**
```text
What you see:
  Title "Master Illumination", flux (e-) against wavelength.
  Blue dots: reference-slice data plus the scaled blue- and red-end
    slice data.
  Red thick line: the spline fit, the illumination spectrum the flat
    is divided by.
  Purple verticals ("Cor lim"): the reference slice's range; outside
    it the end slices are used.
  Green verticals: WAVGOOD.
What to look for:
  Good: the red curve follows the data over the whole range, with no
  jump at the purple lines and no ringing at the ends.
  Bad: a jump or kink at a purple line (end-slice scaling), or wiggles
  where the data are smooth: these go straight into the master flat.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

## MakeMasterSky

### P34 — standard star mask

- Site: `MakeMasterSky.py:347` (sleep :349)
- Helper: `pause(...)`
- Figure: F34 (pl >= 1, standard-star frames)
- Changes reduction: No
- Status:

**what**
```text
Standard star: profile in its brightest slice and the region masked from the sky model.
```

**help**
```text
What you see:
  Title "Std max sl <n>": flux against position along slice n, the
  slice with the largest spread, in a 10 % wavelength window around
  the middle wavelength.
  Blue dots: data.
  Red: Gaussian fit.
  Magenta vertical: peak position.
  Purple verticals: +- 5 sigma of the fit; this range of slice
    positions is excluded from the sky model in all slices.
What to look for:
  Good: the Gaussian fits the star and the purple lines enclose its
  wings.
  Bad: a fit to the wrong feature, or a mask much too narrow (star
  light enters the sky model and is subtracted) or too wide.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P35 — continuum source mask

- Site: `MakeMasterSky.py:435` (sleep :437)
- Helper: `pause(...)`
- Figure: F35 (pl >= 1, `kcwi.sky` "cont" entry without a source position)
- Changes reduction: No
- Status:
- Reviewer note: reached only through a `kcwi.sky` `cont` entry; because of B36 (PR3b), `kcwi.sky` is ignored when the DRP is run with paths, so this prompt is reachable only when it is run in the raw-data directory with bare file names.

**what**
```text
Continuum source: profile in its brightest slice, source region and local sky windows.
```

**help**
```text
What you see:
  Title "Cont source max sl <n>": flux against position along slice n,
  in a 10 % wavelength window around the middle wavelength.
  Blue dots: data.
  Red: Gaussian fit.
  Magenta vertical: peak.
  Purple verticals: +- 7 sigma, the source region.
  Purple dashed verticals: 14 unbinned px further out.
  Only the two windows between the solid and dashed purple lines are
  used for the sky model.
What to look for:
  Good: the solid purple lines enclose the source wings; the windows
  contain only sky.
  Bad: a fit to the wrong feature, or source light inside a window.
  Give the position and width in kcwi.sky to override the fit.
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

### P36 — master sky

- Site: `MakeMasterSky.py:607` (sleep :609)
- Helper: `pause(...)`
- Figure: F36 (pl >= 1; saved as `sky_<frame>_....png` after the prompt)
- Changes reduction: No
- Status:

**what**
```text
Sky model: sky pixels from all slices against wavelength, with the spline fit.
```

**help**
```text
What you see:
  Title "Master Sky", flux (e-) against wavelength.
  Blue dots: unmasked sky pixels from all slices (subsampled).
  Red: the spline sky model subtracted from this frame.
  Green verticals: WAVGOOD.
What to look for:
  Good: the red curve runs through the dots, following sky lines and
  continuum; the dots scatter tightly around it.
  Bad: many dots well above the curve (object light in the sky region:
  mask it with a sky mask or a kcwi.sky entry), or line peaks the model
  misses (over/under-subtracted sky lines in the cube).
Your options:
  <cr>  continue. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for plot_pause
  seconds; at plot_level 0 it is not drawn.
```

## MakeInvsens

### P13 — standard star per slice

- Site: `MakeInvsens.py:196`
- Helper: `ask_next(self.logger, what, help, who="MakeInvsens")`; `QUIT` sets `do_plots = False`.
- Figure: F13 (pl >= 3; one per slice; the line `sl, std: ...` is printed first)
- Changes reduction: No
- Status:

**what**
```text
Standard-star cube: flux summed over wavelength along this slice.
```

**help**
```text
What you see:
  Title "Std Slice <n> (DAR padded)".
  Blue x: flux summed over the good wavelength range, against
    position along slice n.
  The line printed above gives the slice and the standard deviation of
  this profile; the slice with the largest value is taken as the
  star's.
What to look for:
  Good: the star is a clear peak in one slice (fainter in its
  neighbours); other slices are flat sky.
  Bad: no clear peak anywhere, or the largest spread in a slice without
  the star (e.g. an edge artefact): the wrong slices would be summed.
Your options:
  <cr>      show the next slice.
  q / quit  stop showing slices.
  Display only: no answer changes the reduction.
Default:
  <cr> shows the next slice. Below plot_level 3 these are not shown.
```

### P14 — fit wavelength limits

> ⚠ **Check first:** The features suggested for masking (Balmer lines, O2 B/A bands at 6870/7600 Å).

- Site: `MakeInvsens.py:364` (inventory :365)
- Helper: `ask_value(self.logger, prompt, help, parse, default=None, effect, validate, who="MakeInvsens")` in the existing `while not done` loop: `None` (empty answer) ends the loop; a pair sets `wlm0, wlm1` and the loop redraws.
- Figure: F14 (pl >= 2)
- Changes reduction: **Yes**: the wavelength range of the inverse-sensitivity and effective-area fits, so `*_invsens.fits` and `*_ea.fits` (header `INVFW0`, `INVFW1`), and every object later flux-calibrated with them.
- Prompt: `New fit limits? <float> <float> (<cr> = accept)`
- Accepted input (D12): two numbers separated by spaces, each read by Python `float()` (e.g. `3700 5500`, `3.7e3 5500.5`); words after the second number are ignored. Valid if both are finite, `lo <= hi`, `lo >= WAVALL0` and `hi <= WAVALL1` (`nan` is refused from PR3a: today it passes the check, leaves no wavelength in the fit range and MakeInvsens fails; maintainer decision 2026-10-28). One number, a non-number, or failing the check: the reason is printed and the previous limits are kept (B6, fixed in PR3a; today both limits are reset to WAVGOOD).
- Effect: `fit wavelength limits set by user to {value[0]:.1f} - {value[1]:.1f} A (WAVGOOD %.1f - %.1f A)`
- Status:
- Reviewer notes: (1) `nan` passed today's check (every comparison with `nan` is false), emptied the fit range and made MakeInvsens fail; the maintainer classed it as invalid input (2026-10-28), so PR3a refuses it and asks again. (2) The `print` lines before the loop (L337-342) say "Hover the cursor over spectrum to find wavelengths"; after PR1 tooltips are ignored, but the viewer's toolbar shows the cursor position (to confirm in P.6). P.4 should replace those prints with `what`/`help`.

**what**
```text
Check the wavelength limits of the sensitivity fit; enter new limits to avoid bad edges.
```

**help**
```text
What you see:
  Title "<std> Obs Spec". The x range is WAVALL0 - WAVALL1.
  Blue: the standard star's observed spectrum (e-/s, summed over the
    star's slices, extinction corrected).
  Green verticals: WAVGOOD, the range covered by all slices.
  Purple verticals: the current fit limits ("LIMITS").
  Orange dash-dot-dot: CWAVE.
  The viewer shows the cursor position, to read off wavelengths.
What to look for:
  Good: between the purple lines the spectrum is smooth and well
  exposed.
  Bad: the range includes a steep edge, noisy ends, or an artefact
  (e.g. near the dichroic or the end of the detector); a polynomial
  fitted over it will be distorted everywhere. Move the limits inward.
Your options:
  <cr>     accept the limits shown and continue.
  lo hi    new limits in Angstrom, e.g. "3700 5500"; must satisfy
           lo <= hi within WAVALL0 - WAVALL1. The plot is redrawn and
           you are asked again.
  The limits set the fit range of the inverse sensitivity and
  effective area (*_invsens.fits, *_ea.fits, INVFW0/INVFW1).
Default:
  <cr> keeps the current limits: WAVGOOD0 - WAVGOOD1 (WAVGOOD0 raised
  to 3650 A if lower) unless you changed them. Below plot_level 2 the
  WAVGOOD limits are used.
```

### P15 — line masks

> ⚠ **Check first:** The features suggested for masking (Balmer lines, O2 B/A bands at 6870/7600 Å).

- Site: `MakeInvsens.py:468` (inventory :469)
- Helper: `ask_value(self.logger, prompt, help, parse, default=None, effect, validate, who="MakeInvsens")` in the existing loop; `None` ends it; a value `(w0, w1, comment)` is appended to `lmasks`.
- Figure: F15 (pl >= 2)
- Changes reduction: **Yes**: masked ranges get (near) zero weight in the fits, changing `*_invsens.fits` and `*_ea.fits`; all masks are written to `<star>.lmsk` in the current working directory, which later runs read back.
- Prompt: `Mask line: <start> <stop> [comment] in A (<cr> = done)`
- Accepted input (D12): two numbers read by `float()`, then optional comment words. Valid if both lie strictly inside the current fit limits (`wlm0 < w < wlm1`); otherwise "line mask outside wl range" and nothing is added. `start > stop` is accepted (and masks nothing). One number only: re-prompts (B2, fixed in PR3a; today it crashes MakeInvsens).
- Effect: `line mask added by user: {value[0]:.2f} - {value[1]:.2f} A, saved to %s`
- Status:
- Reviewer notes: (1) the existing `print` at L438 says masks are "vertical dashed yellow lines"; the code draws them orange (L462, L464). The help says orange; P.4 should fix or drop that print. (2) B13 (PR3b): header lines of an existing `.lmsk` are written back without a newline, so on the next run the first mask is lost. (3) B32 (PR3b): the shipped `.lmsk` files in `data/stds/` are never found; only a `.lmsk` in the working directory is read. The help therefore says "from a local .lmsk file".

**what**
```text
Mask narrow features (absorption lines, telluric bands, cosmic rays) out of the fit.
```

**help**
```text
What you see:
  Title "<std> Obs Spec", as before.
  Blue: the observed standard spectrum.
  Green verticals: WAVGOOD.
  Purple verticals: the fit limits.
  Orange dash-dot-dot: CWAVE.
  Purple dashed pairs: line masks already defined (from a local .lmsk
    file or entered here).
What to look for:
  Narrow features that the reference spectrum does not reproduce at
  KCWI's resolution: stellar absorption lines (e.g. Balmer lines in hot
  standards), telluric bands (e.g. O2 B and A at 6870 and 7600 A),
  residual cosmic rays. Unmasked, they bend the smooth fit. Mask only
  the feature itself, not broad stretches of continuum.
Your options:
  <cr>               done; use the masks shown.
  start stop [text]  add a mask from start to stop (Angstrom, start <
                     stop), optional comment; both ends must be inside
                     the purple limits. The plot is redrawn with it.
  Masks cannot be removed here: edit or delete <star>.lmsk (in the
  directory you started the DRP from) and rerun. Masked points are left
  out of the fits for *_invsens.fits and *_ea.fits.
Default:
  <cr> keeps the masks shown. Below plot_level 2 only masks from the
  local .lmsk file are used (none if there is no file).
```

### P16 — efficiency

- Site: `MakeInvsens.py:584` (inventory :585; sleep :586 is `2. * plot_pause`)
- Helper: `pause(...)`
- Figure: F16 (pl >= 1; saved as `std_..._eff.png` after the fit loop)
- Changes reduction: No
- Status:
- Reviewer note: the sleep here and at P17, P18, P19 is `2 * plot_pause`; `pause` sleeps `config.plot_pause`. P.3 must keep the doubled pause (e.g. keep the existing branch until P.1 gives `pause` a factor), or the timing changes at plot_level 1.

**what**
```text
Effective efficiency of telescope + KCWI (+ atmosphere) from the standard, with fit.
```

**help**
```text
What you see:
  Title "<std> Efficiency": observed photons / expected photons on the
  Keck aperture, in per cent (not extinction corrected).
  Blue: data.
  Red: polynomial fit.
  Grey x: masked points ("Rej").
  Purple verticals: fit limits.
  Green verticals: WAVGOOD.
  Red dashed horizontal: peak of the fit.
  Red dash-dot horizontal: mean of the fit.
What to look for:
  Good: a smooth throughput curve; the red fit follows it without
  wiggles, also near the purple limits.
  Bad: the fit oscillating at the ends (order too high, or limits too
  wide), or bumps at absorption lines that were not masked.
Your options:
  <cr>  continue to the next plot. The fit order is chosen after the
        residuals plot. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for twice
  plot_pause seconds; at plot_level 0 it is not drawn.
```

### P17 — inverse sensitivity

- Site: `MakeInvsens.py:605` (inventory :606; sleep `2. * plot_pause`)
- Helper: `pause(...)` (see P16 note)
- Figure: F17 (pl >= 1; log y; saved as `std_..._invsens.png`)
- Changes reduction: No
- Status:

**what**
```text
Inverse sensitivity (flux per e-/s) from the standard star, with its polynomial fit.
```

**help**
```text
What you see:
  Title "<std> Inverse sensitivity", logarithmic y axis: reference flux
  divided by observed e-/s, the curve that flux-calibrates later data.
  Blue: data.
  Red: polynomial fit (also drawn outside the limits, where it is
    extrapolated).
  Grey x: masked points.
  Purple verticals: fit limits.
  Green verticals: WAVGOOD.
What to look for:
  Good: the red fit follows the blue data closely between the limits.
  Bad: systematic departures, or the fit turning sharply up or down
  just outside the limits where you still need calibrated data.
Your options:
  <cr>  continue to the next plot. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for twice
  plot_pause seconds; at plot_level 0 it is not drawn.
```

### P18 — calibrated standard

- Site: `MakeInvsens.py:624` (inventory :625; sleep `2. * plot_pause`)
- Helper: `pause(...)` (see P16 note)
- Figure: F18 (pl >= 1; saved as `std_..._cal.png`)
- Changes reduction: No
- Status:

**what**
```text
Standard star calibrated with the fitted inverse sensitivity, over its reference spectrum.
```

**help**
```text
What you see:
  Title "<std> Calibrated", flux in erg/s/cm^2/A.
  Blue: observed spectrum times the fitted inverse sensitivity.
  Black: the star's reference spectrum.
  Purple verticals: fit limits.
  Green verticals: WAVGOOD.
What to look for:
  Good: blue lies on black between the limits, apart from narrow
  features (absorption lines, telluric bands) that differ in resolution.
  Bad: broad differences (fit order or limits) or a mismatch at a
  feature you meant to mask.
Your options:
  <cr>  continue to the residuals. Display only.
Default:
  <cr> continues. At plot_level 1 the plot is shown for twice
  plot_pause seconds; at plot_level 0 it is not drawn.
```

### P19 — residuals and fit order

- Site: `MakeInvsens.py:654` (inventory :655); at pl 1 the `else` sets `done = True` and sleeps `2. * plot_pause` (:685)
- Helper: `ask_value(self.logger, prompt, help, parse=int, default=None, effect, who="MakeInvsens")` inside the existing `plot_level >= 2` branch; `None` ends the loop; an int refits as today.
- Figure: F19 (pl >= 1; F16–F18 were shown just before; saved as `std_..._resid.png`)
- Changes reduction: **Yes**: the polynomial order of the inverse-sensitivity and effective-area fits: `*_invsens.fits`, `*_ea.fits`, header `INVFORD`.
- Prompt: `Current fit order = %d. New fit order? <int> (<cr> = accept)`
- Accepted input (D12): anything Python `int()` reads (e.g. `5`, `+5`, ` 5 `). Text or decimals: "Bad fit order, try again". No bounds: a negative order makes the fit fail with "Bad fit order"; very large orders are accepted.
- Effect: `inverse sensitivity refitted by user with polynomial order {value} (default %d)`
- Status:
- Reviewer note: B3 (PR3b): `ford` is set before the fit, so after a failed fit (e.g. order -1) the recorded order (`INVFORD`, and this log line if logged on acceptance) is the bad value while the fit is the previous one. P.4 should log the effect only after the refit succeeds, or accept that the log line repeats today's mismatch until PR3b.

**what**
```text
Calibration residuals; accept the fit or enter a new polynomial order to refit.
```

**help**
```text
What you see:
  Title "<std> Residuals = <mean> +- <std> (%)": (calibrated -
  reference) / reference in per cent, over the fit range.
  Blue: residuals.
  Grey x: masked points (if any).
  Purple verticals: fit limits.
  Green verticals: WAVGOOD.
  Red: mean.
  Red dashed: mean +- 1 standard deviation.
What to look for:
  Good: residuals scatter around zero with no large-scale wave.
  Bad: a slow wave across the range (order too low), or oscillations
  growing towards the limits (order too high, or limits too wide).
Your options:
  <cr>      accept the current fit order and finish.
  <number>  refit with that polynomial order; the efficiency, inverse
            sensitivity, calibrated and residual plots are shown again
            and you are asked again.
  The order goes into *_invsens.fits and *_ea.fits (keyword INVFORD).
Default:
  <cr> keeps the current order: 9 (7 for BL) times the dichroic
  fraction, at least 3, unless you changed it. Below plot_level 2 the
  default order is used without review.
```

## NandshuffSubtractSky

### P37 — sky scale, blue camera

> ⚠ **Check first:** Why a sky scale is needed in aborted nod-and-shuffle (the agent's own reading).

- Site: `NandshuffSubtractSky.py:185`
- Helper: `ask_value(self.logger, prompt, help, parse=float, default=1.0, effect, validate, who="NandshuffSubtractSky")`, replacing the `while skyscl < 0.` loop (inside its `plot_level >= 2` test).
- Figure: none (this primitive draws nothing; the last figure stays up)
- Changes reduction: **Yes**: the sky panel is multiplied by the factor before it is subtracted: `*_intk.fits`, `*_sky.fits`, header `NASSCL`.
- Prompt: `Enter scale factor for sky to match obj <float> (<cr> = 1.0)`
- Accepted input (D12): anything `float()` reads that is not negative: `0`, `1`, `0.95`, `1e0`, also `nan` and `inf` (B7, PR3b). Negative or non-numeric: warning and asked again. Empty: today re-asks; in PR3a it takes the default 1.0.
- Effect: `nod-and-shuffle sky scale set by user to {value} (default 1.0)`
- Status:
- Reviewer note: reached only for blue frames whose sky panel does not start at the bottom (`skyrow0 >= 10`, "Inverted N&S obs detected", header comment "Aborted nod-and-shuffle observations"). No real data of this kind exists for testing (Q4). The scientific reason for a scale (unequal sky exposure in an aborted script) is my reading; please confirm.

**what**
```text
Aborted (inverted) nod-and-shuffle frame: enter the scale for the sky panel.
```

**help**
```text
What you see:
  No new plot. This blue-camera frame has its nod-and-shuffle panels in
  the inverted order left by an aborted N&S script, so the sky panel
  need not have the same sky exposure as the object panel.
What to look for:
  Estimate the ratio of the sky level in the object panel to that in
  the sky panel, e.g. from a bright sky line in the *_int.fits image
  (ds9). 1.0 means equal sky in both panels.
Your options:
  <number>  multiply the sky panel by it, then subtract it from the
            object panel. Changes *_intk.fits and *_sky.fits; the value
            is recorded as NASSCL. 0 means no sky subtraction.
  <cr>      use 1.0.
  Negative numbers and text are refused.
Default:
  <cr> uses 1.0, as the reduction does below plot_level 2.
```

### P38 — sky scale, red camera

> ⚠ **Check first:** Why a sky scale is needed in aborted nod-and-shuffle (the agent's own reading).

- Site: `NandshuffSubtractSky.py:269`
- Helper: as P37.
- Figure: none
- Changes reduction: **header only today**: in the red non-standard branch the factor is written to `NASSCL` but `sky` is never multiplied by it (compare L205 in the blue branch with L281-286), and the branch then fails at L303/L305 (B12), so no outputs are written.
- Prompt / Accepted input / Effect: as P37.
- Status:
- Reviewer notes: (1) **Code and inventory disagree**: the inventory says "Yes: as P37"; the code does not apply the scale in this branch. This is not in §13; suggest adding it to B12's PR3b fix (apply `sky *= skyscl` and fix the indices). (2) Reached only for red frames with `objrow0 >= 10`, i.e. the inverted-panel case (`NSHFDN != NSHFUP`, L103-111, object panel in the middle). No real data of this kind exists for testing (Q4). (3) Once the fix lands, the help text below must lose its last "Note" paragraph.

**what**
```text
Aborted (inverted) nod-and-shuffle frame: enter the scale for the sky panel.
```

**help**
```text
What you see:
  No new plot. This red-camera frame has its nod-and-shuffle panels in
  a non-standard order, so the sky panel need not have the same sky
  exposure as the object panel.
What to look for:
  Estimate the ratio of the sky level in the object panel to that in
  the sky panel, e.g. from a bright sky line in the *_int.fits image
  (ds9). 1.0 means equal sky in both panels.
Your options:
  <number>  scale factor for the sky panel, recorded as NASSCL.
  <cr>      use 1.0.
  Negative numbers and text are refused.
  Note: in this version the red-channel branch does not apply the
  factor and stops with an error after this prompt (known bug), so no
  *_intk.fits is written for this frame.
Default:
  <cr> uses 1.0, as the reduction does below plot_level 2.
```

## Counts

| File | Prompts | IDs |
|---|---|---|
| scripts/reduce_kcwi.py | 4 | S01–S04 |
| scripts/kcwi_rti.py | 4 | S05–S08 |
| **Scripts** | **8** | |
| SubtractOverscan.py | 2 | P48–P49 |
| MakeMasterBias.py | 1 | P20 |
| FindBars.py | 3 | P05–P07 |
| TraceBars.py | 1 | P51 |
| ExtractArcs.py | 2 | P03–P04 |
| ArcOffsets.py | 2 | P01–P02 |
| ReadAtlas.py | 2 | P39–P40 |
| FitCenter.py | 3 | P08–P10 |
| GetAtlasLines.py | 1 | P11 |
| SolveArcs.py | 6 | P42–P47 |
| SolveAIT.py | 1 | P41 |
| MakeCube.py | 1 | P12 |
| SubtractScatteredLight.py | 1 | P50 |
| MakeMasterFlat.py | 13 | P21–P33 (P25 = 2 `input(` lines) |
| MakeMasterSky.py | 3 | P34–P36 |
| MakeInvsens.py | 7 | P13–P19 |
| NandshuffSubtractSky.py | 2 | P37–P38 |
| **Primitives (17 files)** | **51** | from 52 `input(` lines |
| **Total** | **59** | matches inventory §5 |

By helper: `pause` 35 (P02, P05, P06, P09–P11, P16–P18, P20–P24,
P26–P36, P39, P41, P43, P45–P51), `ask_next` 9 (P01, P03, P04, P07, P08,
P12, P13, P42, P44; P01, P08 and P44 with `n_items`), `ask_value` 7 (P14,
P15, P19, P25, P37, P38, P40): 35 + 9 + 7 = 51. `ask_yes_quit` 8
(S01–S08).

Changes the reduction: P14, P15, P19, P25, P37, P38 (header only, see its
note), P40, and S01–S04, S07, S08 (S05, S06 not today, B11). Matches
inventory §2.7 except P38 and S05/S06, which the notes explain.

## Colours named per site

Every colour named in a help text, in the [scheme](#plot-colour-scheme)
above (files in `kcwidrp/primitives/` unless stated). The plots draw the
"Colour in help" column. The "Code" column gives the line that drew each
element at `cb7444c` and the colour it had then, before the colour scheme;
the WAVGOOD lines of P16-P19, P28 and P33 and the P12 colour bar (drawn by
`_plot_render.py`, `_draw`) are new, as are the rows marked "new".

| ID | Colour in help | Code (line, colour at `cb7444c`) |
|---|---|---|
| S01–S08 | none (no figure) | — |
| P48, P49 | blue circles / red thick | SubtractOverscan.py:110 default / :111 red |
| P20 | blue bars / red thick / red dashed verticals | MakeMasterBias.py:144-145 navy / :146 `#ff8888` / :148, :150 red |
| P05 | blue line / magenta x / grey dashed | FindBars.py:89 / :90-92 red / :93-94 grey |
| P06 | blue line / magenta x / grey dashed | FindBars.py:127 / :128-130 red / :131-132 grey |
| P07 | blue line and circles / magenta vertical | FindBars.py:157 blue, :158 red / :159-160 green |
| P51 | blue x / magenta circles | TraceBars.py:139 blue / :140-141 red |
| P03 | blue / cyan / red | ExtractArcs.py:145 blue / :146 red / :147 green |
| P04 | blue / red | ExtractArcs.py:212 blue / :213 red |
| P01 | black (REFBAR) / blue (this bar) | ArcOffsets.py:86 green / :89 red |
| P02 | blue diamonds / black dash-dot | ArcOffsets.py:113 default / :117 `oplot_slices` |
| P39 | blue / magenta vertical / black dashed vertical | ReadAtlas.py:183 default / :186-187 red / :189-190 green |
| P40 | blue / black / orange dash-dot-dot | ReadAtlas.py:207 default / :211-213 red / :219 magenta dash-dot |
| P08 | blue circles and line / magenta vertical / black dashed vertical | FitCenter.py:331 red, :332 blue / :334 green / :336-337 red |
| P09 | blue x / orange dash-dot-dot / black dash-dot | FitCenter.py:359 default / :360-361 red / :365 |
| P10 | blue x / black dashed / black dash-dot | FitCenter.py:386 default / :387-388 red / :393 |
| P11 | blue / black | GetAtlasLines.py:527-528 lightgray / :529-530 blue |
| P11 | green diamonds (Kept) | GetAtlasLines.py:552-553 green |
| P11 | cyan, pink, olive, brown, grey diamonds | GetAtlasLines.py:534 cyan, :538 red, :542 orange, :546 brown, :550 magenta |
| P11 | purple verticals (WavLim) | GetAtlasLines.py:554-556 brown |
| P42 | black line and circles / black dashed vertical | SolveArcs.py:233-234 blue / :235-236 green / :237 red |
| P42 | blue circles / red line | SolveArcs.py:247 red / :245 black |
| P42 | magenta solid / dashed / dotted verticals | SolveArcs.py:256 black dash-dot / :250 green dashed / :253 magenta dash-dot |
| P43 | blue / gold diamonds / grey x | SolveArcs.py:394 green / :397 orange / :400-401 red |
| P43 | grey dotted / red dashed / orange dash-dot-dot | SolveArcs.py:404 black dotted / :405-406 gray dash-dot / :407-408 magenta |
| P43 | purple verticals (WavLim) | new |
| P44 | blue line and diamonds | SolveArcs.py:421 / :422 darkgrey |
| P44 | black / magenta diamonds / gold / grey x / orange dash-dot-dot | SolveArcs.py:425 blue / :430-431 green / :433-434 orange / :435-436 red / :427-428 magenta |
| P45 | blue diamonds / black dash-dot | SolveArcs.py:485 default / :489 |
| P46 | blue diamonds / red / red dashed / black dash-dot | SolveArcs.py:509 default / :515-516 red / :517-525 green dashed / :528 |
| P46 | black thick (slice mean) | new |
| P47 | blue diamonds / red / red dashed / black dash-dot | SolveArcs.py:552 default / :555-556 red / :557-565 green dashed / :568 |
| P41 | blue / cyan / green / grey | SolveAIT.py:160 blue / :161 red / :162-164 green / :165-167 orange |
| P12 | Spectral11, grey NaN, colour bar | MakeCube.py:276-277; `_plot_render.py:268` (`_image`, no colour bar) |
| P12 | black ticks (pixel boundaries) / green horizontals (WAVGOOD) | new |
| P50 | blue circles / red thick | SubtractScatteredLight.py:101 default / :102 red |
| P21 | blue circles / red | MakeMasterFlat.py:265 default / :266 red |
| P22 | blue circles / red dashed, red / magenta vertical | MakeMasterFlat.py:314 default / :315-316 purple, :317 red / :327 red |
| P22 | purple dashed / purple / purple dotted verticals | MakeMasterFlat.py:319-320 blue / :321-322 green / :323-326 black |
| P23 | blue circles / red | MakeMasterFlat.py:371 default / :372 red |
| P24 | blue circles / red / black dashed vertical | MakeMasterFlat.py:446 green fill / :447 default / :449 red |
| P25 | blue circles / black dashed vertical / magenta verticals | MakeMasterFlat.py:462 default / :463 green / :465 red |
| P26 | blue circles / magenta vertical | MakeMasterFlat.py:495 default / :496 red |
| P27 | blue / cyan circles / red | MakeMasterFlat.py:531 blue / :540 orange / :542 red |
| P27 | red dashed, purple dashed verticals (blueward) | MakeMasterFlat.py:547-548 purple / :543-544 blue |
| P27 | red dotted, purple dotted verticals (redward) | MakeMasterFlat.py:549-550 green / :545-546 black |
| P27 | magenta vertical (Pk) | MakeMasterFlat.py:551 green |
| P28 | black / blue / cyan (dots and lines), green verticals | MakeMasterFlat.py:656-657 black / :658-659 blue / :660-661 red |
| P29 | black / blue (line and circles) | MakeMasterFlat.py:725, :727 black / :728, :730 blue (outline default) |
| P30 | black / red (line and circles) | MakeMasterFlat.py:743, :745-746 black / :747-750 blue |
| P31 | black / cyan (line and circles) | MakeMasterFlat.py:764, :766 black / :767, :769 red |
| P32 | black / red (line and circles) | MakeMasterFlat.py:782, :784 black / :785-787 red |
| P33 | blue dots / red thick / purple verticals, green verticals | MakeMasterFlat.py:884 black / :886 red / :888, :890 orange |
| P34 | blue dots / red / magenta vertical / purple verticals | MakeMasterSky.py:338-339 purple / :341 red / :342-343 green / :344 blue |
| P35 | blue dots / red / magenta vertical / purple, purple dashed | MakeMasterSky.py:420-421 purple / :423 red / :424-431 green / :432 blue |
| P36 | blue dots / red / green verticals | MakeMasterSky.py:598 purple / :600 red / :601-602 green |
| P13 | blue x | MakeInvsens.py:193 default |
| P14 | blue / green / purple / orange dash-dot-dot | MakeInvsens.py:354 black / :355, :357 green / :358, :360 blue / :361 red |
| P15 | as P14, plus purple dashed pairs | MakeInvsens.py:452-459 as P14 / :462, :464 orange dashed |
| P16 | blue / red / grey x / purple, green verticals | MakeInvsens.py:569 black / :571 red / :573 default / :574, :576 green |
| P16 | red dashed / red dash-dot horizontals | MakeInvsens.py:577-578 / :579-580 black |
| P17 | blue / red / grey x / purple, green verticals | MakeInvsens.py:596 black / :597 red / :598 default / :599, :601 green |
| P18 | blue / black / purple, green verticals | MakeInvsens.py:616 black / :617 red / :618, :620 green |
| P19 | blue / grey x / purple, green verticals | MakeInvsens.py:636 black / :639 default / :640, :642 green |
| P19 | red / red dashed | MakeInvsens.py:643-644 red / :645-650 black dashed |
| P37, P38 | none (no figure) | — |
