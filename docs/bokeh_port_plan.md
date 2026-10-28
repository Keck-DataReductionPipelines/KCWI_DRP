# Implementation plan: Bokeh → Matplotlib port

Companion to [bokeh_replacement_proposal.md](bokeh_replacement_proposal.md), which holds the rationale, the requirements (R1–R13) and the options comparison. This document is the execution plan: tasks, the model best suited to each, and how each task is verified.

Date: 2026-10-08. Decision taken: Matplotlib behind a Bokeh-shaped wrapper; call sites keep their current form.

## How to read the model assignments

| Model | ID | Used here for |
|---|---|---|
| Fable 5.1 | `claude-fable-5-1` | The final independent review only. Optional: Opus 5.5 in a fresh session is an acceptable substitute |
| Opus 5.5 | `claude-opus-5-5` | Design work, the live viewer, and anything where a wrong choice propagates to all 48 plots |
| Sonnet 5.5 | `claude-sonnet-5-5` | Well-specified implementation: mechanical edits, tests, scripts, clean-up |
| Haiku 5.5 | `claude-haiku-5-5` | Inventory and sweep checks with a known right answer |

- The assignments are a judgement from each model's published positioning. None of them has been benchmarked on this repository.
- The assignments are chosen to minimise rework: the stronger models go where an error is expensive to find later, not where there is the most typing.
- Tests are deliberately assigned to a different model from the code they test, so the tests are not written from the implementer's own assumptions.
- In Claude Code, switch with `/model`, or have a task run as a subagent with the named model.
- Tasks marked **Human** need a person: they depend on real data, an interactive session, a machine with Firefox, or a judgement call about what "matches".

## One change from the proposal

The proposal suggested shipping PNG output first while Bokeh still drove the live view, which would have meant building each saved plot twice. This plan avoids that. The wrapper lands with a viewer interface from the start, and the live view moves in the same branch. Nothing is released with both libraries active.

## Core design decision: the wrapper records, then renders

The wrapper's `figure()` returns an object that **records** glyph calls and property assignments, and only builds a Matplotlib `Figure` when asked to render (on save, or on display). This one decision resolves several problems:

- **Order independence.** Bokeh lets a primitive set `legend.location` or `x_range` before or after adding glyphs. Rendering at the end means the wrapper never has to care.
- **Thread safety.** The live viewer renders its own `Figure` from the recorded description on its own thread. No Matplotlib object is ever shared between the `action_loop` thread and the viewer thread, which removes the main risk identified in the proposal.
- **Reuse.** `SubtractOverscan` puts the same figures into two `gridplot` calls. Recorded descriptions can be rendered any number of times.
- **Testability.** Tests can assert on the recorded description without rendering anything.

Four rules for the wrapper:

1. **Strict in tests, tolerant in production.** An unsupported method or keyword is never silently ignored, because that would produce a plot that is wrong but looks plausible. What happens instead depends on a module-level `STRICT` flag:
   - `STRICT = True` (set by the test suite, and by the environment variable `KCWI_PLOT_STRICT=1`): raise `TypeError` at the call.
   - `STRICT = False` (the default in a reduction): log one warning naming the primitive file, line and the unsupported name, skip that one call or keyword, and carry on. This is what keeps requirement R10 (a plot failure never stops a reduction) true. An exception raised inside a primitive would abort the action.
   - The guarantee that production never reaches this path comes from the coverage scan in 1.2, not from the runtime behaviour.
2. **No plotting call may raise into a primitive.** `figure()`, every glyph method, every property assignment, `gridplot`, `save_plot`, `bokeh_plot` and `bokeh_clear` catch their own exceptions when `STRICT` is false, log a warning with the traceback at debug level, and return. `figure()` and `gridplot` still return a usable recorder object on failure, so the lines that follow do not hit `AttributeError`.
3. **Snapshot data at record time.** Every array-like passed to a glyph call is copied when recorded (`numpy.array(x, copy=True)`), and scalars and lists are copied likewise. Primitives keep working on their arrays after the call, some of them in place, and the viewer renders later on another thread. Without the copy the thread-safety claim above does not hold. The recorded description is treated as immutable once handed to `bokeh_plot` or `save_plot`; the viewer receives its own deep copy of the description, so a later property assignment on `p` by the primitive cannot race with a render.
4. **No `pyplot`.** Only `matplotlib.figure.Figure` and `FigureCanvasAgg`, so there is no process-global state. For the viewer this means `matplotlib.backends.backend_webagg_core` with our own Tornado application. The higher-level `backend_webagg` module looks figures up in pyplot's global registry and must not be used.

## Task overview

| ID | Task | Model (effort) | Depends on |
|---|---|---|---|
| 0.1 | Capture reference plots from the current pipeline | Human | – |
| 0.2 | Build the plot inventory | Haiku 5.5 | – |
| 0.3 | Prototype the live viewer | Opus 5.5 (high) | – |
| 0.4 | Confirm downstream, remote-viewing and viewer-port assumptions | Human | – |
| 0.5 | Clean the working tree so diffs are readable | Human | – |
| 1.1 | Write the wrapper and PNG path | Opus 5.5 (high) | 0.2, 0.5 |
| 1.2 | Write the wrapper tests | Sonnet 5.5 (medium) | 1.1 |
| 2.1 | Convert the 11 straightforward primitives | Sonnet 5.5 (medium) | 1.1, 1.2 |
| 2.2 | Convert the 5 primitives with special features | Opus 5.5 (medium) | 1.1, 1.2 |
| 3.1 | Implement the live viewer and `StartPlotViewer` | Opus 5.5 (high) | 0.3, 1.1 |
| 3.2 | Update the three launch scripts | Sonnet 5.5 (medium) | 3.1 |
| 4.1 | Write the side-by-side comparison script | Sonnet 5.5 (low) | 0.1 |
| 4.2 | First-pass comparison of old and new plots | Opus 5.5 (high) | 0.1, 2.1, 2.2, 4.1 |
| 4.3 | Human review and sign-off; fix mismatches | Human + Opus 5.5 | 4.2 |
| 5.1 | Remove Bokeh, Selenium, and Firefox-only config | Sonnet 5.5 (medium) | 4.3 |
| 5.1b | Update the ReadTheDocs sources on the `documentation` branch | Sonnet 5.5 (medium) | 5.1 |
| 5.2 | Sweep for leftovers | Haiku 5.5 | 5.1, 5.1b |
| 5.3 | Independent review of the whole branch | Fable 5.1 (high), or Opus 5.5 | 5.2 |
| 5.4 | Clean-machine install and reduction | Human | 5.3 |

Tasks 0.1–0.5 are independent and can run in parallel. So can 2.1 and 2.2, and 3.x alongside 2.x.

## Phase 0: baseline and decisions

### 0.1 Capture reference plots — Human

On a machine where the current pipeline and Firefox work, reduce one blue and one red dataset that include bias, flat, arc, continuum bars and a standard star. Do one unattended pass at the usual `plot_level`, and one pass at `plot_level = 3`. Save each `plots/` directory as `plots_ref_<channel>_<level>/`. During the level-3 pass, screenshot the browser for the plots that are never saved, at least those in `MakeInvsens`, `ReadAtlas`, `SolveArcs` and `MakeMasterFlat`.

The level-3 pass is driven by terminal prompts, and some answers change the reduction (the atlas offset in `ReadAtlas`, the wavelength limits, line masks and fit order in `MakeInvsens`, the peak in `MakeMasterFlat`). It must be made repeatable, or the phase 4 comparison compares two different reductions:

1. Run the level-3 pass with the keyboard input captured, one file per channel:
   ```
   tee answers_<channel>_3.txt | reduce_kcwi <usual arguments>
   ```
   Every line typed, including empty returns, lands in the file in order.
2. Give at least one non-default answer at each prompt that accepts a value, so the replay exercises those paths: a new offset in `ReadAtlas`, new limits and one line mask in `MakeInvsens`, a changed fit order.
3. Immediately replay on the same unported code to prove the file is sufficient, starting from a clean output directory:
   ```
   reduce_kcwi <same arguments> < answers_<channel>_3.txt
   ```
   The replay must run to completion without waiting for the keyboard and without `EOFError`, and must write the same set of PNG names.
4. Keep the reduced output of the level-3 pass (the final cubes and the `.proc` table) alongside the plots, so 4.2 can confirm the ported run reduced the data the same way.

Some PNGs exist only at particular levels. `oscan_<frame>.png` is written only at level 3. `SolveArcs` writes `_nlines.png` and `_resid.png` at every level, including 0. The reference set must therefore also include one short pass at `plot_level = 0`, saved as `plots_ref_<channel>_0/`.

Verification:
- Every file name pattern passed to `save_plot` in the inventory (0.2) appears at least once in the reference set, or the gap is written down with the reason.
- The raw data, config, command line and `answers_*.txt` files are recorded so the same run can be repeated after the port.
- The replay in step 3 completed unattended and produced the same PNG file names as the recorded pass.
- The `plot_level = 0` directory contains the two `SolveArcs` files and nothing unexpected.

### 0.2 Plot inventory — Haiku 5.5

Produce `docs/plot_inventory.md`: one row per `figure()` call with file, line, the variable name the figure is bound to, title string, the `plot_level` condition guarding it, whether it is saved (and the file name pattern and the `plot_level` condition on the `save_plot` call, which is not always the same as the one on the figure) or live-only, and which Bokeh features it uses.

The inventory ends with a **feature summary**: one deduplicated list each of
- every keyword passed to `figure()`,
- every glyph method and every keyword passed to each,
- every attribute read or assigned on a figure (for example `x_range`, `y_range.start`, `y_range.end`, `x_range.range_padding`, `legend.location`, `xgrid.grid_line_color`, `extra_x_ranges`),
- every other name imported from `bokeh` (`Range1d`, `LinearAxis`, `gridplot`, `ColumnDataSource`, `basicConfig`) with the keywords used on it,
- every prompt (`input(` call) with its file, line and guarding condition, marking the ones whose answer changes the reduction.

The feature summary is the specification for 1.1 and the expected set for the coverage scan in 1.2. It must include what `set_plot_lims` and `oplot_slices` in `kcwi_plotting.py` do to a figure, since those run on every figure passed to them.

Why Haiku: this is extraction with a checkable answer.

Verification:
- Row count is 48, across 16 primitives.
- Totals match `grep -c` over `kcwidrp/primitives`: 49 `bokeh_plot(` calls and 20 `save_plot(` calls. (The proposal says 21; the extra one is the `def save_plot` line in `kcwi_plotting.py`.) `StartBokeh.py` contains a 49th `figure(` call that is the placeholder for the Bokeh session, not a plot, and is not a row.
- The feature summary contains at least these items, which are easy to miss: `figure(y_range=...)` (twice in `MakeMasterFlat`), `figure(y_axis_type='log')` (`MakeInvsens`), `line_alpha` on `circle` (`MakeMasterSky`, `MakeMasterFlat`), `scatter` with no `marker` (`FitCenter`), `Range1d(start=, end=)` in keyword form, `y_range.end` (via `set_plot_lims`), and `image(..., palette=, level=)`.
- The summary does **not** list `map` or `close` as figure methods. In `MakeCube` and `FitCenter` the name `p` is also bound to a multiprocessing pool.
- Spot-check five rows by hand against the source.

### 0.3 Live viewer prototype — Opus 5.5 (high)

A throwaway script outside the package that reproduces the framework's threading: the main thread idle, a daemon worker thread that builds plots and blocks in `input()` or `time.sleep()`, and a WebAgg server on a third thread showing whatever the worker last submitted.

The prototype must reproduce these facts about the real process, which were checked against `keckdrpframework` and the primitives:

- **Exit is `os._exit`.** The framework's SIGINT handler only sets a flag; the action thread then calls `on_exit`, whose default is `os._exit(status)`. `atexit` handlers and `finally` blocks on other threads do not run. The prototype must exit the same way, and the viewer design must need no clean-up at exit: daemon threads only, no child processes, no temporary files, no lock files.
- **The process forks while the viewer is live.** `MakeCube` creates a default `multiprocessing.Pool()` on the action thread, which is `fork` on Linux under Python 3.12 and 3.13. The prototype must create a default `Pool()` from the worker thread and map a function over it while the viewer thread is serving a connected browser, and repeat this while a render is in progress.
- **The viewer is built on `backend_webagg_core`**, with our own Tornado application and no `pyplot` (wrapper rule 4).
- **Figures change size and shape.** Consecutive submissions differ in pixel size (a 1000 × 600 single plot, then a 1000 × 600 or 500 × 600 grid from `SubtractOverscan`), in number of axes, and in having a twin axis. The prototype must show how one browser tab follows these changes: either one long-lived `Figure` that is cleared and re-rendered with the toolbar's zoom history reset, or a new figure per submission with the page told to reconnect. The written decision states which.

It must also settle the viewer topology. Today every pipeline process pushes into one shared `bokeh serve` on port 5006 and a browser tab opens automatically. The default design to prototype is:

- A new config key `plot_viewer_port` (default 5006, so existing SSH port-forwards keep working). The viewer binds to `localhost` only.
- If the port is taken, try the next ones up to `plot_viewer_port + 9`. If none is free, log a warning and run with no viewer.
- The chosen URL is logged at INFO and printed to the terminal on its own line, so the user of each terminal knows which tab belongs to it.
- A new config key `plot_viewer_open_browser` (default `True`, matching today's `session.show`). The tab is opened with the standard-library `webbrowser` module. Failure to open one, including on a headless host, is a debug-level log line and nothing more.
- Each pipeline process shows only its own plots. This is a change from today, where all processes share one tab, and it is recorded as an accepted difference in the changelog.

If 0.4 shows observers rely on a single fixed port for several processes, the alternative to prototype is: only the first process to bind the port runs a viewer and the others run with none.

Why Opus at high effort: this is the only part of the port with real concurrency, an event loop on a non-main thread, and behaviour that differs between macOS and Linux, and a wrong conclusion here would misdirect phase 3. The concurrency is modest, though: one Tornado loop on a daemon thread with a single thread-safe hand-off, a pattern Matplotlib documents in its WebAgg embedding example, and the hard parts are spelled out above. The protection against a wrong design is the pass or fail criteria below, run on both operating systems, not the model tier. An earlier draft assigned this to Fable 5.1; that was reviewed and judged unnecessary.

Verification (all must hold for WebAgg to be chosen):
- A plot submitted from the worker appears in the browser while the worker is blocked in `input()`.
- 200 consecutive submissions at 0.1 s intervals produce no exception, no hang, and no growth in open figures or memory.
- Zoom, pan and the x/y cursor readout work on a displayed plot.
- Two instances on one host both come up, on different ports, without interfering. Each prints its own URL, and each tab shows only that instance's plots.
- With all ten ports taken, the instance logs a warning and its submissions are harmless no-ops.
- With no browser connected and `DISPLAY` unset, submissions are harmless no-ops and no browser launch is attempted more than once.
- Ctrl-C exits the process through `os._exit` with no orphaned process, no hang at shutdown, and nothing left on disk. The port can be bound again immediately by a new instance.
- A default (`fork`) `Pool()` created on the worker thread while the viewer is serving completes 50 times in a row with no hang, including when started during a render. The pool's children do not keep the viewer port open after the parent exits. Record the result for `fork` and for `get_context("spawn")` separately.
- Submitting figures of different pixel sizes and layouts in sequence leaves the tab showing each one at the right size, with zoom and pan working on the newest and no stale zoom state from the previous one.
- The worker mutates an array in place immediately after submitting it; the displayed plot shows the values as they were at submission.
- A browser tab that is closed and reopened, or a second tab opened on the same URL, shows the current plot.
- All of the above pass on macOS and on Ubuntu 24.04.

Output: a short written decision added to this document, covering: WebAgg or fallback; one figure reused or one per submission; the port and browser-opening behaviour; and the fork result. If the `fork` criterion fails but `spawn` passes, the decision is to switch `MakeCube` to `get_context("spawn").Pool()`, as `FitCenter` already does, and that change is added to task 2.2. If any other criterion fails and cannot be fixed, phase 3 uses the fallback (an auto-refreshing page showing the latest PNG) and the loss of zoom and cursor readout is raised with the team before continuing. The port, URL-printing, exit and fork criteria apply to the fallback as well.

### 0.4 Confirm assumptions — Human

- Whether anything downstream (RTI, KOA ingestion) reads `plots/*.png`, and whether it depends on more than file name and format.
- How observers view plots today when the pipeline runs remotely (VNC, or an SSH port-forward to port 5006).
- Whether several pipeline processes are ever watched at once, and if so whether one shared tab (today's behaviour) is relied on or one tab per process is acceptable. This decides between the two viewer topologies in 0.3.
- Whether the automatic opening of a browser tab at start-up is wanted on the summit and RTI hosts, or should default to off there (`kcwi_koarti.cfg`).
- Who can merge to the `documentation` branch and how ReadTheDocs is rebuilt, for task 5.1b.

Verification: all answers recorded in this document. The two viewer questions are needed before 3.1 starts; the rest before phase 4.

### 0.5 Clean the working tree — Human

At the time of writing, `git status` on `master` shows every tracked file as modified. The change is a file-mode flip from 644 to 755 on about 300 files, with almost no content change, and `docs/` is untracked. Several verifications below rely on reading `git diff`, so this has to be settled first.

1. Decide whether the mode change is wanted. If not, restore it with `git diff --summary | grep 'mode change'` to confirm the scope, then either `git config core.filemode false` for this clone or reset the modes. If it is wanted, commit it on its own, separately from the port.
2. Review and commit, or stash, the handful of real content changes that are mixed in (`git diff --stat` reports about 27 inserted and 10 deleted lines).
3. Commit `docs/bokeh_replacement_proposal.md` and this plan.
4. Create the port branch from that clean commit and tag the starting point.

Verification:
- `git status --short` on the new branch is empty.
- `git diff --summary <tag>` shows no mode changes at any later point in the port; this is re-checked in 5.3.

## Phase 1: wrapper and PNG output

### 1.1 Wrapper — Opus 5.5 (high)

Rewrite `kcwidrp/core/kcwi_plotting.py` and the internals of `kcwidrp/core/bokeh_plotting.py`. Keep the public names `bokeh_plot`, `bokeh_clear`, `save_plot`, `get_plot_lims`, `set_plot_lims`, `oplot_slices` and `check_running_process`, so primitives do not change how they call them.

The wrapper must support exactly what the feature summary in the inventory (0.2) lists. The table below was checked against the source and is the minimum; if the inventory finds more, the inventory wins and this table is updated.

| Bokeh usage | Notes |
|---|---|
| `figure(title, x_axis_label, y_axis_label, plot_width, plot_height, tooltips)` | `tooltips` accepted and ignored; size maps to pixels at 100 dpi |
| `figure(y_range=[lo, hi])` | Used twice in `MakeMasterFlat`, with a plain two-element list. Accept a list, tuple or `Range1d`; sets fixed y limits |
| `figure(y_axis_type='log')` | Used once in `MakeInvsens` (inverse sensitivity). Maps to a log y scale. Any other value is unsupported. Non-positive values must not raise at render |
| `line`, `circle`, `diamond`, `scatter`, `quad`, `image` | The six glyph methods. Anything else is unsupported |
| Glyph keywords | `color`, `line_color`, `fill_color`, `line_width`, `line_dash`, `alpha`, `line_alpha`, `size`, `legend_label`, `x_range_name` |
| `line_alpha=0.` with `fill_color=` on `circle` | The idiom used in `MakeMasterSky` and `MakeMasterFlat` for dense point clouds: filled marker with no outline. Must map to a marker with no edge, not to an invisible marker |
| `color=` on a marker glyph | Sets both outline and fill, as in Bokeh. `line_color` and `fill_color` override each part |
| `scatter(marker='x')` and `scatter(...)` with no marker | The default marker is a circle (`FitCenter`). `'x'` is the only named marker in use |
| `quad(top=, bottom=, left=, right=)` | Keyword-only, array or scalar per edge (`MakeMasterBias` histogram). Renders as rectangles with `fill_color`, `line_color`, `alpha` |
| `image([array], x=, y=, dw=, dh=, palette=, level=)` | One 2-D array inside a list. `palette="Spectral11"` only; `level` accepted and ignored |
| `p.line('x', 'y', source=ColumnDataSource(...))` | Not supported by the wrapper. The one use, in `MakeInvsens`, is rewritten in 2.2 to pass arrays |
| `Range1d(a, b)` and `Range1d(start=a, end=b)` | Provide a local `Range1d` with both forms |
| `p.x_range = Range1d(...)`, `p.y_range = Range1d(...)` | Replace the range object |
| `p.x_range.start`, `.end`, `p.y_range.start`, `.end` | Each settable on its own, leaving the other end automatic. `set_plot_lims` assigns all four; `SolveArcs` assigns only `y_range.start = 0` |
| `p.x_range.range_padding = p.y_range.range_padding = 0` | Chained assignment in `MakeCube`. Zero margins on that axis |
| `p.legend.location = ...` | Bokeh names (`top_left`, `top_center`, `bottom_center`, …) mapped to Matplotlib's. Assigned before or after the glyphs |
| `p.xgrid.grid_line_color = None` | Turns off the vertical grid lines. Any value other than `None` is unsupported |
| `p.extra_x_ranges = {name: Range1d}` + `p.add_layout(LinearAxis(x_range_name=, axis_label=), 'above')` | Provide a local `LinearAxis`; render with `twiny()`. Only `'above'` is supported |
| `gridplot(plts, ncols=, plot_width=, plot_height=, toolbar_location=)` | Renders the recorded figures into subplots. `plot_width` and `plot_height` are per cell and override each figure's own size. `toolbar_location` accepted and ignored. The same figure list may be passed to several `gridplot` calls |

Appearance defaults, set once: default glyph colour `#1f77b4` with no colour cycling; legend top-right with entries merged by label; marker size equivalent to Bokeh's 4 px; left-aligned bold title; light grey grid; automatic axis limits padded as Bokeh's `DataRange1d` pads them (10 % of the span by default), not with Matplotlib's 5 % margins.

Known mapping traps. Each has a named test in 1.2:

- **Legend across the twin axis.** In `SolveArcs` three glyphs are on the main x-axis and five on the `"pix"` axis, and all share one legend of six entries, with `'Arc'` and `'Atlas'` each used by two glyphs. A Matplotlib legend only sees the artists of its own axes. The wrapper must collect handles from both axes, merge by label, and draw one legend on the axes that is on top.
- **Draw order with a twin axis.** `twiny()` axes are drawn after the base axes, so everything on the `"pix"` range covers everything on the main range regardless of call order. Check the `SolveArcs` plot against the reference and set z-order explicitly if the result differs.
- **`Spectral11` is eleven discrete colours.** Build a `ListedColormap` from Bokeh's eleven hex values, in Bokeh's order (`#5e4fa2` first, `#9e0142` last), mapped linearly from the array minimum to its maximum. Matplotlib's built-in `Spectral` is continuous and runs in the opposite direction, so using it would both smooth and invert the colours.
- **Image orientation.** Bokeh's `image` puts row 0 at the bottom; Matplotlib's `imshow` puts it at the top by default. Use `origin='lower'` and an explicit `extent` of `(x, x + dw, y, y + dh)`.
- **NaN handling.** Bokeh leaves a gap in a line at NaN and drops NaN markers; confirm the wrapper does the same and that an all-NaN series renders an empty plot and does not raise.
- **Marker size units.** Bokeh `size` is a diameter in screen pixels. Matplotlib's `markersize` is a diameter in points, and `scatter`'s `s` is an area in points squared. At 100 dpi one pixel is 0.72 points. `size=1` point clouds must stay one pixel wide.
- **Line width units.** Bokeh `line_width` is in pixels, Matplotlib's in points; same 0.72 factor.
- **Dash names.** `'dashed'`, `'dashdot'`, `'dotted'` map by name; any other value is unsupported.
- **Fixed pixel size.** Save with the figure's own size and dpi and never with a tight bounding box, which would change the pixel size. Long titles and the twin-axis label must fit inside the fixed size; use constrained layout.

`save_plot` keeps today's behaviour apart from the renderer: it writes `plots/<filename>` relative to the current directory, creates `plots/` if missing, accepts a figure or a `gridplot` result, logs the same `>>> Saving to` line, and does not depend on the viewer or on `enable_bokeh`. It is called at `plot_level` 0 by `SolveArcs`, so it must work with no viewer and no `StartPlotViewer` event having run.

Why Opus: every later task and all 48 plots inherit these choices, and the Bokeh-to-Matplotlib mapping has several unit and default traps.

Verification:
- `pytest kcwidrp/tests` passes, including the tests from 1.2.
- Importing `kcwidrp.core.kcwi_plotting` in a fresh interpreter leaves `matplotlib.pyplot`, `bokeh` and `selenium` out of `sys.modules`.
- A saved PNG of a single figure is exactly `plot_width` × `plot_height` pixels. A saved grid is exactly `ncols × plot_width` by `nrows × plot_height` using the sizes given to `gridplot`; confirm this matches the `oscan_*.png` reference from 0.1 and record the reference size here.
- With `STRICT` true, an unsupported method, keyword or attribute raises `TypeError`. With `STRICT` false, the same call logs exactly one warning, raises nothing, and the figure still saves.
- With `STRICT` false, a render that fails inside Matplotlib (forced in a test by patching the canvas to raise) logs a warning and `save_plot` returns normally.
- Mutating an input array in place after the glyph call does not change the saved PNG.
- Every row of the table above and every trap in the list has a passing test.
- Rendering runs with `DISPLAY` unset.

### 1.2 Wrapper tests — Sonnet 5.5 (medium)

Replace `kcwidrp/tests/test_bokeh_plotting.py` with `test_kcwi_plotting.py`. Written from the inventory and the table above, not from the wrapper's source.

- **Structural tests.** For each supported feature, build a figure with synthetic data, render it, and assert on the Matplotlib axes: number of lines, colours, dash styles, limits, labels, title, legend entries and position, the second x-axis, subplot layout, image extent.
- **Coverage scan.** A test that parses every file in `kcwidrp/primitives` and `kcwidrp/core/kcwi_plotting.py` with `ast` and collects each method, keyword and attribute used on a figure object, then asserts the wrapper supports all of them. This is what guarantees that no call site reaches an unsupported feature at 3 a.m. during a reduction. It cannot work from variable names, so it is specified as follows:
  - A name is a figure only from the point where it is assigned the result of a call to `figure(...)`, and stops being one when it is assigned anything else. `p` is rebound many times per function, and in `MakeCube` and `FitCenter` it is first a multiprocessing pool (`p.map`, `p.close`) and later a figure. Tracking is per function, in source order.
  - Figure names in use are `p`, and in `MakeInvsens` also `peff`, `pivs`, `prsd` and `pcal`. The scan must not hard-code these; it discovers them.
  - A figure passed to `set_plot_lims`, `oplot_slices`, `bokeh_plot`, `save_plot` or appended to a list later given to `gridplot` is followed into those helpers: the helpers' own use of the figure is scanned as part of `kcwi_plotting.py`.
  - Collected per figure: keywords of the `figure(...)` call; each `<fig>.<method>(...)` with its keywords; each attribute chain read or assigned (`x_range.start`, `legend.location`, `xgrid.grid_line_color`, …), including both targets of a chained assignment.
  - Also collected: every call to `Range1d`, `LinearAxis` and `gridplot` with its keywords, and the literal values of `marker=`, `line_dash=`, `palette=`, `y_axis_type=`, `legend.location` and the second argument of `add_layout`, which are checked against the wrapper's supported values.
  - A string argument in the x or y position of a glyph call (the `ColumnDataSource` column-name form) is reported as unsupported.
  - Any use the scan cannot classify (a figure passed to an unknown function, an attribute set through `setattr`, keywords passed as `**kwargs`) fails the test with the file and line. Unknown must never pass silently.
  - The scan's result is also compared with the feature summary in `docs/plot_inventory.md`; a difference in either direction fails.
- **Strict and tolerant modes.** The whole suite runs with `STRICT` true. A separate group runs with `STRICT` false and checks, for an unsupported method, an unsupported keyword, an unsupported attribute and an unsupported literal value: one warning is logged with the caller's file and line, nothing is raised, the returned object still accepts further calls, and `save_plot` still writes a file.
- **Snapshot test.** Record a line from an array, overwrite the array in place, render: the rendered data is the original. The same for a figure submitted to a stub viewer and then modified by a further glyph call on the primitive's side.
- **Trap tests.** One named test per item in the 1.1 trap list: twin-axis legend merged across both axes; twin-axis draw order; `Spectral11` endpoint colours and discreteness; image orientation with an asymmetric array; NaN gaps and all-NaN series; one-pixel markers; pixel line width; dash names; fixed pixel size with a very long title.
- **Thread test.** Four threads each render and save ten plots concurrently; all files exist with the right pixel size.
- **Failure isolation.** `save_plot` with an unwritable directory logs a warning and returns; it does not raise. The same for `bokeh_plot` with a viewer whose `show` raises, and for `bokeh_clear` with a viewer of `None` (today's `bokeh_clear` would raise on `None`; the new one must not).
- **Plot level 0.** With no viewer started, `save_plot` on a figure writes its PNG. This protects the two files `SolveArcs` writes at every plot level.

Why Sonnet: well-specified, and independence from the implementer matters more than raw capability.

Verification:
- Tests fail when a deliberately broken wrapper is substituted (for example, colour cycling re-enabled, or `legend.location` ignored). Check at least three such mutations.
- The coverage scan fails if an unsupported call such as `p.wedge(...)` is added to a primitive, and also for each of: an unsupported keyword on a supported glyph; `figure(x_axis_type='log')`; `p.ygrid.grid_line_color = None`; a figure bound to a new name (`q = figure(...)`; `q.wedge(...)`); a figure passed to an unknown helper.
- The coverage scan does **not** flag `p.map(...)` and `p.close()` on the pools in `MakeCube` and `FitCenter`.
- Before phase 2 the scan is run against the unconverted primitives and its only failures are the items that 2.2 is due to remove (`ColumnDataSource` and the column-name `line` call in `MakeInvsens`).
- The suite runs in under 30 seconds with no display and no network.

## Phase 2: call sites

The only intended change in a primitive is its import block: `figure`, `Range1d`, `LinearAxis` and `gridplot` come from `kcwidrp.core.kcwi_plotting` instead of `bokeh`.

### 2.1 Straightforward primitives — Sonnet 5.5 (medium)

`ArcOffsets`, `ExtractArcs`, `FitCenter`, `GetAtlasLines`, `MakeMasterBias`, `MakeMasterFlat`, `MakeMasterSky`, `ReadAtlas`, `SolveAIT`, `SubtractScatteredLight`, `TraceBars`. These use only `figure` and, in some, `Range1d`. `MakeMasterFlat` is large (13 figures) and is the only user of `figure(y_range=...)`; `MakeMasterSky` and `MakeMasterFlat` use the `line_alpha=0.` point-cloud idiom; `FitCenter` uses `scatter` with no marker. All three are wrapper features from 1.1, so these files still need import changes only. `FitCenter` also binds `p` to a multiprocessing pool before it binds it to a figure; that code is not touched.

Verification:
- `git diff` for these 11 files touches import lines only. Any other changed line must be justified individually.
- `pytest kcwidrp/tests/test_import_primitives.py` passes.
- The coverage scan from 1.2 passes.

### 2.2 Primitives with special features — Opus 5.5 (medium)

| Primitive | What needs care |
|---|---|
| `SolveArcs` | The only twin x-axis; glyphs are bound to the second range by `x_range_name` |
| `MakeInvsens` | `ColumnDataSource` with column-name glyph calls and `tooltips`. Pass the arrays directly and drop the data source. Also the only log y-axis (`y_axis_type='log'`), on a plot that draws vertical marker lines from the data minimum to maximum. The printed instruction "Hover the cursor over spectrum to find wavelengths" stays, and must be true of the new viewer's cursor readout |
| `SubtractOverscan` | `gridplot`, with the same figures reused in two grids, shown per amplifier first and then as a grid. Also the only caller of `bokeh_clear` |
| `MakeCube` | The only `image` call, with the `Spectral11` palette and zero range padding. Creates a default (`fork`) `Pool()` shortly before plotting; if the 0.3 decision says so, change it to `get_context("spawn").Pool()` as in `FitCenter`, and confirm `make_cube_helper` and its arguments still pickle |
| `FindBars` | Calls `bokeh.util.logconfig.basicConfig`; remove it |

Why Opus: each of these is a one-off where the right translation needs a decision, and each is in code interleaved with reduction logic.

Verification:
- Same three checks as 2.1; here non-import changes are expected and each is listed in the commit message.
- A dedicated structural test for each of the five, using synthetic data shaped like the real inputs.
- `MakeCube`: image orientation checked against the reference with an asymmetric test array, since `origin` conventions differ between the libraries. Colours at the array minimum and maximum match Bokeh's `Spectral11` first and last entries.
- `MakeCube`, if the pool was switched to `spawn`: one blue and one red cube from the 0.1 data are numerically identical to the reference cubes, and the wall-clock time of the primitive is recorded before and after.
- `SolveArcs`: the legend shows all six labels (`Atlas`, `AtCntr`, `Arc`, `Cntr`, `Gpeak`, `Peak`) once each, and the plot matches the 0.1 screenshot for which series is drawn on top.
- `MakeInvsens`: the inverse-sensitivity plot renders with a log axis when the data include zero or negative values, with a warning at most and no exception.

## Phase 3: live view

### 3.1 Viewer and `StartPlotViewer` — Opus 5.5 (high)

Implement the viewer chosen in 0.3 behind a small interface: `start()`, `show(recorded_figure)`, `clear()`, `stop()`. `bokeh_plot(p, viewer)` and `bokeh_clear(viewer)` delegate to it and do nothing when the viewer is `None`. Rename the `StartBokeh` primitive to `StartPlotViewer`, keeping `start_bokeh` as an alias event name in both pipeline tables, and keeping `context.bokeh_session` as the attribute name so primitives are untouched.

Details the implementation must follow:

- **Hand-off.** `show()` deep-copies the recorded description on the caller's thread and passes the copy to the Tornado loop with `IOLoop.add_callback`, which is the only thread-safe entry point to the loop. All Matplotlib objects of the viewer are created, drawn and destroyed on the viewer thread only. `show()` returns without waiting for the render.
- **`show()` never blocks and never raises.** If the viewer thread has died, `show()` logs one warning, marks the viewer dead, and later calls return at once.
- **Ports and URL.** As decided in 0.3: `plot_viewer_port` with the search range, `localhost` binding, the URL logged and printed, `plot_viewer_open_browser`. `StartPlotViewer` stores the URL on the context so later log lines can repeat it.
- **Start-up failure.** Today `terminate_on_failed_bokeh_start = True` stops the pipeline if the server is unreachable. That key is removed in 5.1; the new behaviour is always to warn and continue (R10). State this in the changelog.
- **No exit-time clean-up.** The framework leaves through `os._exit`, so `stop()` is not called at the end of a normal run and `atexit` does not fire. `stop()` exists for tests and for an orderly restart. The viewer therefore uses daemon threads only and owns no child processes, temporary files or lock files. The listening socket is created with close-on-exec set and is not inheritable.
- **Fork safety.** Whatever 0.3 found. At minimum the viewer thread must not hold the logging lock or any module-level lock across a render, and pool children must not service the viewer socket.
- **Ctrl-C while blocked in `input()`.** The framework's handler only sets `keep_going`; the action thread stays in `input()` until return is pressed. This is today's behaviour and is not changed, but the viewer must keep serving during that wait and must not swallow or re-install the SIGINT handler. Tornado must be started without its own signal handling.
- **Launch scripts.** `enable_bokeh` is read in the three scripts to decide whether to queue the `start_bokeh` event; see 5.1 for the renamed key. The 5-second `time.sleep` that waits for `bokeh serve` goes away with the subprocess in 3.2.

Why Opus at high effort: same reasoning as 0.3, now in production form, including shutdown behaviour and the multi-process case. The implementation details above and the verification list below carry the risk.

Verification:
- The 0.3 criteria, re-run against the real implementation inside `reduce_kcwi`.
- With the viewer port already taken, the pipeline uses the next free one and prints it. With the whole range taken, it logs a warning and reduces normally, writing all PNGs.
- With `enable_bokeh = False` (and with the new key false), no viewer thread or socket is created, and all PNGs for the configured `plot_level` are still written.
- At `plot_level = 0` with the viewer enabled, `SolveArcs` still writes `_nlines.png` and `_resid.png`.
- A full unattended reduction leaves no extra processes, temporary files or bound ports behind (`ps` and `lsof -i` before and after).
- Two `reduce_kcwi` processes on a shared queue (`want_multiprocessing = True`) both complete, with all PNGs written, each showing its own plots at its own printed URL.
- `MakeCube` runs to completion 20 times with a browser connected and the viewer mid-update, on Ubuntu, with the pool start method chosen in 0.3.
- Killing the browser tab, and separately stopping the viewer thread from a debugger hook, mid-run does not stop the reduction; one warning is logged.
- Ctrl-C during `plot_pause` and Ctrl-C followed by return during an `input()` prompt both end the process with nothing left running.
- An interactive `MakeInvsens` session: the cursor readout gives wavelengths precise enough to type the limits the prompt asks for (compare with the tooltip format `0.0` used today).

### 3.2 Launch scripts — Sonnet 5.5 (medium)

In `reduce_kcwi.py`, `kcwi_rti.py` and `smart_reduce_kcwi.py`: remove the `bokeh serve` subprocess launch, the `time.sleep(5)` that waits for it, and the `bokeh_output.txt` log. Fix the broken `check_bokeh_server` import in `smart_reduce_kcwi.py`, and remove its `open http://localhost:5006?...` subprocess call; opening the tab is now the viewer's job. Delete `scripts/boklean`. Leave the `subprocess`, `psutil` and `check_running_process` imports only where something else in the file still uses them.

All three scripts decide whether to queue `start_bokeh` from the config. Route that through one helper that reads the new key and falls back to `enable_bokeh` (see 5.1), so the three scripts cannot drift apart.

Verification:
- All three entry points start with `--help` and exit 0.
- `smart_reduce_kcwi` imports without error (it does not today).
- No `subprocess.Popen` call referring to `bokeh` or to port 5006 remains.
- Start-up of `reduce_kcwi` with the viewer enabled is at least the 5 seconds faster that the removed sleep accounts for.

## Phase 4: confirming the plots match

### 4.1 Comparison script — Sonnet 5.5 (low)

A small script that takes a reference directory and a new directory and writes one HTML page with old and new PNGs side by side, paired by file name, listing files present on only one side.

Verification: run it on a reference directory against itself (every pair identical, nothing unpaired), and against a copy with one file removed and one renamed (both reported).

### 4.2 First-pass comparison — Opus 5.5 (high)

Re-run the 0.1 reductions on the ported branch with the same data, config and command lines: the unattended pass, the `plot_level = 0` pass, and the level-3 pass replayed from the recorded answers with `reduce_kcwi <same arguments> < answers_<channel>_3.txt`. No one types during the level-3 replay.

The replay doubles as a check that the port did not change the prompts. If it stops waiting for input, ends with `EOFError`, or has unread lines left in the answers file, the number or order of prompts has changed. That is a behaviour change in a primitive and is a finding in its own right, to be resolved before comparing any images.

Before looking at plots, confirm the two level-3 runs reduced the data the same way: the final cubes and the `.proc` table match the 0.1 outputs (numerically identical unless 0.3 changed the `MakeCube` pool, in which case within floating-point tolerance with the reason stated).

For each image pair, check and record: same file name; same pixel size; same title and axis labels; same axis limits; same number of series; same colours and dash styles; same legend entries and position; no clipped or overlapping text.

Why Opus: this needs careful reading of plot images against a checklist, and misses here are what reach users.

Verification:
- A table with one row per pair and a pass or a specific stated difference. "Looks similar" is not an acceptable entry.
- Zero unpaired files, or each explained. This holds separately for the level-0, usual-level and level-3 directories.
- The level-3 replay consumed the answers file exactly, and the reduced outputs match as described above.
- This pass is a filter, not the sign-off. It does not replace 4.3.

### 4.3 Review and fixes — Human, with Opus 5.5 for fixes

A maintainer who knows the plots reviews the comparison page and the 4.2 table, then steps through one interactive `plot_level = 3` reduction to check the live-only plots against the 0.1 screenshots, including the interactive prompts in `MakeInvsens` and `ReadAtlas`.

Fixes go into the wrapper wherever possible, so one change corrects every affected plot. A fix made in a primitive is a sign the wrapper mapping is wrong and should be questioned first.

Verification:
- Written sign-off listing accepted differences (expected: fonts, anti-aliasing, no toolbar strip in saved PNGs).
- Every fix is accompanied by a structural test that would have caught it.
- 4.2 is re-run after fixes and is clean.

## Phase 5: removal and release

### 5.1 Remove the old stack — Sonnet 5.5 (medium)

- `pyproject.toml`: remove `bokeh`, `selenium`, `phantomjs`. Add `tornado` explicitly if the WebAgg viewer was chosen; it is installed today only because Bokeh pulls it in. Test whether the `jinja2 < 3.1` pin can go, and whether anything else imports `jinja2` at all.
- `requirements.txt`: regenerate from `pyproject.toml` with the same tool that produced it (its header says it is autogenerated); do not edit by hand.
- `environment.yml`: remove `bokeh`, `selenium`, `geckodriver`, `firefox`, `phantomjs`. This file is already out of step with `pyproject.toml` (it pins `bokeh~=2.0.0` and `matplotlib~=3.1.3`, against 2.4.3 and 3.9 in `pyproject.toml`). Either bring the remaining pins into line or raise with the maintainers whether the file should be deleted; do not leave it half-updated.
- `conda_build_files/meta.yaml`: it lists no run dependencies beyond `python`, so there is nothing to remove. Confirm that and leave it alone.
- `kcwi.cfg`, `kcwi_koarti.cfg`: remove `plot_firefox_compat`, `plot_prewarm_firefox`, `terminate_on_failed_bokeh_start`, and the comment lines that describe them. Add `plot_viewer_port` and `plot_viewer_open_browser` (see 0.3), with `plot_viewer_open_browser` set as decided in 0.4 for the RTI config.
- **The enable switch.** Add `enable_plot_viewer`, not `enable_plots`. Today `enable_bokeh` controls only the live view: PNGs are saved whenever `plot_level` says so, whether or not Bokeh is enabled, and two `SolveArcs` PNGs are saved at every level. The new key keeps exactly that meaning, and its name says so. A key called `enable_plots` would invite someone to set it false and expect, or cause, PNG output to stop, which would change what RTI and KOA receive. Rules:
  - `enable_plot_viewer` present: it decides. `enable_bokeh` present as well: ignored, with a deprecation message logged once.
  - Only `enable_bokeh` present: honoured, with the same deprecation message.
  - Neither present: viewer off.
  - Neither key has any effect on `save_plot`.
  - Update the comment in `kcwi.cfg` to say that the key controls the live viewer only and that saved PNGs are governed by `plot_level`.
- `README.rst`: it has no Firefox or geckodriver steps, so no removal is needed. Add one short paragraph on the plot viewer and its URL if the README is the right place; the full documentation is handled in 5.1b.
- `changelog.txt`: note the change; the removed config keys and the fact that a failed viewer start no longer has an option to stop the pipeline; the new keys; the deprecation of `enable_bokeh`; one viewer tab per pipeline process; and the accepted visual differences.
- `.travis.yml`: it has no browser steps of its own, but it builds its environment from `environment.yml`, so it inherits whatever is decided for that file. It also targets Python 3.7, below the package minimum of 3.12; note this to the maintainers and do not fix it as part of the port.

Verification:
- `pip install .` in a new virtual environment succeeds and `pip list` shows none of the removed packages.
- A user config file that still contains the three removed keys and `enable_bokeh` loads without error.
- A table of the four combinations of `enable_plot_viewer` and `enable_bokeh` (each present or absent, true or false where present) is tested: the viewer starts or not as the rules above say, the deprecation message appears exactly when `enable_bokeh` is present, and in every row the same PNG files are written.
- `kcwidrp/tests/test_parse_configs.py` passes with the shipped configs and is extended to cover the new keys.
- `pytest kcwidrp/tests` passes in that environment.

### 5.1b Documentation branch — Sonnet 5.5 (medium)

The user documentation is not on `master`. It is built by ReadTheDocs from the `documentation` branch (`origin/documentation`, sources under `docs/source/`). Five pages there mention Bokeh, Firefox or geckodriver and must be updated on that branch, by whoever 0.4 identified as able to merge to it:

| File | Change |
|---|---|
| `docs/source/installing.rst` | Remove the `conda install -c conda-forge selenium geckodriver firefox phantomjs` step and the text around it |
| `docs/source/configuration.rst` | Remove the three deleted keys; document `enable_plot_viewer`, `plot_viewer_port`, `plot_viewer_open_browser`; note the `enable_bokeh` deprecation; state that PNG saving follows `plot_level` only |
| `docs/source/quickstart.rst` | Replace the Bokeh server and browser description with the viewer: where the URL is printed, one tab per process, how to forward the port over SSH |
| `docs/source/issues.rst` | Remove the "Firefox/geckodriver cannot be found" section and its WebDriver workaround; add an entry for the viewer port being unreachable or taken |
| `docs/source/api.rst` | `StartBokeh` becomes `StartPlotViewer`; `bokeh_plotting` entries updated to whatever module names 1.1 kept |

The `documentation` branch also carries its own copies of package files. Do not port code there; update only `docs/source/` and whatever the ReadTheDocs build needs to import the renamed primitive.

Verification:
- The docs build locally with `sphinx-build` with no new warnings, in particular none from `api.rst` about a missing `StartBokeh`.
- A search of `docs/source/` on the branch for the 5.2 terms returns only intentional hits.
- The ReadTheDocs build for the branch succeeds after the merge, and the published installation page no longer mentions a browser.

### 5.2 Sweep — Haiku 5.5

Search the repository for `bokeh`, `selenium`, `firefox`, `geckodriver`, `marionette`, `phantomjs`, `export_png`, `5006`, `boklean`, `bokeh_output`, `ColumnDataSource`, `terminate_on_failed`. Run the same search over `docs/source/` on the `documentation` branch. Classify every hit as intentional or leftover.

Expected intentional hits, and nothing else: the public function names `bokeh_plot` and `bokeh_clear` and the module `bokeh_plotting`; the `context.bokeh_session` attribute; the `start_bokeh` event alias in the two pipeline tables; the `enable_bokeh` fallback and its deprecation message; `5006` as the default of `plot_viewer_port`; the changelog; and the two planning documents in `docs/`.

Verification: the list of intentional hits is short, explicit, and matches what phases 3 and 5.1 said they would keep. Anything else is fixed.

### 5.3 Independent review — Fable 5.1 (high), or Opus 5.5 (high) in a fresh session

A review of the full branch diff in a fresh session with no access to the implementation conversations, focused on: thread and process safety of the viewer, including the fork in `MakeCube` and reliance on any exit-time clean-up that `os._exit` would skip; any path where a plotting failure could stop a reduction, including every place the wrapper can raise when `STRICT` is false; recorded data that is not copied; behaviour changes in primitives beyond imports, including any change to the number or order of prompts; any way the viewer switch could affect PNG saving; wrapper mappings that are plausible but wrong (units, defaults, axis orientation, the items in the 1.1 trap list); and stray file-mode changes in the diff (0.5).

Why a fresh context: the reviewer should not share the implementers' assumptions, and this is the last automated check before release. Most of that independence comes from the fresh session with no access to the implementation conversations. Using Fable adds a little more, because Opus wrote both the wrapper and the viewer, and it is a single pass. It is a nicety, not a requirement: if Fable is not available, run this with Opus 5.5 in a fresh session and accept slightly less independence from the author.

Verification: every finding is either fixed with a test or closed with a written reason.

### 5.4 Clean-machine acceptance — Human

On an Ubuntu 24.04 machine with **no Firefox, no Chrome and no geckodriver installed**, install from the branch and run the blue and red reference reductions unattended, then one interactive session over the remote-viewing route confirmed in 0.4.

Verification, mapped to the proposal's requirements:

| Check | Requirement |
|---|---|
| Install and full reduction succeed with no browser on the machine | R1 |
| All expected PNGs are written from the `action_loop` thread | R2 |
| Two concurrent pipeline processes complete with no stray processes, each with its own printed viewer URL | R3 |
| PNGs are written with the viewer disabled, and at `plot_level = 0` | R6 |
| Unattended run succeeds with `DISPLAY` unset | R4 |
| PNG names and sizes match the reference set | R6 |
| Live view updates during `input()` and `plot_pause`, viewed remotely | R7 |
| Zoom, pan and cursor readout work in the interactive steps | R8 |
| Killing the viewer mid-run does not stop the reduction | R10 |
| A config with an unsupported plot call injected into one primitive (test build only) logs a warning and completes | R10 |
| Time per saved plot is recorded and under 1 s | R11 |

R5, R9, R12 and R13 are covered by the phase 1, 2 and 5.1 checks.

## Rollback

Until 5.1, the old code is one `git revert` of the branch merge away, and no dependency has been removed. Tag the commit before the merge. After release, anyone who needs the Bokeh path can pin the previous version; the saved PNG names are unchanged, so mixed-version outputs remain comparable.

## Open items carried from the proposal

- The WebAgg-or-fallback decision (0.3) is the only point where the plan can branch. Everything in phases 1, 2 and 4 is unaffected by which way it goes.
- Downstream use of `plots/*.png` and the remote-viewing route (0.4) are still assumptions until confirmed.
- Viewer topology is a default, not yet a decision: one viewer per process on `plot_viewer_port` and the nine ports above it, with the URL printed. 0.4 confirms or replaces it before 3.1.
- Whether `MakeCube` moves from a `fork` pool to a `spawn` pool depends on the 0.3 fork test. It is the one reduction-code change this port may need.
- `environment.yml` is stale independently of this port. 5.1 must either fix it or get a decision to delete it.
- The proposal's count of 21 `save_plot` calls is one too high; the correct figure is 20. The proposal has not been edited.
- The time estimates in the proposal (roughly 7–11 working days) were for a person doing the work. With the model assignments above, elapsed time will be set mainly by the human tasks: capturing references, the interactive review, and the clean-machine run.
