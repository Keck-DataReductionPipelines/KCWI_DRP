# Proposal: replace Bokeh with Matplotlib to drop the Firefox dependency

Status: proposal, not yet implemented. Date: 2026-10-08.

## Summary

**Recommendation: port plotting to Matplotlib**, using the Agg backend for saved PNGs and the WebAgg backend for the live browser view, behind a small Bokeh-shaped compatibility layer in `kcwidrp/core/kcwi_plotting.py`.

- Firefox is needed for exactly one thing: `bokeh.io.export_png`, which screenshots the plot in a headless browser through Selenium and geckodriver. Bokeh has no browser-free image export, so the browser dependency cannot be removed while Bokeh writes the PNGs.
- Matplotlib writes PNGs in-process with no external binary. It is already a declared dependency (3.9.x), so the port removes five dependencies (bokeh, selenium, geckodriver, firefox, phantomjs) and adds none.
- The Bokeh surface the pipeline uses is small (six glyph methods and a handful of properties), so a thin wrapper can keep the ~48 plot call sites almost unchanged. That is what keeps the "do the new plots match the old ones" effort low: appearance is defined once, not per plot.
- Rendering from worker threads was tested and works. The live view from a worker thread is the one part not yet proven; it needs a short prototype before committing (see [Risks](#risks-and-open-questions)).

## What the pipeline uses today

| Aspect | Current state |
|---|---|
| Call sites | 16 primitives, 48 `figure()` calls, 49 `bokeh_plot()` (live) calls, 21 `save_plot()` (PNG) calls. `MakeMasterFlat` (13) and `MakeInvsens` (7) are the largest. |
| Glyphs | `line` (124), `circle` (28), `diamond` (17), `scatter` (8), `quad` (1, bias histogram), `image` (1, `MakeCube`) |
| Other features | `legend_label`, `legend.location` (7), `line_dash`, `Range1d` axis limits, `xgrid` settings (6), `gridplot` (`SubtractOverscan`), one twin x-axis (`extra_x_ranges` + `LinearAxis` in `SolveArcs`), hover tooltips (2, `MakeInvsens`) |
| Live display | Scripts launch `bokeh serve` as a subprocess; `StartBokeh` connects with `pull_session()` and opens a browser tab; `bokeh_plot()` swaps the figure and calls `session.push()` |
| Saved plots | `export_png` → Selenium → geckodriver → headless Firefox, writing `plots/<name>.png` |
| User interaction | Always terminal `input()` (accept, new offset, new wavelength limits, mask a line). The plot is only looked at, never clicked. `MakeInvsens` asks the user to hover over the spectrum to read off wavelengths. |
| Threading | The framework runs every primitive on a daemon thread named `action_loop`, not the main thread. With `want_multiprocessing = True`, several pipeline processes share one event queue. |
| Sizing | `plot_width` × `plot_height` from config (1000 × 600 px) |

Costs of the current design that the port should remove:

- About 130 lines of `kcwi_plotting.py` exist only to find and manage Firefox: snap-wrapper detection, a geckodriver fallback list, a cached global driver, and an optional pre-warm step.
- `_save_plot_firefox` calls `os.chdir(~)` around each export. The working directory is process-wide, so this is unsafe if anything else on another thread touches relative paths at that moment.
- Three config switches exist only for the browser path: `plot_firefox_compat`, `plot_prewarm_firefox`, `terminate_on_failed_bokeh_start`.
- `scripts/boklean` exists to kill leftover bokeh, geckodriver and marionette processes.
- Bokeh is pinned to 2.4.x (2021), which in turn holds `jinja2` below 3.1.

## Minimum requirements

Options should be evaluated against these. "Must" items are pass/fail.

| # | Requirement | Level |
|---|---|---|
| R1 | Writes PNGs with no web browser, WebDriver, or other external binary. Installable with `pip` alone. | Must |
| R2 | Renders correctly from a non-main thread (the `action_loop` thread), with no GUI event loop required on the main thread. | Must |
| R3 | Safe with several pipeline processes on one host (`want_multiprocessing`): no shared global state, no port or temp-file collisions, no orphaned helper processes. | Must |
| R4 | Works headless (no `DISPLAY`), for RTI/KOA-style unattended runs. | Must |
| R5 | Covers every feature in the table above: the six glyph types, legends with placement, dash styles, explicit axis limits, a twin x-axis, a grid of subplots, an image with a palette. | Must |
| R6 | Saved plots keep their file names, PNG format, and pixel size (`plot_width` × `plot_height`), so anything that collects `plots/*.png` is unaffected. | Must |
| R7 | Live view that updates while the pipeline is blocked in `input()` or sleeping for `plot_pause`, viewable remotely (browser, VNC or SSH port-forward). | Must |
| R8 | Live view supports zoom/pan and a cursor coordinate readout, which the interactive steps in `MakeInvsens` and `ReadAtlas` rely on. | Must for `plot_level` ≥ 2 |
| R9 | Low port effort: call sites change mechanically, and visual matching is controlled in one place rather than per plot. | Should |
| R10 | Plot failures never stop a reduction: a failed render or missing viewer logs a warning and continues. | Should |
| R11 | Rendering is no slower than today; well under 1 s per plot. | Should |
| R12 | Actively maintained, compatible with Python ≥ 3.12, and does not force pins on other packages. | Should |
| R13 | No net new heavy dependencies. | Should |

## Options considered

Only Matplotlib was exercised in this repository's environment. The assessments of the other libraries come from their documented behaviour, not from tests here.

| Option | R1 no browser | R2 threads | R5 features | R7/R8 live view | R9 port effort | Verdict |
|---|---|---|---|---|---|---|
| **Matplotlib** (Agg + WebAgg) | Yes | Yes (tested) | Yes | Yes, via WebAgg in a browser tab; needs prototype | Low with a wrapper | **Recommended** |
| Stay on Bokeh, newer version | No. Export still needs Selenium and a browser in every version | n/a | Yes | Yes | None | Fails R1 |
| Bokeh with Chrome instead of Firefox | No. Swaps one browser for another | n/a | Yes | Yes | Very low | Fails R1 |
| Plotly + Kaleido | No in practice. Kaleido 1.x requires an installed Chrome; the old 0.2 line bundled Chromium but is unmaintained | Export runs in a helper process | Yes | Needs Dash or regenerated HTML; no push-to-session equivalent | High (different API) | Fails R1 |
| Altair + vl-convert | Yes | Likely | Partial. Declarative grammar is a poor fit for imperative overlays and the image plot | No built-in live view | High | Fails R9 |
| PyQtGraph / VisPy | Yes | No. Qt requires the main thread, which the framework owns | Yes | Native window only; poor over SSH | High | Fails R2, R4 |
| HoloViews / hvPlot | Only through its Matplotlib backend | Same as Matplotlib | Yes | Through Bokeh or Panel | High, plus a new dependency | No gain over Matplotlib |

## Recommended design

### 1. PNG output: Matplotlib Agg, object-oriented API only

Build each plot with `matplotlib.figure.Figure` and `FigureCanvasAgg` directly, never `pyplot`. `pyplot` keeps a process-global figure registry and a "current figure", which is where Matplotlib's threading problems come from. A `Figure` that is created, drawn and saved on one thread has no shared state.

Tested in the `kcwidrp` conda environment (Python 3.12, Matplotlib 3.9.4, macOS): four threads concurrently rendering 1000 × 600 plots with lines, markers, a legend, a grid and a second x-axis produced all 20 PNGs correctly, at about 0.17 s per plot. The current Selenium export time was not measured because Firefox is not installed on the test machine.

Pixel size maps exactly: `Figure(figsize=(plot_width/100, plot_height/100), dpi=100)`.

### 2. A Bokeh-shaped wrapper, so call sites barely change

Add a small `figure()` to `kcwidrp/core/kcwi_plotting.py` that accepts the arguments the primitives already pass and returns an object exposing only what they use:

| Bokeh call used today | Matplotlib equivalent inside the wrapper |
|---|---|
| `figure(title, x_axis_label, y_axis_label, plot_width, plot_height)` | `Figure` + `set_title` / `set_xlabel` / `set_ylabel` |
| `p.line(x, y, color/line_color, line_width, line_dash, alpha, legend_label)` | `ax.plot(...)` |
| `p.circle`, `p.diamond`, `p.scatter(marker='x')` | `ax.plot(..., linestyle='none', marker='o' / 'D' / 'x')` |
| `p.quad(top, bottom, left, right)` | `ax.bar(..., align='edge')` |
| `p.image([arr], x, y, dw, dh, palette)` | `ax.imshow(arr, extent=..., origin='lower', cmap='Spectral')` |
| `p.x_range = Range1d(a, b)`, `p.y_range.start = 0`, `set_plot_lims` | `ax.set_xlim` / `ax.set_ylim` |
| `p.legend.location = "top_left"` | `ax.legend(loc='upper left')` |
| `p.extra_x_ranges` + `LinearAxis(..., 'above')` | `ax.twiny()` |
| `gridplot(plts, ncols)` | `Figure.subplots(nrows, ncols)` |
| `ColumnDataSource` + `tooltips` | Dropped; replaced by the viewer's cursor readout |

Primitives then change their import line (`from kcwidrp.core.kcwi_plotting import figure`) and little else. `bokeh_plot(p, session)` and `save_plot(p, filename)` keep their signatures; the `session` argument becomes the viewer handle.

Trade-off: the wrapper is a small private API to maintain, and it is deliberately not idiomatic Matplotlib. The alternative, rewriting roughly 200 glyph calls directly as Matplotlib, gives cleaner code but means every plot has to be visually re-checked individually. Given the requirement to minimise matching effort, the wrapper is the better first step; individual plots can be converted to plain Matplotlib later with no deadline.

### 3. Matching the old plots

Put Bokeh's look in the wrapper once:

- **Default colour.** Bokeh draws every glyph without an explicit colour in the same blue, `#1f77b4`. Matplotlib would cycle through colours. The wrapper must default to that blue with no cycling, or multi-line plots will change colour.
- **Legend.** Default position top-right, entries merged by `legend_label` (several primitives reuse one label across glyphs, for example `'Arc'` in `SolveArcs`).
- **Markers.** Bokeh's default marker size is 4 screen pixels; set the equivalent Matplotlib marker size.
- **Title and grid.** Left-aligned bold title, light grey grid on by default.
- **Dash names.** `'dashed'`, `'dashdot'`, `'dotted'` map by name.

Expected, accepted differences: fonts and anti-aliasing will not be pixel-identical, and the Bokeh toolbar strip will no longer appear in saved PNGs.

Verification, in order of cost:

1. **Reference set.** Before the port, run one blue and one red standard reduction at `plot_level = 3` on a machine that still has Firefox, and keep the `plots/` directory.
2. **Side-by-side sheet.** After the port, run the same data and generate one HTML page pairing old and new PNGs by file name. Review once against a checklist: same file name, same data, same axis limits, same labels and title, same legend entries, same colours.
3. **Structural tests.** Add pytest cases that drive the wrapper with synthetic data and assert on the resulting axes (number of lines, limits, labels, legend text, image size in pixels). These replace `test_bokeh_plotting.py` and run with no display or browser.

Pixel-comparison tests are not recommended; they are brittle across Matplotlib and FreeType versions.

### 4. Live view: WebAgg in a background thread

Matplotlib's native windows (Tk, Qt, macOS) must run on the main thread, which the framework occupies, so they are not usable from `action_loop`. WebAgg serves the figure to a browser tab over Tornado, which is the same delivery model users have today, including over an SSH port-forward. It provides zoom, pan and a live x/y cursor readout, which covers the hover-to-read-wavelength step in `MakeInvsens`.

Design:

- A renamed `StartBokeh` primitive (`StartPlotViewer`) starts the Tornado loop on its own daemon thread and opens the tab. No subprocess is launched, so the scripts drop their `bokeh serve` `Popen` calls and `boklean` becomes unnecessary.
- `bokeh_plot(p, viewer)` hands a fully built figure to the viewer under a lock. Primitives already build the whole figure before calling `bokeh_plot`, so the hand-off point exists.
- If the viewer cannot start, log a warning and carry on saving PNGs (R10).
- Under `want_multiprocessing`, each process either takes the next free port or only the process that owns the terminal starts a viewer.

Fallback if the prototype shows WebAgg is unreliable from a worker thread: a standard-library HTTP thread that serves the most recent PNG on an auto-refreshing page. That is trivially thread-safe and meets R7, but not R8, so interactive steps would lose zoom and cursor readout.

### 5. Clean-up that comes with the port

- **Dependencies** (`pyproject.toml`, `requirements.txt`, `environment.yml`, `conda_build_files/meta.yaml`): remove `bokeh`, `selenium`, `phantomjs`, `geckodriver`, `firefox`. Add `tornado` explicitly, since it is currently installed only because Bokeh pulls it in. Re-check whether the `jinja2 < 3.1` pin is still needed.
- **Config**: remove `plot_firefox_compat`, `plot_prewarm_firefox`, `terminate_on_failed_bokeh_start`. Keep `enable_bokeh` working as a deprecated alias for a new `enable_plots` so existing user configs do not break.
- **Code**: delete the Firefox and geckodriver discovery code, the cached driver, and the `os.chdir` workaround. Keep the `start_bokeh` event name as an alias in both pipeline event tables.
- **Docs**: update the installation and plotting sections that tell users to install Firefox and geckodriver.

## Migration plan

| Phase | Work | Rough size |
|---|---|---|
| 0 | Prototype WebAgg from the `action_loop` thread, including a blocked `input()` and two concurrent processes. Capture the reference PNG set. | 1–2 days |
| 1 | Write the wrapper and structural tests. Port PNG output first, behind a config switch, leaving Bokeh live display in place. This alone removes Firefox. | 2–3 days |
| 2 | Port the live view; convert `StartBokeh`; update the three launch scripts. | 2–3 days |
| 3 | Side-by-side review against the reference set; fix mismatches in the wrapper. | 1–2 days |
| 4 | Remove Bokeh, Selenium and the config switches; update dependencies and docs. | 1 day |

Sizes are estimates, not measurements. Phase 1 on its own meets the stated goal of removing Firefox, so it can ship first if phase 2 slips. It works because the wrapper can produce a Matplotlib figure for saving while the existing Bokeh figure is still used for display during the transition, at the cost of building those 21 saved plots twice for a release.

## Risks and open questions

- **WebAgg off the main thread is unproven here.** Matplotlib is not thread-safe, and with WebAgg the Tornado thread draws a figure that the action thread built. The hand-off-under-lock design should be sufficient, but this is the main technical risk and the reason for phase 0.
- **Downstream consumers of `plots/*.png`.** It is assumed, not confirmed, that RTI/KOA ingestion depends only on file names and format. Worth confirming with whoever owns that ingestion before phase 3.
- **Remote viewing.** WebAgg binds to localhost by default, like `bokeh serve`. Confirm how observers reach the plots today (VNC or port-forward) so the same route keeps working.
- **Hover tooltips.** The two `MakeInvsens` tooltips become a cursor readout in the viewer's toolbar rather than a pop-up at the data point. The printed prompt ("Hover the cursor over spectrum...") stays accurate.

## Unrelated issue found while surveying

`kcwidrp/scripts/smart_reduce_kcwi.py` imports `check_bokeh_server` from `kcwidrp.core.bokeh_plotting`, but that function does not exist there (the module defines `check_running_process`). The script will fail at import as written. The port touches this file anyway and should fix it.
