"""Plot colors by role (docs/interactive_prompts.md, "Plot color scheme").

The primitives name a role here instead of a color, so a scheme change is
one edit::

    from kcwidrp.core import kcwi_plot_colors as pc
    p.line(x, y, color=pc.DATA)

Nothing enforces this: any color the wrapper accepts may still be given
directly, and a plot may use its own colors where no role fits (the P11
rejection classes, the P28-P32 slice series).  The help texts name colors
in words and are not derived from these values; change them by hand.
"""

DATA = 'blue'               # data being inspected
SECOND = 'cyan'             # second data series of the same kind
REFERENCE = 'black'         # atlas, REFBAR arc, reference spectrum or slice
FIT = 'red'                 # fit, model, smoothing
MEAN = 'red'                # mean (solid) and +-sigma, +-RMS (dashed)
SUBSET_MEAN = 'black'       # mean of a subset, e.g. per slice (thick)
CWAVE = 'orange'            # central wavelength (dash-dot-dot)
WAVGOOD = 'green'           # good wavelength range
WAVALL = 'gray'             # full wavelength range
FOUND = 'magenta'           # found or adopted value (peak, centroid)
PREDICTED = 'black'         # predicted or calculated value (dashed)
LIMIT = 'purple'            # limits and windows
LINE_MASK = 'purple'        # line masks (dashed pairs)
THRESHOLD = 'gray'          # thresholds (dashed), zero lines (dotted)
MASKED = 'gray'             # masked points, rejected before a fit (x)
OUTLIER = 'gold'            # rejected as outliers in residual clipping
SLICE_BOUNDARY = 'black'    # slice boundaries (dash-dot)
