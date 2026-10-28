"""Config schema (plan T.6): every key in the shipped config files has the
type and allowed values stated in SCHEMA below.

The files are parsed with the loader the pipeline uses,
``ConfigClass(path, default_section='KCWI')`` (reduce_kcwi.py, kcwi_rti.py).
That loader treats the ``[KCWI]`` section as the default section, lower-cases
keys, and converts values with ``ConfigClass._getType`` (``True``/``False``/
``None``, then ``eval``, then int, float, str). The channel sections
``[BLUE]``/``[RED]`` are read by the scripts as raw strings and cast there
(``int(...)``, ``float(...)``), so their values are checked after
``_getType`` as well.

The key list of each file comes from a plain ``configparser`` read of the same
file (the framework's ConfigClass mixes its own defaults into ``properties``).
A key that is in a config file but not in SCHEMA fails the test: add it to
the table with its type, allowed values and a source note.

Every SCHEMA entry carries a ``source``: either the comment next to the key in
kcwi.cfg ("cfg comment") or the code that restricts the value ("code
path:function"). Allowed values and ranges that appear only in a source note
as "documented" are not enforced by the code.
"""
import configparser

import pytest
from keckdrpframework.config.framework_config import ConfigClass

from kcwidrp.core.kcwi_pkg_resources import get_resource_path

CONFIG_FILES = ["configs/kcwi.cfg", "configs/kcwi_koarti.cfg"]
CHANNEL_SECTIONS = ["BLUE", "RED"]


def _entry(kind, source, choices=None, lo=None, hi=None, lo_open=False):
    return dict(kind=kind, source=source, choices=choices, lo=lo, hi=hi,
                lo_open=lo_open)


def _bool(source):
    return _entry("bool", source)


def _int(source, **kw):
    return _entry("int", source, **kw)


def _num(source, **kw):
    return _entry("num", source, **kw)


def _str(source, **kw):
    return _entry("str", source, **kw)


# Keys are lower case, as ConfigClass stores them.
SCHEMA = {
    # --- switches -------------------------------------------------------
    "skipscat": _bool("cfg comment 'Skip subtracting scattered light?'; "
                      "used as a truth value in SubtractScatteredLight.py"),
    "skipsky": _bool("cfg comment 'Skip sky subtraction?'; reduce_kcwi.py "
                     "sets it from -k (a bool)"),
    "air_to_vacuum": _bool("cfg comment 'Correct to vacuum wavelengths?'; "
                           "truth value in WavelengthCorrections.py"),
    "enable_plot_viewer": _bool("cfg comment 'Show plots in the live plot "
                                "viewer'; plot_config.viewer_requested()"),
    "plot_viewer_open_browser": _bool("cfg comment; StartPlotViewer.py "
                                      "truth value"),
    "saveintims": _bool("cfg comment 'Save intermediate images?'; truth "
                        "value in TrimOverscan.py"),
    "clobber": _bool("cfg comment 'Overwrite outputs?'; truth value in the "
                     "pipelines' action_planner and TraceBars.py"),
    # --- wavelength correction -------------------------------------------
    "radial_velocity_correction": _str(
        "cfg comment 'Options are [heliocentric, barycentric, none]'; "
        "WavelengthCorrections.py options list",
        choices={"heliocentric", "barycentric", "none"}),
    # --- plotting -------------------------------------------------------
    "plot_level": _int(
        "cfg comment 'Plot level: 0 - no plots nor interaction ... 3'; "
        "primitives compare it with >= 1, 2, 3",
        choices={0, 1, 2, 3}),
    "plot_viewer_port": _int("cfg comment; first port the viewer tries "
                             "(then the next 9)", lo=1, hi=65526),
    "plot_pause": _num("cfg comment 'pause between plots (seconds)'; "
                       "passed to time.sleep()", lo=0),
    "verbose": _int("cfg comment 'Wavelength fitting verbosity: 1 - normal, "
                    "2 - verbose'; SolveArcs.py and GetAtlasLines.py test "
                    "`> 1`", choices={1, 2}),
    "plot_width": _int("cfg comment 'Bokeh plot dimensions in pixels'; "
                       "passed as plot_width to the plotting code", lo=1),
    "plot_height": _int("cfg comment 'Bokeh plot dimensions in pixels'; "
                        "passed as plot_height to the plotting code", lo=1),
    # --- output ---------------------------------------------------------
    "output_directory": _str("cfg: parameters that should be the same for "
                             "both channels; join()ed into paths, must be "
                             "non-empty", choices=None),
    "procfile": _str("cfg: 'procfile = kcwib.proc' (a bare file name; "
                     "_getType leaves it a str); Proctab.read_proctab "
                     "opens it"),
    # --- wavelength fitting ----------------------------------------------
    "default_arc_lamp": _str(
        "cfg comment 'choices are ThAr and FeAr'; compared with the "
        "frame's calibration lamp in ProcessArc.py",
        choices={"ThAr", "FeAr"}),
    "taperfrac": _num("cfg comment '(0. - 1.)'", lo=0.0, hi=1.0),
    "midfrac": _num("cfg comment '-1 for default of 0.33'; ReadAtlas.py "
                    "uses it only if 0 < mf <= 1, else its default. "
                    "Allowed: -1 or (0, 1]", lo=-1.0, hi=1.0),
    "atoff": _int("cfg comment 'pixels to offset ref bar (0 - for "
                  "default)'; -o is parsed as int in reduce_kcwi.py"),
    "linelist": _str("cfg comment 'Optional line list ... instead of Atlas "
                     "generated'; empty string means none"),
    "linethresh": _num("cfg comment 'Line threshold in electrons for "
                       "fitting'", lo=0.0),
    "tukeyalpha": _num("cfg comment 'Tukey alpha value'; scipy tukey window "
                       "alpha is in [0, 1]", lo=0.0, hi=1.0),
    "fracmax": _num("cfg comment 'Fraction of max of a given line'; a "
                    "fraction in (0, 1]", lo=0.0, hi=1.0, lo_open=True),
    # --- instrument constants ---------------------------------------------
    "nbars": _int("cfg: number of bars, 120 for KCWI", lo=1),
    "refbar": _int("cfg: reference bar index, must be < NBARS", lo=0),
    "pixscale": _num("cfg comment 'unbinned pixel size in degrees'",
                     lo=0.0, lo_open=True),
    "slicescale": _num("cfg comment 'Large slicer slice size in degrees'",
                       lo=0.0, lo_open=True),
    "pix": _num("cfg comment 'pixel size in mm'", lo=0.0, lo_open=True),
    "fcam": _num("cfg comment 'FOCAL LENGTH OF THE CAMERA'", lo=0.0,
                 lo_open=True),
    "gamma": _num("cfg comment 'mean out-of-plane angle for diffraction "
                  "(deg)'"),
    "wavefid": _num("cfg comment 'Fiducial wavelength for wavelength bins'",
                    lo=0.0, lo_open=True),
    "rotoff": _num("cfg comment 'Rotator/IFU offset angle in degrees'"),
    "knotspp": _num("cfg comment 'Knots per Y image pixel for sky "
                    "modeling'", lo=0.0, lo_open=True),
    "psfwid": _num("cfg comment 'Nominal window for pt. source (unbinned "
                   "px)'; MakeInvsens.py int()s psfwid / xbinsize", lo=1),
    # --- minimum frame counts ---------------------------------------------
    "bias_min_nframes": _int("code: min_files in action_planner", lo=1),
    "flat_min_nframes": _int("code: min_files in action_planner", lo=1),
    "dome_min_nframes": _int("code: min_files in action_planner", lo=1),
    "twiflat_min_nframes": _int("code: min_files in action_planner", lo=1),
    "dark_min_nframes": _int("code: min_files in action_planner", lo=1),
    "arc_min_nframes": _int("code: min_files in action_planner; scripts "
                            "int() the [BLUE]/[RED] values", lo=1),
    "contbars_min_nframes": _int("code: min_files in action_planner; "
                                 "scripts int() the [BLUE]/[RED] values",
                                 lo=1),
    "object_min_nframes": _int("code: min_files in action_planner; scripts "
                               "int() the [BLUE]/[RED] values", lo=1),
    "minoscanpix": _int("code: scripts int() the [BLUE]/[RED] values; "
                        "SubtractOverscan.py", lo=1),
    "oscanbuf": _int("code: scripts int() the [BLUE]/[RED] values; "
                     "SubtractOverscan.py", lo=0),
    # --- site (radial velocity correction) --------------------------------
    "latitude": _num("cfg comment 'Keck II Location'; degrees",
                     lo=-90.0, hi=90.0),
    "longitude": _num("cfg comment 'Keck II Location'; degrees",
                      lo=-360.0, hi=360.0),
    "altitude": _num("cfg comment 'Keck II Location'; metres"),
    # --- cosmic ray rejection (astroscrappy.detect_cosmics) -----------------
    "crr_minexptime": _num("cfg: frames shorter than this are not cleaned "
                           "(RemoveCosmicRays.py `exptime >= ...`)", lo=0.0),
    "crr_pssl": _num("cfg: astroscrappy pssl, previously subtracted sky "
                     "level", lo=0.0),
    "crr_gain": _num("cfg: astroscrappy gain", lo=0.0, lo_open=True),
    "crr_readnoise": _num("cfg: astroscrappy readnoise", lo=0.0),
    "crr_sigclip": _num("cfg: astroscrappy sigclip", lo=0.0, lo_open=True),
    "crr_sigfrac": _num("cfg: astroscrappy sigfrac", lo=0.0, lo_open=True),
    "crr_objlim": _num("cfg: astroscrappy objlim", lo=0.0, lo_open=True),
    "crr_psffwhm": _num("cfg: astroscrappy psffwhm", lo=0.0, lo_open=True),
    "crr_fsmode": _str("astroscrappy documented fsmode values",
                       choices={"median", "convolve"}),
    "crr_psfmodel": _str("astroscrappy documented psfmodel values",
                         choices={"gauss", "gaussx", "gaussy", "moffat"}),
    "crr_satlevel": _num("cfg: astroscrappy satlevel (counts)", lo=0.0,
                         lo_open=True),
    "crr_verbose": _bool("cfg: astroscrappy verbose"),
    "crr_sepmed": _bool("cfg: astroscrappy sepmed"),
    "crr_cleantype": _str("astroscrappy documented cleantype values",
                          choices={"median", "medmask", "meanmask", "idw"}),
    "crr_niter": _int("cfg: astroscrappy niter", lo=1),
    # --- present only in kcwi_koarti.cfg ------------------------------------
    "inter": _int("kcwi_koarti.cfg only ('inter = 1', next to the "
                  "commented-out 'interactive = 1'); nothing reads it"),
}


def _plain_parse(relpath):
    """Raw configparser read; the [KCWI] section is the default section."""
    path = get_resource_path("kcwidrp", relpath)
    parser = configparser.ConfigParser(default_section="KCWI")
    parser.read(path)
    return parser


def _own_keys(relpath):
    """Return {section: [keys]} the file itself defines (lower case)."""
    parser = _plain_parse(relpath)
    out = {"KCWI": list(parser.defaults())}
    for sec in parser.sections():
        out[sec] = list(parser._sections[sec])
    return out


def _loaded(relpath):
    path = get_resource_path("kcwidrp", relpath)
    return ConfigClass(path, default_section="KCWI")


def _value(config, section, key):
    """The value the pipeline sees for key in the given section."""
    if section == "KCWI":
        return config.properties[key]
    # Channel sections are read raw by the scripts and cast there.
    return config._getType(dict(config.items(section))[key])


def _check(key, value):
    spec = SCHEMA[key]
    kind = spec["kind"]
    ok_type = {
        "bool": lambda v: isinstance(v, bool),
        "int": lambda v: isinstance(v, int) and not isinstance(v, bool),
        "num": lambda v: (isinstance(v, (int, float))
                          and not isinstance(v, bool)),
        "str": lambda v: isinstance(v, str),
    }[kind]
    problems = []
    if not ok_type(value):
        problems.append(f"{key} = {value!r} has type {type(value).__name__},"
                        f" schema says {kind} ({spec['source']})")
        return problems
    if spec["choices"] is not None and value not in spec["choices"]:
        problems.append(f"{key} = {value!r} is not one of "
                        f"{sorted(spec['choices'], key=str)} "
                        f"({spec['source']})")
    lo, hi = spec["lo"], spec["hi"]
    if lo is not None:
        if value < lo or (spec["lo_open"] and value == lo):
            problems.append(f"{key} = {value!r} is below the allowed "
                            f"minimum {lo}{' (exclusive)' if spec['lo_open'] else ''}"
                            f" ({spec['source']})")
    if hi is not None and value > hi:
        problems.append(f"{key} = {value!r} is above the allowed maximum "
                        f"{hi} ({spec['source']})")
    if key == "output_directory" and not value:
        problems.append("output_directory must be non-empty")
    return problems


def _cases():
    cases = []
    for relpath in CONFIG_FILES:
        for section in ["KCWI"] + CHANNEL_SECTIONS:
            if section != "KCWI" and relpath != "configs/kcwi.cfg":
                continue
            cases.append(pytest.param(relpath, section,
                                      id=f"{relpath.split('/')[-1]}-{section}"))
    return cases


@pytest.mark.parametrize("relpath, section", _cases())
def test_every_key_matches_schema(relpath, section):
    config = _loaded(relpath)
    keys = _own_keys(relpath)[section]
    assert keys, f"{relpath} [{section}] has no keys"
    unknown = [k for k in keys if k not in SCHEMA]
    # source: specification (plan T.6): every key has an entry in SCHEMA
    assert not unknown, (
        f"{relpath} [{section}]: keys with no SCHEMA entry (add type, "
        f"allowed values and source): {unknown}")
    problems = []
    for key in keys:
        problems.extend(_check(key, _value(config, section, key)))
    # source: specification (plan T.6): type and allowed values per key
    assert not problems, "\n".join(problems)


def test_channel_sections_exist_in_main_config():
    # reduce_kcwi.py reads kcwi_config.BLUE[...] / .RED[...] for -b / -r.
    parser = _plain_parse("configs/kcwi.cfg")
    # source: reduce_kcwi.py main(), channel specific parameters
    assert parser.sections() == CHANNEL_SECTIONS
    needed = ["arc_min_nframes", "contbars_min_nframes", "object_min_nframes",
              "minoscanpix", "oscanbuf", "procfile", "linethresh",
              "tukeyalpha", "fracmax"]
    keys = _own_keys("configs/kcwi.cfg")
    for sec in CHANNEL_SECTIONS:
        missing = [k for k in needed if k not in keys[sec]]
        assert not missing, f"[{sec}] lacks keys the scripts read: {missing}"


def test_schema_has_no_stale_entries():
    present = set()
    for relpath in CONFIG_FILES:
        for keys in _own_keys(relpath).values():
            present.update(keys)
    stale = sorted(set(SCHEMA) - present)
    # source: specification: the table describes the shipped files
    assert not stale, f"SCHEMA entries for keys in no config file: {stale}"


def test_schema_rejects_string_plot_level():
    # The schema must tell the int 1 from the string "1" (plan T.6
    # validation: `plot_level = "1"` fails). Checked on the checker itself.
    assert _check("plot_level", 1) == []
    assert _check("plot_level", "1")
    assert _check("plot_level", True)
    assert _check("plot_level", 4)
    assert _check("radial_velocity_correction", "geocentric")
