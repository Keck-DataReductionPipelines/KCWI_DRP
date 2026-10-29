"""Run one primitive without the framework queue (task T.2; rules in docs/testing.md).

``run_primitive(PrimitiveClass, ccddata, ...)`` builds the ``action`` and
``context`` a primitive receives from ``Framework.execute``, calls
``PrimitiveClass(action, context).apply()`` and returns its output.  Tests
normally use the ``run_primitive`` fixture (``conftest.py``), which supplies
the instrument config (``kcwi_config``), the Proctab (``tmp_redux``) and runs
with the current directory set to the ``tmp_redux`` root, so files the
primitive writes (``redux/``, ``plots/``, the proc file) land in ``tmp_path``.

To chain primitives, pass the previous output (an ``Arguments``) as
``ccddata``: its attributes become the next action's args, as in the
framework.  Building args again from a processed frame would be wrong, e.g.
``map_ccd`` needs the raw ASECn/BSECn/DSECn keywords that TrimOverscan removes.

Framework objects: real or stand-in
-----------------------------------
Checked against keckdrpframework 1.0.0 (site-packages/keckdrpframework).

================  ===========================================================
Object            What the harness uses
================  ===========================================================
action            Real ``models.action.Action(event, event_info, args)``, with
                  a real ``models.event.Event`` and ``event_info =
                  (<class name>, None, None)`` (action name, next state, new
                  event), as ``Framework.event_to_action`` builds it from the
                  pipeline's event table.
action.args       ``_HarnessArguments``: a subclass of the real
                  ``models.arguments.Arguments`` that only adds
                  ``__getattr__``, so reading a missing attribute raises
                  ``HarnessAttributeError`` naming the attribute and the
                  primitive.  Filled as ``ingest_file._perform`` fills its
                  ``out_args`` (see ARGS_PROVIDED), by calling the real
                  ``ingest_file`` helper methods on ``ccddata``.
context           ``_HarnessContext``: a subclass of the real
                  ``models.processing_context.ProcessingContext`` (built with
                  two real ``core.queues.SimpleEventQueue`` objects, so
                  ``push_event``/``append_event`` work and the pushed events
                  can be inspected), adding only ``__getattr__`` as above.
                  Attributes are set as ``scripts/reduce_kcwi.py`` sets them.
context.config    Real ``config.framework_config.ConfigClass`` read from
                  ``kcwidrp/configs/framework.cfg`` (as ``Framework.__init__``
                  reads it), with ``.instrument`` = the ``kcwi.cfg``
                  ConfigClass.  Note: ``ConfigClass.__getattr__`` returns
                  ``None`` for an unknown key; that is the real behaviour and
                  the harness keeps it.
context.logger    ``logging.getLogger("DRPF")`` and ``logging.getLogger(
pipeline_logger   "KCWI")``: the logger names reduce_kcwi uses, but without
                  ``logging.config.fileConfig(logger.cfg)``, which would add
                  file handlers under ``logs/``.  Records reach ``caplog``.
================  ===========================================================

Why subclasses rather than plain real objects: the real ``Arguments`` and
``ProcessingContext`` raise a bare ``AttributeError`` for a missing
attribute, deep inside the primitive.  The subclasses change nothing else.
``HarnessAttributeError`` is an ``AttributeError`` subclass, so ``getattr(x,
name, default)`` and ``hasattr`` keep working.  No stand-in is hand-written:
``Framework`` itself is not built because its constructor needs a pipeline
module, a logger config and installs signal handlers; ``BasePrimitive`` reads
only ``action.args``, ``context.logger`` and ``context.config``.

What the harness does around ``apply()``
----------------------------------------
* ``config.instrument``: plot_level = 0 and plot_pause = 0 by default; the
  channel-specific values (``arc_min_nframes``, ``contbars_min_nframes``,
  ``object_min_nframes``, ``minoscanpix``, ``oscanbuf``, ``procfile``) are
  copied from the ``[BLUE]`` or ``[RED]`` section chosen by the frame's
  CAMERA keyword, as ``reduce_kcwi -b/-r`` does; ``procfile`` then points into
  the temporary redux directory.  ``config_overrides`` are applied last and
  win over all of these.
* ``builtins.input`` is replaced: with ``answers=None`` any call fails the
  test, naming the prompt; with a list, answers are returned in order, a call
  after the list is used up fails ("too few answers"), and answers left over
  when ``apply()`` returns fail the test ("too many answers").
* ``time.sleep`` is a no-op (calls are recorded in ``RunRecord.sleeps``).
* Exceptions from the primitive propagate unchanged: ``BasePrimitive.apply``
  logs and re-raises them; ``Framework.execute`` would catch them, the harness
  does not.

Attributes primitives read
--------------------------
``scan_primitive_attributes()`` lists every ``self.context.<attr>`` and
``self.action.args.<attr>`` read in ``kcwidrp/primitives`` (see its docstring
for the patterns covered).  CONTEXT_PROVIDED / ARGS_PROVIDED list what the
harness sets; CONTEXT_NOT_PROVIDED / ARGS_NOT_PROVIDED list every other
attribute found, with the step that sets it in a real reduction.  A test
supplies those with ``context_extra=`` or ``action_extra=``.
``unit/test_harness.py`` fails when the scan finds an attribute in neither
list.
"""
import ast
import builtins
import logging
import pathlib
import time
import warnings
from dataclasses import dataclass, field
from unittest import mock

import pytest
from keckdrpframework.config.framework_config import ConfigClass
from keckdrpframework.core.queues import SimpleEventQueue
from keckdrpframework.models.action import Action
from keckdrpframework.models.arguments import Arguments
from keckdrpframework.models.event import Event
from keckdrpframework.models.processing_context import ProcessingContext

from kcwidrp.core.kcwi_pkg_resources import get_resource_path

PRIMITIVES_DIR = pathlib.Path(__file__).resolve().parent.parent / "primitives"

# ---------------------------------------------------------------------------
# What the harness provides
# ---------------------------------------------------------------------------

# Set on every context.  The first five are what reduce_kcwi.py sets (D11/T.11:
# bokeh_session = None, as init_plot_session does with the viewer off); the
# rest are set by the real ProcessingContext constructor.
CONTEXT_PROVIDED = frozenset({
    "config", "logger", "pipeline_logger", "proctab", "bokeh_session",
    "name", "state", "event_queue", "event_queue_hi", "data_set", "debug",
})

# Context attributes primitives read that a single-primitive run does not
# have: they are set by an earlier primitive in the arc/geometry chain.
CONTEXT_NOT_PROVIDED = {
    "arcs": "ExtractArcs",
    "bars": "ExtractArcs",
    "trace": "TraceBars / ExtractArcs",
    "bar_offsets": "ArcOffsets",
    "prelim_disp": "CalcPrelimDisp",
}

# The args ingest_file._perform puts in out_args, in its order, with the
# ingest_file method that computes each one (None = computed inline below).
# unit/test_harness.py checks this list against ingest_file._perform.
INGEST_ARGS = (
    ("in_proctab", None), ("last_suffix", None), ("ccddata", None),
    ("table", None), ("name", None), ("imtype", None), ("groupid", None),
    ("camera", "camera"), ("dich", "dich"), ("camangle", "camang"),
    ("filter", "filter"), ("grangle", "grangle"), ("grating", "grating"),
    ("adjang", "adjang"), ("rho", "rho"), ("cwave", "cwave"),
    ("resolution", "resolution"), ("atsig", "atsig"),
    ("dwout", "delta_wave_out"), ("namps", "namps"), ("nasmask", "nasmask"),
    ("shufrows", "shufrows"), ("numopen", "numopen"), ("ampmode", "ampmode"),
    ("xbinsize", None), ("ybinsize", None), ("ifunum", None),
    ("ifuname", None), ("plotlabel", "plotlabel"), ("stdlabel", "stdlabel"),
    ("illum", "illum"), ("map_ccd", None),
    ("calibration_lamp", "calibration_lamp"), ("ttime", None),
)

ARGS_PROVIDED = frozenset(name for name, _ in INGEST_ARGS)

# Args primitives read that ingest does not set, with what sets them in a
# real reduction.  A test passes them with action_extra=.
_PIPELINE = "kcwi_pipeline action (e.g. action_bias)"
ARGS_NOT_PROVIDED = {
    # set by kcwi_pipeline before a stacking primitive
    "min_files": _PIPELINE,
    "new_type": _PIPELINE,
    "stack_type": _PIPELINE,
    "want_type": _PIPELINE,
    # set by earlier primitives
    "ait_file": "SolveAIT",
    "arc_image": "ExtractArcs",
    "arc_number": "ExtractArcs",
    "at_flux": "GetAtlasLines",
    "at_wave": "GetAtlasLines",
    "atmaxrow": "GetAtlasLines",
    "atmaxwave": "GetAtlasLines",
    "atminrow": "GetAtlasLines",
    "atminwave": "GetAtlasLines",
    "av_bar_nls": "SolveArcs",
    "av_bar_sig": "SolveArcs",
    "bar_avg": "FindBars",
    "bar_id": "ExtractArcs",
    "cont_source_pos": "MakeMasterSky",
    "cont_source_width": "MakeMasterSky",
    "contbar_image": "FindBars / ExtractArcs",
    "contbar_image_number": "FindBars / ExtractArcs",
    "contsky": "MakeMasterSky",
    "cube_size": "CubeImage",
    "destination_control_points": "ExtractArcs",
    # Nothing sets it: SolveAIT reads it inside try/except AttributeError and
    # falls back to 1.0 (HarnessAttributeError is an AttributeError, so the
    # fallback still works under the harness).
    "dich_frac": "nothing (SolveAIT falls back to 1.0)",
    "dichroic_fraction": "GetAtlasLines",
    "fincoeff": "SolveArcs",
    "geometry_file": "SolveGeom / SolveAIT",
    "invsname": "FluxCalibrate",
    "master_flat": "CorrectIllumination",
    "maxrow": "ReadAtlas",
    "middle_centers": "FindBars",
    "middle_row": "FindBars",
    "minrow": "ReadAtlas",
    "offset_wave": "ReadAtlas",
    "original_filename": "ExtractArcs",
    "reference_bar_separation": "ExtractArcs",
    "reference_delta_x": "FindBars",
    "refdisp": "ReadAtlas",
    "reflux": "ReadAtlas",
    "refwave": "ReadAtlas",
    "skyfile": "MakeMasterSky / SubtractSky",
    "skymask": "MakeMasterSky / SubtractSky",
    "slice_id": "ExtractArcs",
    "source_control_points": "ExtractArcs",
    "st_bar_nls": "SolveArcs",
    "st_bar_sig": "SolveArcs",
    "stdfile": "MakeInvsens / MakeMasterSky",
    "stdname": "MakeInvsens / MakeMasterSky",
    "twkcoeff": "FitCenter",
    "wave0out": "SolveGeom / SolveAIT",
    "wave1out": "SolveGeom / SolveAIT",
    "waveall0": "SolveGeom / SolveAIT",
    "waveall1": "SolveGeom / SolveAIT",
    "wavegood0": "SolveGeom / SolveAIT",
    "wavegood1": "SolveGeom / SolveAIT",
    "wavemid": "SolveGeom / SolveAIT",
    "wavensall0": "SolveGeom / SolveAIT",
    "wavensall1": "SolveGeom / SolveAIT",
    "wavensgood0": "SolveGeom / SolveAIT",
    "wavensgood1": "SolveGeom / SolveAIT",
    "wavensmid": "SolveGeom / SolveAIT",
    "window": "FindBars",
    "x0": "ReadAtlas",
    "x0out": "SolveGeom / SolveAIT",
    "xsvals": "SolveArcs",
    "xvals": "ReadAtlas",
    # read in SendHTTP only (the RTI pipeline sets it)
    "ingest_time": "kcwi_rti pipeline (not in reduce_kcwi)",
}

# reduce_kcwi.py copies these from the [BLUE] / [RED] section of kcwi.cfg
# into config.instrument (cast to int, as it does).
CHANNEL_KEYS = ("arc_min_nframes", "contbars_min_nframes",
                "object_min_nframes", "minoscanpix", "oscanbuf")


# ---------------------------------------------------------------------------
# Errors and the framework subclasses
# ---------------------------------------------------------------------------

class HarnessAttributeError(AttributeError):
    """A primitive read a context/args attribute the harness does not set."""


def _missing(kind, owner, name, extra_kw):
    hint = ARGS_NOT_PROVIDED.get(name) if kind == "action.args" else \
        CONTEXT_NOT_PROVIDED.get(name)
    where = (" (in a reduction it is set by %s)" % hint) if hint else \
        " (not in the harness's documented lists)"
    return HarnessAttributeError(
        "run_primitive: %s read %s.%s, which the harness does not provide%s; "
        "pass it with %s={%r: ...}" % (owner, kind, name, where, extra_kw,
                                       name))


class HarnessIngestError(Exception):
    """The frame's header does not let ingest_file's helpers build an arg."""


class _HarnessArguments(Arguments):
    """Real Arguments; a missing attribute fails naming it and the primitive."""
    _harness_primitive = "<primitive>"

    def __getattr__(self, name):
        # Only called when normal lookup fails.  Dunder lookups (copy, pickle)
        # keep the normal AttributeError.
        if name.startswith("__"):
            raise AttributeError(name)
        raise _missing("action.args", type(self)._harness_primitive, name,
                       "action_extra")


class _HarnessContext(ProcessingContext):
    """Real ProcessingContext; a missing attribute fails naming it."""
    _harness_primitive = "<primitive>"

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        raise _missing("context", type(self)._harness_primitive, name,
                       "context_extra")


@dataclass
class RunRecord:
    """What happened during one run (filled when passed as ``record=``)."""
    action: object = None
    context: object = None
    prompts: list = field(default_factory=list)   # input() prompt texts
    sleeps: list = field(default_factory=list)    # time.sleep() arguments
    output: object = None


# ---------------------------------------------------------------------------
# Building action and context
# ---------------------------------------------------------------------------

def framework_config():
    """The framework ConfigClass, read as Framework.__init__ reads it."""
    return ConfigClass(get_resource_path("kcwidrp", "configs/framework.cfg"))


def apply_channel_parameters(instrument, camera, procfile=None):
    """Copy [BLUE]/[RED] values into ``instrument`` as reduce_kcwi -b/-r does.

    ``procfile`` replaces the section's proc file name when given.
    """
    section = getattr(instrument, camera.upper())
    if section is None:
        raise ValueError("kcwi.cfg has no [%s] section" % camera)
    for key in CHANNEL_KEYS:
        setattr(instrument, key, int(section[key]))
    instrument.procfile = procfile if procfile is not None \
        else section["procfile"]


def ingest_args(ccddata, context, args_class=_HarnessArguments, skip=()):
    """Return Arguments as ingest_file._perform builds them, from ``ccddata``.

    Uses the real ``ingest_file`` helper methods on ``ccddata`` (no file is
    read).  Names in ``skip`` are not computed (the caller supplies them).
    """
    from kcwidrp.primitives.kcwi_file_primitives import ingest_file
    stub = Action(None, ("ingest_file", None, None), Arguments())
    ing = ingest_file(stub, context)
    ing.ccddata = ccddata
    args = args_class()
    for name, method in INGEST_ARGS:
        if name in skip:
            continue
        try:
            value = _ingest_value(name, method, ing, ccddata, context)
        except Exception as exc:
            raise HarnessIngestError(
                "run_primitive: could not build action.args.%s from the frame "
                "header as ingest_file does (%s: %s); fix the frame or pass "
                "it with action_extra={%r: ...}"
                % (name, type(exc).__name__, exc, name)) from exc
        setattr(args, name, value)
    return args


def _ingest_value(name, method, ing, ccddata, context):
    """One out_args value, computed as ingest_file._perform computes it."""
    if method is not None:
        return getattr(ing, method)()
    if name == "in_proctab":
        return context.proctab.in_proctab(frame=ccddata)
    if name == "last_suffix":
        return context.proctab.last_suffix(frame=ccddata) \
            if context.proctab.in_proctab(frame=ccddata) else ""
    if name == "ccddata":
        return ccddata
    if name == "table":
        return None           # kcwi_fits_reader: no 'Exposure Events' HDU
    if name == "name":
        # ingest: the file name the event carried; here the OFNAME keyword.
        return ccddata.header.get("OFNAME", "harness.fits")
    if name == "imtype":
        return ing.get_keyword("IMTYPE")
    if name == "groupid":
        gid = ing.get_keyword("GROUPID")
        return "NONE" if gid is None or len(gid) <= 0 else gid
    if name in ("xbinsize", "ybinsize"):
        xb, yb = map(int, ing.get_keyword("BINNING").split(","))
        return xb if name == "xbinsize" else yb
    if name == "ifunum":
        return int(ing.get_keyword("IFUNUM"))
    if name == "ifuname":
        return ing.get_keyword("IFUNAM")
    if name == "map_ccd":
        xb, yb = map(int, ing.get_keyword("BINNING").split(","))
        return ing.map_ccd(xb, yb)
    if name == "ttime":
        return ing.get_keyword("TTIME")
    # INGEST_ARGS and this function disagree.
    raise RuntimeError("harness: no rule for ingest arg %r" % name)


def build_action_context(PrimitiveClass, ccddata, *, instrument_config,
                         proctab, procfile=None, config_overrides=None,
                         context_extra=None, action_extra=None):
    """Return (action, context) for running ``PrimitiveClass`` on ``ccddata``.

    ``ccddata`` is a frame (args are built as ingest_file builds them) or the
    Arguments a previous primitive returned (its attributes are reused).

    ``instrument_config`` is the kcwi.cfg ConfigClass (plot_level and
    plot_pause already 0, as the kcwi_config fixture makes it).
    """
    pname = PrimitiveClass.__name__
    args_class = type("_HarnessArguments", (_HarnessArguments,),
                      {"_harness_primitive": pname})
    context_class = type("_HarnessContext", (_HarnessContext,),
                         {"_harness_primitive": pname})

    overrides = {k.lower(): v for k, v in (config_overrides or {}).items()}
    instrument = instrument_config
    frame = ccddata.ccddata if isinstance(ccddata, Arguments) else ccddata
    camera = frame.header.get("CAMERA", "BLUE")
    apply_channel_parameters(instrument, camera, procfile=procfile)
    for key, value in overrides.items():
        # Attributes set as reduce_kcwi sets them shadow ConfigClass
        # properties, so an override must replace both.
        instrument.properties[key] = value
        if key in instrument.__dict__:
            setattr(instrument, key, value)

    config = framework_config()
    config.instrument = instrument

    context = context_class(SimpleEventQueue(), SimpleEventQueue(),
                            logging.getLogger("DRPF"), config)
    context.pipeline_logger = logging.getLogger("KCWI")
    context.proctab = proctab
    context.bokeh_session = None
    for key, value in (context_extra or {}).items():
        setattr(context, key, value)

    action_extra = dict(action_extra or {})
    if isinstance(ccddata, Arguments):
        # Chaining: the previous primitive's output becomes this action's
        # args, as Framework.execute passes action output to the next event.
        # Its attributes (including map_ccd, computed once at ingest from the
        # raw header) are carried over; nothing is re-ingested.
        args = args_class(*list(ccddata))
        args.__dict__.update({k: v for k, v in ccddata.__dict__.items()
                              if k != "_pos_args"})
    else:
        args = ingest_args(ccddata, context, args_class=args_class,
                           skip=set(action_extra))
    for key, value in action_extra.items():
        setattr(args, key, value)

    with warnings.catch_warnings():
        # keckdrpframework 1.0.0 Event() calls datetime.utcnow().
        warnings.simplefilter("ignore", DeprecationWarning)
        event = Event(pname, args)
    action = Action(event, (pname, None, None), args)
    return action, context


def run_primitive(PrimitiveClass, ccddata, config_overrides=None,
                  context_extra=None, answers=None, action_extra=None, *,
                  instrument_config, proctab, procfile=None, record=None):
    """Run ``PrimitiveClass(action, context).apply()`` and return its output.

    See the module docstring.  ``record`` (a RunRecord) receives the action,
    context, prompts and sleep calls of the run.
    """
    pname = PrimitiveClass.__name__
    action, context = build_action_context(
        PrimitiveClass, ccddata, instrument_config=instrument_config,
        proctab=proctab, procfile=procfile,
        config_overrides=config_overrides, context_extra=context_extra,
        action_extra=action_extra)
    rec = record if record is not None else RunRecord()
    rec.action, rec.context = action, context
    queue = None if answers is None else list(answers)
    n_given = 0 if answers is None else len(answers)

    def fake_input(prompt=""):
        rec.prompts.append(prompt)
        if queue is None:
            pytest.fail("run_primitive: %s called input(%r), but the test "
                        "gave no answers (pass answers=[...])"
                        % (pname, prompt), pytrace=False)
        if not queue:
            pytest.fail("run_primitive: too few answers: %s called input(%r) "
                        "after all %d answer(s) were used"
                        % (pname, prompt, n_given), pytrace=False)
        return queue.pop(0)

    def fake_sleep(seconds):
        rec.sleeps.append(seconds)

    with mock.patch.object(builtins, "input", fake_input), \
            mock.patch.object(time, "sleep", fake_sleep):
        output = PrimitiveClass(action, context).apply()
    rec.output = output
    if queue:
        pytest.fail("run_primitive: too many answers: %s returned with %d of "
                    "%d answer(s) unused: %r" % (pname, len(queue), n_given,
                                                 queue), pytrace=False)
    return output


# ---------------------------------------------------------------------------
# Attribute scan
# ---------------------------------------------------------------------------

def _chain(node):
    """['self', 'action', 'args', 'x'] for self.action.args.x, else []."""
    parts = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if not isinstance(node, ast.Name):
        return []
    parts.append(node.id)
    return parts[::-1]


def _is_args(chain, args_aliases, action_aliases):
    return (chain[:3] == ["self", "action", "args"] and len(chain) == 3) or \
        (len(chain) == 2 and chain[0] in action_aliases
         and chain[1] == "args") or \
        (len(chain) == 1 and chain[0] in args_aliases)


def _is_context(chain, context_aliases):
    return (chain == ["self", "context"]) or \
        (len(chain) == 1 and chain[0] in context_aliases)


def _reads_in_function(func, args_aliases, context_aliases, action_aliases):
    """Yield ('args'|'context', attr) for loads in one function body."""
    args_aliases = set(args_aliases)
    context_aliases = set(context_aliases)
    # Local aliases: x = self.action.args / x = self.context
    for node in ast.walk(func):
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and \
                isinstance(node.targets[0], ast.Name):
            ch = _chain(node.value)
            if _is_args(ch, args_aliases, action_aliases):
                args_aliases.add(node.targets[0].id)
            elif _is_context(ch, context_aliases):
                context_aliases.add(node.targets[0].id)
    aug_targets = {id(n.target) for n in ast.walk(func)
                   if isinstance(n, ast.AugAssign)}
    for node in ast.walk(func):
        if not isinstance(node, ast.Attribute):
            continue
        if not (isinstance(node.ctx, ast.Load) or id(node) in aug_targets):
            continue
        base = _chain(node.value)
        if not base:
            continue
        if _is_args(base, args_aliases, action_aliases):
            yield "args", node.attr
        elif _is_context(base, context_aliases):
            yield "context", node.attr


def scan_primitive_attributes(directory=PRIMITIVES_DIR):
    """AST scan of ``directory/*.py`` for context and action.args reads.

    Returns ``{"context": {attr: {module, ...}}, "args": {attr: {...}}}``.

    Counted as a read (attribute load, or augmented assignment target):
      * ``self.context.X`` and ``self.action.args.X``;
      * ``context.X`` / ``action.args.X`` where ``context`` / ``action`` is a
        parameter of the function (``__init__(self, action, context)`` reads
        ``context.pipeline_logger``);
      * ``v.X`` where ``v = self.action.args`` or ``v = self.context`` earlier
        in the same function;
      * reads of parameter ``p`` inside a module-level function ``f`` defined
        in ``directory``, when some primitive calls ``f(..., self.action.args,
        ...)`` with ``self.action.args`` in the position of ``p`` (this is
        how ``plotlabel(args)`` reads ``args.filter``).
    Plain stores (``self.action.args.X = ...``) are not reads.
    """
    found = {"context": {}, "args": {}}
    trees = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SyntaxWarning)  # '\*' in docstrings
        for path in sorted(pathlib.Path(directory).glob("*.py")):
            trees[path.stem] = ast.parse(path.read_text(), filename=str(path))

    def add(kind, attr, module):
        found[kind].setdefault(attr, set()).add(module)

    module_funcs = {}   # name -> (module, FunctionDef)
    for module, tree in trees.items():
        for node in tree.body:
            if isinstance(node, ast.FunctionDef):
                module_funcs[node.name] = (module, node)

    helper_params = set()   # (func name, param index) passed self.action.args
    for module, tree in trees.items():
        for func in ast.walk(tree):
            if not isinstance(func, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            params = [a.arg for a in func.args.args]
            ctx_al = {"context"} if "context" in params else set()
            act_al = {"action"} if "action" in params else set()
            for kind, attr in _reads_in_function(func, (), ctx_al, act_al):
                add(kind, attr, module)
        for call in ast.walk(tree):
            if isinstance(call, ast.Call) and isinstance(call.func, ast.Name) \
                    and call.func.id in module_funcs:
                for i, arg in enumerate(call.args):
                    if _chain(arg) == ["self", "action", "args"]:
                        helper_params.add((call.func.id, i))
    for fname, index in sorted(helper_params):
        module, func = module_funcs[fname]
        params = [a.arg for a in func.args.args]
        if index < len(params):
            for kind, attr in _reads_in_function(func, {params[index]}, (),
                                                 ()):
                if kind == "args":
                    add("args", attr, "%s.%s" % (module, fname))
    return found


def ingest_out_args_written(directory=PRIMITIVES_DIR):
    """Names assigned as ``out_args.X = ...`` in ingest_file._perform."""
    path = pathlib.Path(directory) / "kcwi_file_primitives.py"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SyntaxWarning)
        tree = ast.parse(path.read_text())
    names = []
    for cls in tree.body:
        if isinstance(cls, ast.ClassDef) and cls.name == "ingest_file":
            for func in cls.body:
                if isinstance(func, ast.FunctionDef) and \
                        func.name == "_perform":
                    for node in ast.walk(func):
                        targets = []
                        if isinstance(node, ast.Assign):
                            targets = node.targets
                        for t in targets:
                            elts = t.elts if isinstance(t, ast.Tuple) else [t]
                            for e in elts:
                                if _chain(e)[:1] == ["out_args"] and \
                                        len(_chain(e)) == 2:
                                    names.append(_chain(e)[1])
    return names
