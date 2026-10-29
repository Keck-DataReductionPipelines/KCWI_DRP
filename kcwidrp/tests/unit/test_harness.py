"""Tests for the primitive test harness (task T.2, kcwidrp/tests/harness.py).

Validation list of T.2 in docs/drp_upgrade_plan.md:
  * SubtractOverscan, TrimOverscan and CorrectGain run on synthetic data;
  * with plot_level = 0 a run calls no input(), time.sleep() or plotting
    function;
  * answers with one too many / one too few fail with a clear message;
  * a missing context/args attribute fails naming it and the primitive;
  * the AST scan of context/args reads is covered by the harness's lists.
"""
import os
from unittest import mock

import pytest
from keckdrpframework.models.action import Action
from keckdrpframework.models.arguments import Arguments
from keckdrpframework.models.processing_context import ProcessingContext
from keckdrpframework.primitives.base_primitive import BasePrimitive

from kcwidrp.primitives import SubtractOverscan as subtract_overscan_module
from kcwidrp.primitives.CorrectGain import CorrectGain
from kcwidrp.primitives.SubtractOverscan import SubtractOverscan
from kcwidrp.primitives.TrimOverscan import TrimOverscan
from kcwidrp.tests import harness

# SubtractOverscan fits the overscan only when its width exceeds minoscanpix
# (75 for BLUE in kcwi.cfg) and uses oscanbuf (20) columns of margin on each
# side.  The fixture's default overscan is 8 px, so these tests ask for a
# 100 px overscan on a 256 x 280 frame (4 amps of 128 rows x 140 columns:
# 40 data + 100 overscan columns).
WIDE_OSCAN = dict(overscan=100, shape=(256, 280))


# -- small test primitives ---------------------------------------------------

class _AsksTwice(BasePrimitive):
    """Calls input() twice and returns the answers."""

    def _perform(self):
        first = input("first? ")
        second = input("second? ")
        self.action.args.got = (first, second)
        return self.action.args


class _ReadsBarAvg(BasePrimitive):
    """Reads an args attribute that ingest does not set (FindBars sets it)."""

    def _perform(self):
        return self.action.args.bar_avg


class _ReadsArcs(BasePrimitive):
    """Reads a context attribute set by ExtractArcs in a reduction."""

    def _perform(self):
        return self.context.arcs


class _GetattrDefault(BasePrimitive):
    """Reads an optional attribute the way SolveAIT reads dich_frac."""

    def _perform(self):
        try:
            return self.action.args.dich_frac
        except AttributeError:
            return "fallback"


class _Raises(BasePrimitive):
    def _perform(self):
        raise ValueError("boom from primitive")


# -- the three primitives of the T.2 validation list -------------------------

def test_subtract_overscan_runs(run_primitive, synthetic_ccddata):
    ccd = synthetic_ccddata(**WIDE_OSCAN)
    out = run_primitive(SubtractOverscan, ccd)
    assert isinstance(out, Arguments)
    # source: analytic - overscan subtraction does not change the frame size.
    assert out.ccddata.shape == (256, 280)
    # source: analytic - SubtractOverscan sets OSCANSUB True when it fits.
    assert out.ccddata.header["OSCANSUB"] is True


def test_trim_overscan_runs(run_primitive, synthetic_ccddata):
    out = run_primitive(TrimOverscan, synthetic_ccddata())
    # source: analytic - default fixture frame 136 x 316, 4 amps of 32 data
    # columns x 64 data rows: trimmed width 2 * 32 = 64, height 2 * 64 = 128.
    assert out.ccddata.shape == (128, 64)
    assert out.ccddata.header["OSCANTRM"] is True     # source: analytic
    assert "BSEC1" not in out.ccddata.header          # source: analytic


def test_correct_gain_runs_on_trimmed_output(run_primitive,
                                            synthetic_ccddata):
    trimmed = run_primitive(TrimOverscan, synthetic_ccddata())
    # Chained: the TrimOverscan output is the CorrectGain action's args.
    out = run_primitive(CorrectGain, trimmed)
    assert out.ccddata.shape == (128, 64)             # source: analytic
    assert out.ccddata.header["GAINCOR"] is True      # source: analytic
    assert out.ccddata.header["BUNIT"] == "electron"  # source: analytic


# -- the objects the harness builds -------------------------------------------

def test_action_and_context_are_real_framework_objects(run_primitive,
                                                       synthetic_ccddata,
                                                       tmp_redux):
    run_primitive(TrimOverscan, synthetic_ccddata())
    action, context = run_primitive.last.action, run_primitive.last.context
    assert type(action) is Action
    assert isinstance(action.args, Arguments)
    assert isinstance(context, ProcessingContext)
    assert action.name == "TrimOverscan"
    assert context.bokeh_session is None              # D11 / T.11
    assert context.proctab is tmp_redux.proctab
    inst = context.config.instrument
    assert inst.plot_level == 0                       # source: T.2 default
    assert inst.plot_pause == 0                       # source: T.2 default
    assert inst.minoscanpix == 75   # source: kcwi.cfg [BLUE] minoscanpix
    assert inst.oscanbuf == 20      # source: kcwi.cfg [BLUE] oscanbuf
    assert inst.procfile == str(tmp_redux.proctab_file)
    # The framework config is framework.cfg, as Framework.__init__ reads it.
    assert context.config.name == "KCWI-Pipeline"     # source: framework.cfg
    # Files are written relative to the tmp_redux root, not the repository.
    assert os.getcwd() == str(tmp_redux.root)


def test_red_channel_values_and_overrides(run_primitive, synthetic_ccddata):
    # map_ccd is supplied because the RED synthetic frame numbers its amp
    # sections 1-4, while ingest_file.map_ccd reads RED sections 0-3.
    run_primitive(_GetattrDefault, synthetic_ccddata(camera="RED"),
                  config_overrides={"oscanbuf": 3, "plot_pause": 7},
                  action_extra={"map_ccd": None})
    inst = run_primitive.last.context.config.instrument
    assert inst.minoscanpix == 20   # source: kcwi.cfg [RED] minoscanpix
    assert inst.oscanbuf == 3       # source: override above
    assert inst.plot_pause == 7     # source: override above


def test_ingest_args_match_ingest_file_perform():
    # Every out_args attribute ingest_file._perform writes is provided, and
    # the harness provides nothing it does not write.
    assert set(harness.ingest_out_args_written()) == harness.ARGS_PROVIDED


def test_unbuildable_ingest_arg_names_it(run_primitive, synthetic_ccddata):
    ccd = synthetic_ccddata()
    del ccd.header["IFUNUM"]
    with pytest.raises(harness.HarnessIngestError,
                       match=r"could not build action\.args\.ifunum .*"
                             r"action_extra=\{'ifunum'"):
        run_primitive(TrimOverscan, ccd)


def test_exception_from_primitive_propagates(run_primitive,
                                            synthetic_ccddata):
    with pytest.raises(ValueError, match="boom from primitive"):
        run_primitive(_Raises, synthetic_ccddata())


# -- plot_level = 0: no input(), sleep() or plotting --------------------------

def _spy_plotting(monkeypatch):
    spies = {name: mock.MagicMock(name=name) for name in
             ("bokeh_plot", "save_plot", "figure", "bokeh_clear", "gridplot")}
    for name, spy in spies.items():
        monkeypatch.setattr(subtract_overscan_module, name, spy)
    return spies


def test_plot_level_0_calls_no_input_sleep_or_plot(run_primitive,
                                                   synthetic_ccddata,
                                                   monkeypatch):
    spies = _spy_plotting(monkeypatch)
    run_primitive(SubtractOverscan, synthetic_ccddata(**WIDE_OSCAN))
    assert run_primitive.last.prompts == []
    assert run_primitive.last.sleeps == []
    for name, spy in spies.items():
        assert spy.call_count == 0, name


def test_plot_spies_see_calls_at_plot_level_2(run_primitive,
                                              synthetic_ccddata, monkeypatch):
    # Shows the spies above are on the names SubtractOverscan calls, so the
    # plot_level = 0 test is not vacuous.
    spies = _spy_plotting(monkeypatch)
    run_primitive(SubtractOverscan, synthetic_ccddata(**WIDE_OSCAN),
                  config_overrides={"plot_level": 2})
    assert spies["figure"].call_count == 4      # source: one plot per amp
    assert spies["bokeh_plot"].call_count == 4  # source: one plot per amp
    # source: plot_pause defaults to 0 in the harness; one pause per plot.
    assert run_primitive.last.sleeps == [0, 0, 0, 0]
    assert run_primitive.last.prompts == []


# -- answers ------------------------------------------------------------------

def test_answers_fed_in_order(run_primitive, synthetic_ccddata):
    out = run_primitive(_AsksTwice, synthetic_ccddata(), answers=["a", "b"])
    assert out.got == ("a", "b")
    assert run_primitive.last.prompts == ["first? ", "second? "]


def test_one_answer_too_many_fails(run_primitive, synthetic_ccddata):
    with pytest.raises(pytest.fail.Exception,
                       match=r"too many answers: _AsksTwice returned with 1 "
                             r"of 3 answer\(s\) unused: \['c'\]"):
        run_primitive(_AsksTwice, synthetic_ccddata(),
                      answers=["a", "b", "c"])


def test_one_answer_too_few_fails(run_primitive, synthetic_ccddata):
    with pytest.raises(pytest.fail.Exception,
                       match=r"too few answers: _AsksTwice called "
                             r"input\('second\? '\) after all 1 answer"):
        run_primitive(_AsksTwice, synthetic_ccddata(), answers=["a"])


def test_input_without_answers_fails_naming_prompt(run_primitive,
                                                   synthetic_ccddata,
                                                   monkeypatch):
    _spy_plotting(monkeypatch)
    with pytest.raises(pytest.fail.Exception,
                       match=r"SubtractOverscan called input\('Next\? <cr>: "
                             r"'\), but the test gave no answers"):
        run_primitive(SubtractOverscan, synthetic_ccddata(**WIDE_OSCAN),
                      config_overrides={"plot_level": 3})


def test_unused_answers_fail_even_with_no_prompt(run_primitive,
                                                 synthetic_ccddata):
    with pytest.raises(pytest.fail.Exception, match="too many answers"):
        run_primitive(TrimOverscan, synthetic_ccddata(), answers=[""])


# -- missing attributes ---------------------------------------------------------

def test_missing_args_attribute_names_it(run_primitive, synthetic_ccddata):
    with pytest.raises(harness.HarnessAttributeError,
                       match=r"_ReadsBarAvg read action\.args\.bar_avg, which "
                             r"the harness does not provide \(in a reduction "
                             r"it is set by FindBars\); pass it with "
                             r"action_extra="):
        run_primitive(_ReadsBarAvg, synthetic_ccddata())


def test_missing_context_attribute_names_it(run_primitive,
                                            synthetic_ccddata):
    with pytest.raises(harness.HarnessAttributeError,
                       match=r"_ReadsArcs read context\.arcs, .*ExtractArcs"
                             r".*context_extra="):
        run_primitive(_ReadsArcs, synthetic_ccddata())


def test_extras_supply_missing_attributes(run_primitive, synthetic_ccddata):
    assert run_primitive(_ReadsBarAvg, synthetic_ccddata(),
                         action_extra={"bar_avg": 5}) == 5
    assert run_primitive(_ReadsArcs, synthetic_ccddata(),
                         context_extra={"arcs": [1]}) == [1]


def test_missing_attribute_is_still_an_attribute_error(run_primitive,
                                                       synthetic_ccddata):
    # try/except AttributeError in a primitive behaves as in a reduction.
    assert run_primitive(_GetattrDefault, synthetic_ccddata()) == "fallback"


# -- attribute scan -------------------------------------------------------------

def test_scan_is_covered_by_harness_lists():
    found = harness.scan_primitive_attributes()
    ctx, args = set(found["context"]), set(found["args"])
    # Report (shown with pytest -rP or -s).
    print("context reads (%d): %s" % (len(ctx), sorted(ctx)))
    print("action.args reads (%d): %s" % (len(args), sorted(args)))
    ctx_known = harness.CONTEXT_PROVIDED | set(harness.CONTEXT_NOT_PROVIDED)
    args_known = harness.ARGS_PROVIDED | set(harness.ARGS_NOT_PROVIDED)
    assert ctx - ctx_known == set(), \
        "context attributes read by primitives but not documented in " \
        "harness (add to CONTEXT_NOT_PROVIDED with what sets them)"
    assert args - args_known == set(), \
        "action.args attributes read by primitives but not documented in " \
        "harness (add to ARGS_NOT_PROVIDED with what sets them)"
    # No stale "not provided" entries.
    assert set(harness.CONTEXT_NOT_PROVIDED) - ctx == set()
    assert set(harness.ARGS_NOT_PROVIDED) - args == set()


def test_scan_finds_each_read_pattern():
    found = harness.scan_primitive_attributes()
    # self.action.args.X
    assert "SubtractOverscan" in found["args"]["map_ccd"]
    # context.X in __init__(self, action, context)
    assert "SubtractOverscan" in found["context"]["pipeline_logger"]
    # v = self.action.args; v.X
    assert "MakeMasterDark" in found["args"]["new_type"]
    # helper called with self.action.args: plotlabel(args) reads args.filter
    assert "kcwi_file_primitives.plotlabel" in found["args"]["filter"]
    # self.context.X
    assert "ExtractArcs" in found["context"]["trace"]


def test_scan_reports_a_new_attribute(tmp_path):
    (tmp_path / "NewPrim.py").write_text(
        "class NewPrim:\n"
        "    def __init__(self, action, context):\n"
        "        self.logger = context.pipeline_logger\n"
        "    def _perform(self):\n"
        "        self.action.args.ccddata = 1\n"        # a store, not a read
        "        x = self.action.args.brand_new_arg\n"
        "        return self.context.brand_new_ctx\n")
    found = harness.scan_primitive_attributes(tmp_path)
    assert set(found["args"]) == {"brand_new_arg"}
    assert set(found["context"]) == {"pipeline_logger", "brand_new_ctx"}
    assert "brand_new_arg" not in harness.ARGS_PROVIDED | set(
        harness.ARGS_NOT_PROVIDED)
