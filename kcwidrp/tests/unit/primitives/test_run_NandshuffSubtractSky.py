"""NandshuffSubtractSky prompt site run through ``run_primitive`` (T.4).

No real nod-and-shuffle data exist for this upgrade (Q4), so the frames are
synthetic.  Panel layout (class code, BLUE): with NSHFDN != NSHFUP + 1 the
observation was aborted and the panels are inverted, sky in rows SHUFROWS
... 2*SHUFROWS-1, object above it; the sky row start (>= 10) makes this the
"non-standard" configuration, where the user is asked for a sky scale
factor (docs/plot_inventory.md P37/P38):

    ``Enter scale factor for sky to match obj (float): ``

at plot_level >= 2 (1.0 is used below that).  The object panel is moved onto
the sky panel and ``obj - scale * sky`` is formed there; NASSCL records the
scale.  Negative or non-numeric answers are asked again.

Known bugs pinned as strict xfails (D8): B7 (0, nan and inf are accepted as
the scale) and B12 (the RED non-standard branch indexes 2-D arrays with
three indices and raises after the scale is entered).
"""
import numpy as np
import pytest

from kcwidrp.primitives.NandshuffSubtractSky import NandshuffSubtractSky

HEAD = "c0fed5c"
SHUF = 20
PROMPT = "Enter scale factor for sky to match obj (float): "
NS = {"nasmask": True, "numopen": 2}     # set by ingest from BNASNAM/NUMOPEN
SKY, OBJ = 100.0, 350.0


def blue_frame(synthetic_ccddata):
    """BLUE, inverted (aborted) N&S: sky rows 20-39, object rows 40-59."""
    f = synthetic_ccddata(SHUFROWS=SHUF, NSHFUP=0, NSHFDN=0)
    f.data[:] = 0.0
    f.data[SHUF:2 * SHUF, :] = SKY
    f.data[2 * SHUF:3 * SHUF, :] = OBJ
    return f


def run(run_primitive, frame, plot_level, answers=None):
    return run_primitive(NandshuffSubtractSky, frame, action_extra=NS,
                         config_overrides={"plot_level": plot_level},
                         answers=answers).ccddata


def sky_panel(out):
    # rows SHUF+1 .. 2*SHUF-1 hold obj - scale*sky (row SHUF is cleared with
    # the rows below it)
    return out.data[SHUF + 1:2 * SHUF, :]


def test_plot_level_1_uses_scale_one_without_prompt(run_primitive,
                                                    synthetic_ccddata):
    out = run(run_primitive, blue_frame(synthetic_ccddata), 1)
    assert run_primitive.last.prompts == []
    # source: characterization @ c0fed5c (inventory P37: pl < 2 -> 1.0)
    assert out.header["NASSCL"] == 1.0
    # source: analytic @ c0fed5c - object panel minus sky panel
    assert np.all(sky_panel(out) == OBJ - SKY)


def test_valid_scale_is_applied(run_primitive, synthetic_ccddata):
    out = run(run_primitive, blue_frame(synthetic_ccddata), 2, ["1.5"])
    assert run_primitive.last.prompts == [PROMPT]
    # source: characterization @ c0fed5c (inventory P37) and analytic
    assert out.header["NASSCL"] == 1.5
    assert np.all(sky_panel(out) == OBJ - 1.5 * SKY)


def test_negative_and_non_numeric_are_asked_again(run_primitive,
                                                  synthetic_ccddata):
    # source: characterization @ c0fed5c (inventory P37)
    out = run(run_primitive, blue_frame(synthetic_ccddata), 2,
              ["-1", "abc", "", "0.5"])
    assert run_primitive.last.prompts == [PROMPT] * 4
    assert out.header["NASSCL"] == 0.5


B7_ZERO = pytest.mark.xfail(strict=True, raises=pytest.fail.Exception,
                            reason=(
    "B7: the sky scale loop only rejects negative numbers "
    "(NandshuffSubtractSky.py:183-192), so 0 is accepted (no sky "
    "subtraction, no warning) and the next answer is left over"))
B7_NONFINITE = pytest.mark.xfail(strict=True, raises=ValueError, reason=(
    "B7: nan and inf pass the 'skyscl < 0.' check "
    "(NandshuffSubtractSky.py:183-192); the subtraction goes ahead and "
    "writing NASSCL = nan/inf to the header raises ValueError"))


@pytest.mark.parametrize("bad", [pytest.param("0", marks=B7_ZERO),
                                 pytest.param("nan", marks=B7_NONFINITE),
                                 pytest.param("inf", marks=B7_NONFINITE)])
def test_zero_nan_inf_scale_is_asked_again(run_primitive, synthetic_ccddata,
                                           bad):
    # source: analytic @ c0fed5c - plan section 13 B7 - a scale of 0, nan or
    # inf is not a valid sky scale; it must be asked again
    out = run(run_primitive, blue_frame(synthetic_ccddata), 2, [bad, "1.5"])
    assert run_primitive.last.prompts == [PROMPT, PROMPT]
    assert out.header["NASSCL"] == 1.5


@pytest.mark.xfail(strict=True, raises=IndexError, reason=(
    "B12: RED non-standard N&S branch indexes 2-D arrays with three indices "
    "(sky[0, objrow0, :], NandshuffSubtractSky.py:303, 305) and raises "
    "IndexError after the scale is entered"))
def test_red_non_standard_configuration(run_primitive, synthetic_ccddata):
    # RED, 1x1, inverted (NSHFDN != NSHFUP): object rows 20-39, sky 40-59
    f = synthetic_ccddata(camera="RED", binning=(1, 1), SHUFROWS=SHUF,
                          NSHFUP=0, NSHFDN=1)
    f.data[:] = 0.0
    f.data[SHUF:2 * SHUF, :] = OBJ
    f.data[2 * SHUF:3 * SHUF, :] = SKY
    out = run(run_primitive, f, 2, ["1.0"])
    # source: analytic @ c0fed5c - the sky panel is moved onto the object panel
    # and subtracted there
    assert out.header["NASSCL"] == 1.0
    assert np.all(out.data[SHUF:2 * SHUF, :] == OBJ - SKY)
