"""CalcPrelimDisp run through ``run_primitive`` (task T.4, numeric case).

Specification (class docstring): from the grating angle and the camera
(articulation) angle in the header,

    BLUE: alpha = grating_angle - 13 - adjustment (180 for BH, else 0)
    RED:  alpha = 155.892 - grating_angle
    beta  = camera_angle - alpha

and the preliminary dispersion in A per binned pixel follows from the
grating equation.  ``context.prelim_disp`` receives it.

Expected values, computed here without the primitive's closed form:
the first-order grating equation with an out-of-plane angle gamma
(e.g. Schroeder, "Astronomical Optics", 2nd ed., eq. 14.2.1),

    lambda(beta) = d * cos(gamma) * (sin(alpha) + sin(beta)),   d = 1/rho,

differentiated numerically (central difference) in beta, times the angle a
binned pixel subtends at the camera, PIX * ybin / FCAM (rad).  PIX = 0.015
mm, FCAM = 305 mm and GAMMA = 4 deg are the documented constants in
kcwi.cfg; the groove densities rho (lines per micrometre) are the per-grating
table documented in ``ingest_file.rho``.  A second test checks the angle
convention against the real instrument: for two real 2026-01-15 arc
headers, the code's alpha and beta put the central wavelength CWAVE where
the grating equation says it should be.
"""
import logging
import math
import re

import numpy as np
import pytest

from kcwidrp.primitives.CalcPrelimDisp import CalcPrelimDisp

HEAD = "c0fed5c"
# source: published @ c0fed5c - kcwi.cfg
PIX, FCAM, GAMMA = 0.015, 305.0, 4.0

# grating: (camera, rho [lines/um], central wavelength used to pick a
# plausible geometry [A])
# source: characterization @ c0fed5c - rho from the table documented in
# ingest_file.rho ("lines/mm divided by 1000")
GRATINGS = {
    "BH1": ("BLUE", 3.6, 3600.0), "BH2": ("BLUE", 3.2805, 4200.0),
    "BH3": ("BLUE", 2.8017, 5000.0), "BM": ("BLUE", 1.901, 4500.0),
    "BL": ("BLUE", 0.870, 4600.0),
    "RH1": ("RED", 2.360, 6000.0), "RH2": ("RED", 2.020, 7000.0),
    "RH3": ("RED", 1.735, 8000.0), "RH4": ("RED", 1.490, 9500.0),
    "RM1": ("RED", 1.300, 7200.0), "RM2": ("RED", 0.945, 9000.0),
    "RL": ("RED", 0.530, 7500.0),
}


def grating_lambda(alpha_deg, beta_deg, rho):
    """First-order grating equation, wavelength in Angstrom."""
    a, b, g = (math.radians(x) for x in (alpha_deg, beta_deg, GAMMA))
    return 1.0e4 / rho * math.cos(g) * (math.sin(a) + math.sin(b))


def expected_dispersion(alpha_deg, beta_deg, rho, ybin):
    # source: analytic @ c0fed5c - numerical d(lambda)/d(beta) of the grating
    # equation times the angle of one binned pixel at the camera
    h = 1.0e-4                                       # deg
    dl_dbeta = (grating_lambda(alpha_deg, beta_deg + h, rho) -
                grating_lambda(alpha_deg, beta_deg - h, rho)) / \
        math.radians(2 * h)                          # A per radian
    return dl_dbeta * PIX * ybin / FCAM


def geometry(grating):
    """(alpha, beta, grating angle, camera angle) for a plausible setting:
    near-Littrow for the grating's typical central wavelength, with the
    camera 3 deg short of 2 * alpha so that alpha != beta."""
    camera, rho, cwave = GRATINGS[grating]
    s = cwave * rho / 1.0e4 / math.cos(math.radians(GAMMA))
    alpha = math.degrees(math.asin(s / 2.0)) + 1.5
    beta = alpha - 3.0
    if camera == "BLUE":
        adj = 180.0 if grating.startswith("BH") else 0.0
        grangle = alpha + 13.0 + adj
    else:
        grangle = 155.892 - alpha
    return alpha, beta, grangle, alpha + beta


def frame_for(synthetic_ccddata, grating, ybin):
    camera, _, _ = GRATINGS[grating]
    _, _, grangle, camangle = geometry(grating)
    c = "B" if camera == "BLUE" else "R"
    hdr = {c + "GRATNAM": grating, c + "GRANGLE": grangle,
           c + "ARTANG": camangle, "IMTYPE": "MARC"}
    return synthetic_ccddata(camera=camera, binning=(ybin, ybin), **hdr)


@pytest.mark.parametrize("ybin", [1, 2])
@pytest.mark.parametrize("grating", sorted(GRATINGS))
def test_preliminary_dispersion_per_grating(run_primitive, synthetic_ccddata,
                                            grating, ybin):
    alpha, beta, _, _ = geometry(grating)
    rho = GRATINGS[grating][1]
    run_primitive(CalcPrelimDisp, frame_for(synthetic_ccddata, grating, ybin))
    got = run_primitive.last.context.prelim_disp
    expected = expected_dispersion(alpha, beta, rho, ybin)
    # Tolerance: central difference with h = 1e-4 deg has truncation error
    # ~h^2/6 ~ 5e-10 relative and rounding ~1e-16 / 1.7e-6 rad ~ 1e-10;
    # 1e-8 relative is far below any change of angle convention or constant
    # (a 0.01 deg change of beta alters the dispersion by >= 1e-5).
    assert got == pytest.approx(expected, rel=1e-8, abs=0)
    # the geometry is a real one: it reaches the grating's central wavelength
    # (source: analytic @ c0fed5c - grating equation; Littrow +-1.5 deg)
    assert grating_lambda(alpha, beta, rho) == pytest.approx(
        GRATINGS[grating][2], rel=0.01)


# Real arc headers (instrument state of 2026-01-15)
# source: published @ c0fed5c - kb260115_00011 / kr260115_00013 headers
# (BL/RM1)
REAL = [
    ("BLUE", dict(BGRATNAM="BL", BGRANGLE=28.7556705, BARTANG=23.1307564),
     0.870, 4599.8271484),
    ("RED", dict(RGRATNAM="RM1", RGRANGLE=123.0774384, RARTANG=56.022934),
     1.300, 7199.9848633),
]


@pytest.mark.parametrize("camera, hdr, rho, cwave", REAL,
                         ids=["kb260115_00011", "kr260115_00013"])
def test_angle_convention_matches_real_central_wavelength(
        run_primitive, synthetic_ccddata, caplog, camera, hdr, rho, cwave):
    caplog.set_level(logging.INFO, logger="KCWI")
    frame = synthetic_ccddata(camera=camera, binning=(2, 2), IMTYPE="MARC",
                              **hdr)
    run_primitive(CalcPrelimDisp, frame)
    msg = [r.getMessage() for r in caplog.records
           if "Initial alpha, beta" in r.getMessage()]
    alpha, beta = map(float, re.findall(r"(-?\d+\.\d+)", msg[0]))
    # source: analytic @ c0fed5c - grating equation with the code's angles
    # gives the central wavelength recorded by the instrument.  Tolerance 1 %:
    # the detector centre is not exactly the optical axis and the angle
    # encoders' zero points are approximate (BL: 4585 vs 4600 A, RM1: 7182 vs
    # 7200 A, i.e. <= 0.3 %); dropping the 13 deg offset moves BL by 4.5 % and
    # flipping the sign of beta by > 60 %.
    assert grating_lambda(alpha, beta, rho) == pytest.approx(cwave, rel=0.01)
    disp = run_primitive.last.context.prelim_disp
    assert disp == pytest.approx(expected_dispersion(alpha, beta, rho, 2),
                                 rel=1e-4)    # Tolerance: angles logged
    # to 0.001 deg -> <= 3e-5 relative in cos(beta)


def test_not_a_master_arc_is_skipped(run_primitive, synthetic_ccddata):
    # source: analytic @ c0fed5c - docstring of _pre_condition - only MARC
    # frames
    run_primitive(CalcPrelimDisp, synthetic_ccddata(IMTYPE="OBJECT"))
    assert "prelim_disp" not in vars(run_primitive.last.context)
