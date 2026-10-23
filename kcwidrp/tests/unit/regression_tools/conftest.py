"""Fixtures for the regression-tooling tests (T.7): a synthetic world.

``make_world(root)`` builds a tiny data directory (one raw dataset of two
dummy "FITS" files, a file list, a sky file, a 3-line answers file, a drafted
manifest). ``FAKE_SCRIPT`` is a stand-in for ``reduce_kcwi``: it writes a few
small outputs (two FITS with several HDUs, a text file, a procfile, two PNGs,
DRPF/KCWI logs) and calls ``input()`` ``--inputs N`` times. It exits 0 even
with ``--fail``, which writes the framework's failure lines to the log, like
the real framework (plan section 6.1).
"""
import shutil
import sys
from types import SimpleNamespace

import pytest

from kcwidrp.tests.regression import manifest as mf
from kcwidrp.tests.regression import runner

FAKE_SCRIPT = r'''
import os, struct, sys, time, zlib
import numpy as np
from astropy.io import fits

args = sys.argv[1:]
n_inputs = int(args[args.index("--inputs") + 1]) if "--inputs" in args else 0
fail = "--fail" in args
proc = "kcwir.proc" if "-r" in args else "kcwib.proc"
for d in ("redux", "plots", "logs"):
    os.makedirs(d, exist_ok=True)
rng = np.random.default_rng(7)
now = time.strftime("%Y-%m-%dT%H:%M:%S")


def png(path, w, h):
    def chunk(tag, body):
        return (struct.pack(">I", len(body)) + tag + body
                + struct.pack(">I", zlib.crc32(tag + body)))
    raw = b"".join(b"\x00" + b"\x00" * w for _ in range(h))
    with open(path, "wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n"
                + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 0, 0, 0, 0))
                + chunk(b"IDAT", zlib.compress(raw)) + chunk(b"IEND", b""))


prim = fits.PrimaryHDU()
prim.header["CWAVE"] = 4500.5
prim.header["OBJECT"] = "star"
# A fixed observation time, copied from the raw header as the DRP does.
prim.header["DATE-OBS"] = "2024-04-06T05:39:35.591"
prim.header["INFILE"] = os.getcwd() + "/files.txt"
prim.header.add_history("run at " + now)
flux = fits.ImageHDU(rng.normal(size=(3, 4, 5)).astype("f4"), name="FLUX")
mask = fits.ImageHDU(np.zeros((3, 4, 5), "u1"), name="MASK")
tab = fits.BinTableHDU.from_columns(
    [fits.Column(name="WAVE", format="D", array=np.linspace(4000, 5000, 6)),
     fits.Column(name="NAME", format="8A", array=list("abcdef"))], name="TAB")
fits.HDUList([prim, flux, mask, tab]).writeto("redux/ds_icubes.fits",
                                              overwrite=True)
p2 = fits.PrimaryHDU(rng.normal(size=(6, 6)))
p2.header["EXPTIME"] = 12.0
p2.writeto("redux/ds_int.fits", overwrite=True)
open("redux/ds_atlines.txt", "w").write("line list\n")
open(proc, "w").write(
    "| FRAMENO | TYPE | filename |\n"
    "|       1 | BIAS | %s/lev0/KB.1.fits |\n"
    "|       2 | ARCLAMP | %s/lev0/KB.2.fits |\n" % (os.getcwd(), os.getcwd()))
png("plots/p1.png", 8, 6)
png("plots/p2.png", 9, 7)
with open("logs/DRPF.log", "w") as f:
    for ts, act in (("10:00:00", "start"), ("10:00:05", "make_cube"),
                    ("10:00:12", "make_cube"), ("10:00:14", "other")):
        f.write("2026-01-01 %s:DRPF:INFO: Event completed: name n, action %s,"
                " arg name x, recurr False\n" % (ts, act))
    if fail:
        f.write("2026-01-01 10:00:15:DRPF:ERROR: Exception in a primitive "
                "while invoking MakeCube\n")
with open("logs/KCWI.log", "w") as f:
    f.write("2026-01-01 10:00:01:KCWI:INFO: started\n")
    f.write("2026-01-01 10:00:02:KCWI:ERROR:CorrectDefects.py: Defect list not"
            " found for /x/data/defect_L2U2L1U1_2x2.dat\n")
    f.write("2026-01-01 10:00:03:KCWI:ERROR:MakeCube.py: No reference "
            "geometry, cannot make cube!\n")
    if fail:
        f.write("2026-01-01 10:00:15:KCWI:ERROR:Failed executing primitive "
                "MakeCube: boom\n")
for i in range(n_inputs):
    input("prompt %d: " % i)
'''

DATASET = "ds9_fake"


def write_fake_script(root):
    path = root / "fake_reduce.py"
    path.write_text(FAKE_SCRIPT)
    return path


def make_command(script, n_inputs_by_level=None, fail_ids=(), calls=None):
    """A reduction-command factory like ``runner.default_command``."""
    n_inputs_by_level = n_inputs_by_level or {3: 3}

    def command(entry, run_dir):
        if calls is not None:
            calls.append(entry["id"])
        argv = [sys.executable, str(script),
                mf.CHANNELS[entry["channel"]][0], "--inputs",
                str(n_inputs_by_level.get(entry["plot_level"], 0))]
        if entry["id"] in fail_ids:
            argv.append("--fail")
        return argv
    return command


def build_raw(root):
    """data dir with raw dataset, answers and a drafted manifest (no references)."""
    data = root / "data"
    ds = data / "raw" / DATASET
    (ds / "lev0").mkdir(parents=True)
    (ds / "lev0" / "KB.1.fits").write_bytes(b"raw frame one " * 50)
    (ds / "lev0" / "KB.2.fits").write_bytes(b"raw frame two " * 50)
    (ds / "files_flatlamp_blue.txt").write_text(
        "lev0/KB.1.fits\nlev0/KB.2.fits\n")
    (ds / "kcwi.sky").write_text("# sky\nKB.1.fits KB.2.fits\n")
    (data / "answers").mkdir()
    (data / "answers" / ("answers_%s_blue_3.txt" % DATASET)).write_text(
        "a\nb\nc\n")
    manifest, _notes = mf.draft_manifest(data)
    mf.save_manifest(data, manifest)
    return data


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """Raw data, manifest, fake command and blessed-looking references."""
    root = tmp_path_factory.mktemp("regworld")
    data = build_raw(root)
    script = write_fake_script(root)
    mfst = mf.load_manifest(data)
    cmd = make_command(script)
    for eid in ("ds9_blue_L1", "ds9_blue_L3"):
        entry = mf.get_entry(mfst, eid)
        res = runner.run_entry(data, mfst, entry, root / "refrun" / eid,
                               command=cmd)
        assert res.status == "ok"
        ref = data / "reference" / eid
        ref.mkdir(parents=True)
        run = root / "refrun" / eid
        for name in ("redux", "plots"):
            shutil.copytree(run / name, ref / name)
        shutil.copy(run / "kcwib.proc", ref)
        if (run / "prompt_record.jsonl").is_file():
            shutil.copy(run / "prompt_record.jsonl", ref)
    return SimpleNamespace(root=root, data=data, mfst=mfst, script=script,
                           command=cmd)


@pytest.fixture
def make_world(tmp_path):
    """Factory for a private world without references."""
    def _make():
        data = build_raw(tmp_path)
        script = write_fake_script(tmp_path)
        return SimpleNamespace(root=tmp_path, data=data,
                               mfst=mf.load_manifest(data), script=script,
                               command=make_command(script))
    return _make


@pytest.fixture(scope="module")
def baseline_run(world):
    """A second run of ds9_blue_L1, to be compared with its reference."""
    entry = mf.get_entry(world.mfst, "ds9_blue_L1")
    run_dir = world.root / "baseline" / "ds9_blue_L1"
    res = runner.run_entry(world.data, world.mfst, entry, run_dir,
                           command=world.command)
    assert res.status == "ok"
    return run_dir


@pytest.fixture
def command_factory():
    return make_command
