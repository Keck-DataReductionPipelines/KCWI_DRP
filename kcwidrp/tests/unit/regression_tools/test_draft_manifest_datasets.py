"""draft-manifest skips datasets without usable files (written from the spec only).

A dataset is kept only if at least one file named in its
files_flatlamp_<channel>.txt exists; its entries (<short>_<channel>_L<level>,
levels 0, 1, 3) are added only if the dataset is kept. The draft is always
valid, and the CLI validates it before writing (invalid: print errors, exit 2,
write nothing).
"""
import copy
import subprocess
import sys
from pathlib import Path

import kcwidrp
from kcwidrp.tests.regression import __main__ as cli
from kcwidrp.tests.regression import manifest as mf

GOOD = "ds01_good"
BAD = "ds02_allmissing"
NOLIST = "ds03_nolist"
REPO = Path(kcwidrp.__file__).resolve().parent.parent


def make_data(root, with_nolist=False):
    data = root / "data"
    good = data / "raw" / GOOD
    (good / "lev0").mkdir(parents=True)
    (good / "lev0" / "KB.1.fits").write_bytes(b"raw one " * 20)
    (good / "files_flatlamp_blue.txt").write_text("lev0/KB.1.fits\n")
    bad = data / "raw" / BAD
    bad.mkdir(parents=True)
    (bad / "files_flatlamp_blue.txt").write_text(
        "lev0/KB.7.fits\nlev0/KB.8.fits\n")
    if with_nolist:
        (data / "raw" / NOLIST / "lev0").mkdir(parents=True)
        (data / "raw" / NOLIST / "lev0" / "KB.9.fits").write_bytes(b"x")
    (data / "answers").mkdir()
    for ds in (GOOD, BAD):
        (data / "answers" / ("answers_%s_blue_3.txt" % ds)).write_text("a\n")
    return data


def test_all_missing_dataset_contributes_nothing(tmp_path):
    data = make_data(tmp_path)
    manifest, notes = mf.draft_manifest(data)
    assert list(manifest["datasets"]) == [GOOD]  # source: spec, bad dataset dropped
    # source: spec, levels 0, 1, 3 for the one channel with a list
    assert sorted(e["id"] for e in manifest["entries"]) == [
        "ds01_blue_L0", "ds01_blue_L1", "ds01_blue_L3"]
    assert all(e["dataset"] == GOOD for e in manifest["entries"])
    assert mf.validate_manifest(manifest) == []  # source: spec, always valid
    assert any("every file listed" in n and "skipped" in n for n in notes)


def test_dataset_without_list_gets_a_note(tmp_path):
    data = make_data(tmp_path, with_nolist=True)
    manifest, notes = mf.draft_manifest(data)
    assert NOLIST not in manifest["datasets"]  # source: spec, no list, no dataset
    assert not [e for e in manifest["entries"] if e["id"].startswith("ds03_")]
    assert mf.validate_manifest(manifest) == []
    assert any("no files_flatlamp" in n for n in notes)


def test_cli_writes_valid_manifest_without_bad_dataset(tmp_path):
    data = make_data(tmp_path)
    proc = subprocess.run(
        [sys.executable, "-m", "kcwidrp.tests.regression",
         "--data-dir", str(data), "draft-manifest"],
        cwd=str(REPO), capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    loaded = mf.load_manifest(data)  # raises if invalid
    assert list(loaded["datasets"]) == [GOOD]
    ids = sorted(e["id"] for e in loaded["entries"])
    assert ids == ["ds01_blue_L0", "ds01_blue_L1", "ds01_blue_L3"]
    assert not [i for i in ids if i.startswith("ds02")]


def test_cli_refuses_to_write_invalid_draft(tmp_path, monkeypatch, capsys):
    data = make_data(tmp_path)
    good, _ = mf.draft_manifest(data)
    bad = copy.deepcopy(good)
    bad["entries"][0]["dataset"] = "no_such_dataset"
    assert mf.validate_manifest(bad) != []  # source: spec, unknown dataset
    monkeypatch.setattr(mf, "draft_manifest", lambda d: (bad, []))
    out = tmp_path / "out.json"
    rc = cli.main(["--data-dir", str(data), "draft-manifest",
                   "--output", str(out)])
    assert rc == 2  # source: spec, invalid draft exits 2
    assert not out.exists()
    assert not (data / "manifest.json").exists()
    cap = capsys.readouterr()  # spec: errors printed, stream unspecified
    assert "no_such_dataset" in cap.out + cap.err
