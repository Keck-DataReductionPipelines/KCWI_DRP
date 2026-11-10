"""plot_compare_page: side-by-side PNG page (bokeh plan 4.1, amendment 8)."""
import re
import shutil
import struct
import subprocess
import sys
import zlib
from html import unescape
from urllib.parse import unquote

import pytest

from kcwidrp.tests.regression import plot_compare_page as pcp


def write_png(path, w, h, shade=0):
    def chunk(tag, body):
        return (struct.pack(">I", len(body)) + tag + body
                + struct.pack(">I", zlib.crc32(tag + body)))
    raw = b"".join(b"\x00" + bytes([shade]) * w for _ in range(h))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\x89PNG\r\n\x1a\n"
                     + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 0, 0,
                                                  0, 0))
                     + chunk(b"IDAT", zlib.compress(raw)) + chunk(b"IEND", b""))


@pytest.fixture
def ref(tmp_path):
    d = tmp_path / "ref" / "plots"
    write_png(d / "a_nlines.png", 100, 60, 10)
    write_png(d / "b_resid.png", 120, 80, 20)
    write_png(d / "oscan_1.png", 50, 40, 30)
    return d


def page(tmp_path, ref_dir, new_dir, name="out/compare.html"):
    out = tmp_path / name
    result = pcp.write_page(ref_dir, new_dir, out)
    return result, out.read_text()


def test_directory_against_itself_all_paired_none_unpaired(tmp_path, ref):
    result, text = page(tmp_path, ref, ref)
    # source: analytic -- a directory against itself pairs every file
    assert [p["name"] for p in result["pairs"]] == [
        "a_nlines.png", "b_resid.png", "oscan_1.png"]
    assert result["only_ref"] == [] and result["only_new"] == []
    assert all(p["size_match"] and p["identical"] for p in result["pairs"])
    assert "SIZE MISMATCH" not in text
    for name in ("a_nlines.png", "b_resid.png", "oscan_1.png"):
        assert name in text
    assert "100 x 60" in text and "120 x 80" in text and "50 x 40" in text
    assert "3 paired, 0 with a size mismatch" in text
    assert "<li>" not in text                   # no unpaired lists


def test_removed_and_renamed_file_both_reported(tmp_path, ref):
    new = tmp_path / "new" / "plots"
    shutil.copytree(ref, new)
    (new / "a_nlines.png").unlink()                       # one removed
    (new / "b_resid.png").rename(new / "b_resid_new.png")   # one renamed
    result, text = page(tmp_path, ref, new)
    # source: analytic -- removed file is only in ref; renamed file is
    # only-in-ref under the old name and only-in-new under the new name.
    assert result["only_ref"] == ["a_nlines.png", "b_resid.png"]
    assert result["only_new"] == ["b_resid_new.png"]
    assert [p["name"] for p in result["pairs"]] == ["oscan_1.png"]
    assert "<li>a_nlines.png</li>" in text
    assert "<li>b_resid.png</li>" in text
    assert "<li>b_resid_new.png</li>" in text


def test_size_mismatch_flagged_with_both_sizes(tmp_path, ref):
    new = tmp_path / "new" / "plots"
    shutil.copytree(ref, new)
    write_png(new / "b_resid.png", 130, 80, 20)
    result, text = page(tmp_path, ref, new)
    bad = [p for p in result["pairs"] if not p["size_match"]]
    assert [p["name"] for p in bad] == ["b_resid.png"]
    assert bad[0]["ref_size"] == (120, 80) and bad[0]["new_size"] == (130, 80)
    assert "SIZE MISMATCH: reference 120 x 80, new 130 x 80" in text
    assert text.count("SIZE MISMATCH") == 1
    assert "1 with a size mismatch" in text


def test_same_size_different_pixels_is_not_a_mismatch(tmp_path, ref):
    new = tmp_path / "new" / "plots"
    shutil.copytree(ref, new)
    write_png(new / "oscan_1.png", 50, 40, 99)
    result, text = page(tmp_path, ref, new)
    p = [q for q in result["pairs"] if q["name"] == "oscan_1.png"][0]
    assert p["size_match"] and not p["identical"]
    assert "SIZE MISMATCH" not in text and "pixels differ" in text


def test_image_references_are_relative_and_resolve(tmp_path, ref):
    new = tmp_path / "new" / "plots"
    shutil.copytree(ref, new)
    out = tmp_path / "deep" / "er" / "page.html"
    pcp.write_page(ref, new, out)
    srcs = re.findall(r"<img src='([^']*)'", out.read_text())
    assert len(srcs) == 6
    for s in srcs:
        s = unescape(s)
        assert not s.startswith(("/", "file:"))
        target = (out.parent / unquote(s)).resolve()
        assert target.is_file()
        assert target.suffix == ".png"


def test_file_uri_mode(tmp_path, ref):
    out = tmp_path / "page.html"
    pcp.write_page(ref, ref, out, file_uri=True)
    srcs = re.findall(r"<img src='([^']*)'", out.read_text())
    assert srcs and all(unescape(s).startswith("file://") for s in srcs)


def test_names_with_spaces_are_escaped(tmp_path):
    d = tmp_path / "r"
    write_png(d / "my plot & more.png", 10, 10)
    out = tmp_path / "p.html"
    pcp.write_page(d, d, out)
    text = out.read_text()
    assert "my plot &amp; more.png" in text
    src = unescape(re.findall(r"<img src='([^']*)'", text)[0])
    assert (out.parent / unquote(src)).is_file()


def test_non_png_files_ignored_and_subdirs_paired(tmp_path):
    r, n = tmp_path / "r", tmp_path / "n"
    for d in (r, n):
        write_png(d / "sub" / "x.png", 8, 8)
        (d / "notes.txt").write_text("hi")
    (r / "only.txt").write_text("x")
    result = pcp.compare_dirs(r, n)
    assert [p["name"] for p in result["pairs"]] == ["sub/x.png"]
    assert result["only_ref"] == [] and result["only_new"] == []


def test_empty_directories(tmp_path):
    (tmp_path / "r").mkdir()
    (tmp_path / "n").mkdir()
    result, text = page(tmp_path, tmp_path / "r", tmp_path / "n")
    assert result == {"pairs": [], "only_ref": [], "only_new": []}
    assert "0 paired" in text


def test_command_line_entry_point(tmp_path, ref):
    out = tmp_path / "cli.html"
    proc = subprocess.run(
        [sys.executable, "-m", "kcwidrp.tests.regression.plot_compare_page",
         str(ref), str(ref), str(out)], capture_output=True, text=True,
        timeout=120)
    assert proc.returncode == 0, proc.stderr
    assert out.is_file()
    assert "3 paired" in proc.stdout


def test_command_line_rejects_missing_directory(tmp_path, ref):
    proc = subprocess.run(
        [sys.executable, "-m", "kcwidrp.tests.regression.plot_compare_page",
         str(ref), str(tmp_path / "nope"), str(tmp_path / "o.html")],
        capture_output=True, text=True, timeout=120)
    assert proc.returncode != 0
    assert not (tmp_path / "o.html").exists()
