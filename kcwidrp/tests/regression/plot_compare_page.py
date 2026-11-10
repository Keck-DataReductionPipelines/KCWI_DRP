"""Side-by-side plot comparison page (bokeh plan 4.1, via upgrade amendment 8).

    python -m kcwidrp.tests.regression.plot_compare_page REF_DIR NEW_DIR OUT.html

Writes one HTML page with the reference and the new PNG of every plot side by
side, paired by relative file name, with each image's pixel size and a flag
where the sizes differ.  Files present on only one side are listed.  The page
refers to the images by path relative to ``OUT.html`` (or by ``file://`` URI
when no relative path exists, or with ``--file-uri``); it embeds nothing, so
keep it next to the directories it compares.

The pairing and the size rule are the ones of regression tier 6
(``compare.compare_outputs``): same names, same PNG pixel sizes.  The helpers
``compare._files_under`` and ``compare.png_size`` are reused.  Pass the
``plots/`` directories of a run and of its reference.
"""
import argparse
import html
import os
import sys
from pathlib import Path
from urllib.parse import quote

from .compare import _files_under, png_size

CSS = """
body { font-family: sans-serif; margin: 1em 2em; }
table { border-collapse: collapse; }
td, th { border: 1px solid #bbb; padding: 4px 8px; vertical-align: top; }
img { max-width: 48vw; height: auto; border: 1px solid #ddd; }
tr.bad td.status { background: #fcc; font-weight: bold; }
tr.ok td.status { background: #cfc; }
td.size { font-family: monospace; white-space: nowrap; }
"""


def _pngs(root):
    return {name: p for name, p in _files_under(root).items()
            if name.lower().endswith(".png")}


def _ref(path, out_dir, file_uri=False):
    """Relative URL of ``path`` from ``out_dir``, or a file:// URI."""
    path = Path(path).resolve()
    if not file_uri:
        try:
            rel = os.path.relpath(path, Path(out_dir).resolve())
            return quote(Path(rel).as_posix(), safe="/.-_~")
        except ValueError:          # another drive on Windows
            pass
    return path.as_uri()


def _fmt(size):
    return "unreadable" if size is None else "%d x %d" % size


def compare_dirs(ref_dir, new_dir):
    """Pair the PNGs of two directories by relative name."""
    ref, new = _pngs(ref_dir), _pngs(new_dir)
    pairs = []
    for name in sorted(set(ref) & set(new)):
        sr, sn = png_size(ref[name]), png_size(new[name])
        pairs.append({"name": name, "ref": ref[name], "new": new[name],
                      "ref_size": sr, "new_size": sn,
                      "size_match": sr == sn and sr is not None,
                      "identical": ref[name].read_bytes()
                      == new[name].read_bytes()})
    return {"pairs": pairs,
            "only_ref": sorted(set(ref) - set(new)),
            "only_new": sorted(set(new) - set(ref))}


def render_page(result, ref_dir, new_dir, out_dir, file_uri=False):
    e = html.escape
    pairs = result["pairs"]
    n_bad = sum(not p["size_match"] for p in pairs)
    n_same = sum(p["identical"] for p in pairs)
    out = ["<!DOCTYPE html>", "<html><head><meta charset='utf-8'>",
           "<title>Plot comparison</title><style>%s</style></head><body>" % CSS,
           "<h1>Plot comparison</h1>",
           "<p>Reference: <code>%s</code><br>New: <code>%s</code></p>" % (
               e(str(ref_dir)), e(str(new_dir))),
           "<p>%d paired, %d with a size mismatch, %d byte-identical; "
           "%d only in the reference, %d only in the new run.</p>" % (
               len(pairs), n_bad, n_same, len(result["only_ref"]),
               len(result["only_new"]))]
    for title, key in (("Only in the reference (missing from the new run)",
                        "only_ref"), ("Only in the new run", "only_new")):
        out.append("<h2>%s</h2>" % title)
        if result[key]:
            out.append("<ul>%s</ul>" % "".join(
                "<li>%s</li>" % e(n) for n in result[key]))
        else:
            out.append("<p>None.</p>")
    out.append("<h2>Paired plots</h2><table><tr><th>file</th>"
               "<th>reference</th><th>new</th><th>status</th></tr>")
    for p in pairs:
        if not p["size_match"]:
            status = "SIZE MISMATCH: reference %s, new %s" % (
                _fmt(p["ref_size"]), _fmt(p["new_size"]))
        elif p["identical"]:
            status = "same size, identical bytes"
        else:
            status = "same size, pixels differ"
        out.append(
            "<tr class='%s'><td>%s</td>"
            "<td><img src='%s' alt='reference'><br>"
            "<span class='size'>%s</span></td>"
            "<td><img src='%s' alt='new'><br>"
            "<span class='size'>%s</span></td>"
            "<td class='status'>%s</td></tr>" % (
                "ok" if p["size_match"] else "bad", e(p["name"]),
                e(_ref(p["ref"], out_dir, file_uri)), e(_fmt(p["ref_size"])),
                e(_ref(p["new"], out_dir, file_uri)), e(_fmt(p["new_size"])),
                e(status)))
    out.append("</table></body></html>")
    return "\n".join(out)


def write_page(ref_dir, new_dir, out_html, file_uri=False):
    """Write the page; returns the pairing result (see ``compare_dirs``)."""
    out_html = Path(out_html)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    result = compare_dirs(ref_dir, new_dir)
    out_html.write_text(render_page(result, ref_dir, new_dir, out_html.parent,
                                    file_uri), encoding="utf-8")
    return result


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python -m kcwidrp.tests.regression.plot_compare_page",
        description="Write an HTML page comparing two directories of PNGs.")
    ap.add_argument("ref_dir")
    ap.add_argument("new_dir")
    ap.add_argument("out_html")
    ap.add_argument("--file-uri", action="store_true",
                    help="refer to the images by file:// URI, not by "
                         "relative path")
    args = ap.parse_args(argv)
    for d in (args.ref_dir, args.new_dir):
        if not Path(d).is_dir():
            ap.error("not a directory: %s" % d)
    r = write_page(args.ref_dir, args.new_dir, args.out_html, args.file_uri)
    print("%s: %d paired (%d size mismatches), %d only in the reference, "
          "%d only in the new run" % (
              args.out_html, len(r["pairs"]),
              sum(not p["size_match"] for p in r["pairs"]),
              len(r["only_ref"]), len(r["only_new"])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
