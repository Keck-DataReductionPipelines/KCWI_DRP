"""The regression comparison tiers (T.7).

Tier 0 checksums are in ``manifest.verify_checksums``; tier 1 completion,
tier 5 science metrics and tier 7 timings are here with tiers 2-4 and 6.

Rules (docs/testing.md sections 6 and 10):

* header cards whose keyword starts with HISTORY, COMMENT, CHECKSUM or
  DATASUM are ignored;
* a string value that is a file path (``_is_path``: absolute, or ending in
  a file extension the DRP writes) only has to agree in its file name, so
  the same input read from another directory is not a change; every other
  string (e.g. ``BUNIT = '1e+16 erg / (Angstrom s cm2)'``, ``OBJECT =
  'DOME/Dark/32.5 sec.'``) is compared exactly.
  DATE* cards are compared: in KCWI DRP outputs they are all observation
  times copied from the raw header (DATE, DATE-OBS, DATE-BEG, DATEREND, ...);
  the DRP writes no write-time date (checked on the S.3 outputs);
* ``*.pkl`` outputs (the geometry) are compared by content, with the same
  path rule (a stored file path only has to agree in its file name); other
  non-FITS files byte for byte;
* integers, booleans and strings are compared exactly; floats exactly unless
  the manifest gives a tolerance (with a reason) for the product;
* the ``.proc`` table is compared on every column; the ``filename`` column by
  base name;
* timings never fail.
"""
import fnmatch
import json
import math
import os
import pickle
import re
import struct
from pathlib import Path

import numpy as np

IGNORE_CARD_PREFIXES = ("HISTORY", "COMMENT", "CHECKSUM", "DATASUM")
MAX_FINDINGS = 200

TIER_NAMES = {
    0: "checksums", 1: "completion", 2: "structure", 3: "pixels",
    4: "header_values", 5: "metrics", 6: "plots", 7: "timings",
}

# Benign ERROR lines (docs/baseline/README.md), each allowed only for the
# channel that produces it: the red channel has no defect list for L2U2L1U1
# at 2x2 (B35), and its arcs before the last of a group reach MakeCube before
# the master arc exists (B34). Elsewhere these lines are failures. Tier 1 also
# checks that every object frame and master arc has its cube, so an allowed
# line cannot hide a missing cube.
ALLOWED_ERRORS = {
    "red": [
        r"CorrectDefects\.py: Defect list not found for "
        r".*defect_L2U2L1U1_2x2\.dat",
        r"MakeCube\.py: No reference geometry, cannot make cube!",
    ],
}
PROCFILES = {"blue": "kcwib.proc", "red": "kcwir.proc"}
FAILURE_PATTERNS = ("Failed executing primitive", "while invoking")


class TierResult:
    """Outcome of one tier for one entry."""

    def __init__(self, tier, status, summary="", details=None, data=None):
        self.tier = tier
        self.name = TIER_NAMES[tier]
        self.status = status          # pass | fail | info | skipped
        self.summary = summary
        self.details = list(details or [])
        self.data = data or {}

    def to_dict(self):
        return {"tier": self.tier, "name": self.name, "status": self.status,
                "summary": self.summary, "details": self.details[:50],
                "n_details": len(self.details), "data": self.data}


def apply_informational(result, entry):
    """For ``compare.mode == against``: tiers listed there are reported only."""
    cmp_ = entry.get("compare", {})
    if (cmp_.get("mode") == "against" and result.status == "fail"
            and result.tier in cmp_.get("informational_tiers", [])):
        result.status = "info"
        result.summary = "informational (not failed): " + result.summary
    return result


# -- tolerances -----------------------------------------------------------------

def tolerance_for(tolerances, filename, extname=None):
    """(rtol, atol) for a product; exact (0, 0) when no pattern matches.

    Patterns are globs matched against the file's base name or against
    ``<base name>::<extension name>``.
    """
    base = os.path.basename(filename)
    keys = [base] + (["%s::%s" % (base, extname)] if extname else [])
    for pattern, tol in (tolerances or {}).items():
        if any(fnmatch.fnmatchcase(k, pattern) for k in keys):
            return float(tol.get("rtol", 0.0)), float(tol.get("atol", 0.0))
    return 0.0, 0.0


# -- array comparison -------------------------------------------------------------

def diff_arrays(a, b, rtol=0.0, atol=0.0):
    """Compare two same-shape arrays; return None when they agree.

    Otherwise a dict with the number of differing values, the fraction, the
    largest absolute and relative differences and the location of the largest
    absolute difference.
    """
    a, b = np.asarray(a), np.asarray(b)
    if a.dtype.kind not in "iufcb" or b.dtype.kind not in "iufcb":
        # strings and objects: exact only, no numeric statistics
        bad = a != b
        if not bad.any():
            return None
        first = tuple(int(i) for i in np.argwhere(bad)[0])
        return {"n_diff": int(bad.sum()), "fraction": float(bad.mean()),
                "max_abs": float("nan"), "max_rel": float("nan"),
                "location": first, "n_nan_mismatch": 0}
    inexact = a.dtype.kind in "fc" or b.dtype.kind in "fc"
    if not inexact:
        if np.array_equal(a, b):
            return None
        bad = a != b
        af, bf = _as_float(a), _as_float(b)
    else:
        if rtol == 0.0 and atol == 0.0:
            if np.array_equal(a, b, equal_nan=True):
                return None
            with np.errstate(invalid="ignore"):
                bad = ~((a == b) | (np.isnan(a) & np.isnan(b)))
        else:
            with np.errstate(invalid="ignore"):
                close = np.abs(a - b) <= atol + rtol * np.abs(b)
                bad = ~(close | (a == b) | (np.isnan(a) & np.isnan(b)))
            if not bad.any():
                return None
        af, bf = a.astype(np.float64), b.astype(np.float64)
    n_bad = int(np.count_nonzero(bad))
    with np.errstate(invalid="ignore", divide="ignore"):
        absdiff = np.where(bad, np.abs(af - bf), -1.0)
        absdiff = np.where(np.isnan(absdiff), -1.0, absdiff)
        flat = int(np.argmax(absdiff))
        loc = tuple(int(i) for i in np.unravel_index(flat, absdiff.shape))
        max_abs = float(absdiff.flat[flat])
        rel_mask = bad & (bf != 0) & np.isfinite(af - bf)
        max_rel = (float(np.max(np.abs(af - bf)[rel_mask]
                               / np.abs(bf)[rel_mask]))
                   if rel_mask.any() else float("nan"))
    n_nan = int(np.count_nonzero(np.isnan(af) != np.isnan(bf))) if inexact \
        else 0
    return {"n_diff": n_bad, "fraction": n_bad / max(a.size, 1),
            "max_abs": max_abs if max_abs >= 0 else float("nan"),
            "max_rel": max_rel, "location": loc, "n_nan_mismatch": n_nan}


def _as_float(x):
    return x.astype(np.float64) if x.dtype.kind in "iub" else x


def _fmt_stats(s):
    return ("%d of the values differ (fraction %.3g), max |diff| %.6g, "
            "max rel diff %.3g, at index %s" % (
                s["n_diff"], s["fraction"], s["max_abs"], s["max_rel"],
                s["location"])
            + (", %d NaN mismatches" % s["n_nan_mismatch"]
               if s["n_nan_mismatch"] else ""))


def _is_int(v):
    return isinstance(v, int) and not isinstance(v, bool)


def numbers_agree(va, vb, rtol, atol):
    """Header values: two integers (or bools, strings, ...) exactly; the
    tolerance applies only when at least one value is a float."""
    if _is_int(va) and _is_int(vb):
        return va == vb
    if isinstance(va, float) and isinstance(vb, float):
        if math.isnan(va) and math.isnan(vb):
            return True
        if va == vb:
            return True
        return abs(va - vb) <= atol + rtol * abs(vb)
    if (isinstance(va, (int, float)) and isinstance(vb, (int, float))
            and not isinstance(va, bool) and not isinstance(vb, bool)):
        return (va == vb) or abs(va - vb) <= atol + rtol * abs(vb)
    return va == vb


# -- FITS and tables ----------------------------------------------------------------

def _header_keys(header):
    return {k for k in header
            if k and not k.startswith(IGNORE_CARD_PREFIXES)}


def _hdu_label(i, hdu):
    return hdu.name if hdu.name else "HDU%d" % i


def compare_fits_pair(run_path, ref_path, relname, tolerances, out,
                      diag_dir=None):
    """Append structure, pixel and header findings to ``out``."""
    from astropy.io import fits
    with fits.open(run_path, memmap=True) as ha, \
            fits.open(ref_path, memmap=True) as hb:
        if len(ha) != len(hb):
            out.structure.append("%s: %d HDUs, reference has %d" % (
                relname, len(ha), len(hb)))
        for i, (a, b) in enumerate(zip(ha, hb)):
            label = _hdu_label(i, b)
            where = "%s[%s]" % (relname, label)
            if a.name != b.name:
                out.structure.append("%s: HDU %d is named %r, reference %r"
                                     % (relname, i, a.name, b.name))
            rtol, atol = tolerance_for(tolerances, relname, b.name)
            ka, kb = _header_keys(a.header), _header_keys(b.header)
            if ka != kb:
                out.structure.append(
                    "%s: header keywords differ; only in run %s, only in "
                    "reference %s" % (where, sorted(ka - kb)[:8],
                                      sorted(kb - ka)[:8]))
            for k in sorted(ka & kb):
                va, vb = a.header[k], b.header[k]
                if _is_path(va) and _is_path(vb):
                    if not _same_file_name(va, vb):
                        out.headers.append("%s: %s = %r, reference %r" % (
                            where, k, va, vb))
                    continue          # same file read from another directory
                if not numbers_agree(va, vb, rtol, atol):
                    out.headers.append("%s: %s = %r, reference %r" % (
                        where, k, va, vb))
            _compare_hdu_data(a, b, relname, label, where, rtol, atol, out,
                              diag_dir)


def _compare_hdu_data(a, b, relname, label, where, rtol, atol, out,
                      diag_dir):
    da, db = a.data, b.data
    if (da is None) != (db is None):
        out.structure.append("%s: data present in one only" % where)
        return
    if da is None:
        return
    if hasattr(da, "columns") != hasattr(db, "columns"):
        out.structure.append("%s: image in one, table in the other" % where)
        return
    if hasattr(da, "columns"):
        na, nb = list(da.columns.names), list(db.columns.names)
        if na != nb or len(da) != len(db):
            out.structure.append(
                "%s: table columns/rows differ (%s x %d vs %s x %d)" % (
                    where, na, len(da), nb, len(db)))
            return
        for name in na:
            ca, cb = np.asarray(da[name]), np.asarray(db[name])
            if ca.dtype != cb.dtype or ca.shape != cb.shape:
                out.structure.append("%s: column %s dtype/shape %s %s vs "
                                     "%s %s" % (where, name, ca.dtype,
                                                ca.shape, cb.dtype, cb.shape))
                continue
            _record_pixel(out, "%s column %s" % (where, name), relname,
                          label, ca, cb, rtol, atol, None)
        return
    if da.shape != db.shape or da.dtype != db.dtype:
        out.structure.append("%s: shape/dtype %s %s, reference %s %s" % (
            where, da.shape, da.dtype, db.shape, db.dtype))
        return
    _record_pixel(out, where, relname, label, da, db, rtol, atol, diag_dir)


def _record_pixel(out, where, relname, label, a, b, rtol, atol, diag_dir):
    stats = diff_arrays(a, b, rtol, atol)
    if stats is None:
        return
    finding = {"file": relname, "ext": label, "where": where, **stats,
               "rtol": rtol, "atol": atol}
    if diag_dir is not None and a.ndim >= 2 and a.dtype.kind in "iufb":
        try:
            from . import report
            loc = stats["location"]
            plane = loc[0] if a.ndim >= 3 else None
            sl = (loc[0],) if a.ndim >= 3 else ()
            ia = np.asarray(a[sl], dtype=np.float64)
            ib = np.asarray(b[sl], dtype=np.float64)
            while ia.ndim > 2:           # 4-D and up: flatten leading axes
                ia, ib = ia[0], ib[0]
            png = Path(diag_dir) / (re.sub(r"[^\w.-]+", "_", where) + ".png")
            report.write_diff_image(ia, ib, png, "%s%s" % (
                where, "" if plane is None else ", plane %d" % plane))
            finding["diff_png"] = str(png)
        except Exception as exc:        # a diagnostic must never hide a failure
            finding["diff_png_error"] = repr(exc)
    out.pixels.append(finding)
    out.pixel_messages.append("%s: %s" % (where, _fmt_stats(stats)))


def read_proc(path):
    from astropy.table import Table
    return Table.read(str(path), format="ascii.fixed_width")  # as the DRP does


def compare_proc_tables(run_path, ref_path, name, out):
    ta, tb = read_proc(run_path), read_proc(ref_path)
    if ta.colnames != tb.colnames:
        out.structure.append("%s: columns %s vs %s" % (
            name, ta.colnames, tb.colnames))
        return
    if len(ta) != len(tb):
        out.structure.append("%s: %d rows, reference %d" % (
            name, len(ta), len(tb)))
        return
    for col in ta.colnames:
        va, vb = [str(x) for x in ta[col]], [str(x) for x in tb[col]]
        if col == "filename":
            va = [os.path.basename(x) for x in va]
            vb = [os.path.basename(x) for x in vb]
        if va != vb:
            idx = next(i for i in range(len(va)) if va[i] != vb[i])
            n = sum(x != y for x, y in zip(va, vb))
            msg = "%s: column %s differs in %d row(s); first row %d: %r vs " \
                  "reference %r" % (name, col, n, idx, va[idx], vb[idx])
            out.pixels.append({"file": name, "ext": "table", "where": msg,
                               "n_diff": n})
            out.pixel_messages.append(msg)


class OutputComparison:
    """Tiers 2, 3, 4 and 6 for one run against one reference."""

    def __init__(self):
        self.structure = []
        self.pixels = []
        self.pixel_messages = []
        self.headers = []
        self.plots = []
        self.n_files = 0
        self.n_plots = 0


PATH_SUFFIXES = (".fits", ".fit", ".txt", ".dat", ".pkl", ".proc", ".cfg")


def _is_path(value):
    """True for a string that is a file path: absolute, or ending in one of
    the file extensions the DRP writes. A slash alone does not make a path
    (units such as ``erg / (Angstrom s cm2)``, object names such as
    ``DOME/Dark/32.5 sec.``)."""
    if not isinstance(value, str):
        return False
    v = value.strip()
    return (v.startswith("/") or v.lower().endswith(PATH_SUFFIXES)
            or re.match(r"^[A-Za-z]:\\", v) is not None)


def _same_file_name(va, vb):
    return os.path.basename(va.strip().replace("\\", "/")) == \
        os.path.basename(vb.strip().replace("\\", "/"))


def diff_objects(a, b, where="", limit=20):
    """Differences between two unpickled DRP objects, as strings.

    Dicts by key, sequences by element, arrays exactly (NaN equal to NaN),
    other objects (e.g. skimage transforms) by their ``__dict__``, numbers and
    other values with ``==``. A string that is a file path only has to agree
    in its file name, the same rule as for header cards: the same input read
    from another directory (e.g. through the regression data directory's
    symlinks) is not a change.
    """
    found = []

    def walk(x, y, path):
        if len(found) >= limit:
            return
        if type(x) is not type(y):
            found.append("%s: type %s, reference %s"
                         % (path, type(x).__name__, type(y).__name__))
        elif isinstance(x, dict):
            for k in sorted(set(x) | set(y), key=repr):
                if k not in x or k not in y:
                    found.append("%s[%r]: only in the %s" % (
                        path, k, "run" if k in x else "reference"))
                else:
                    walk(x[k], y[k], "%s[%r]" % (path, k))
        elif isinstance(x, (list, tuple)):
            if len(x) != len(y):
                found.append("%s: length %d, reference %d"
                             % (path, len(x), len(y)))
            else:
                for i, (u, v) in enumerate(zip(x, y)):
                    walk(u, v, "%s[%d]" % (path, i))
        elif isinstance(x, np.ndarray):
            if x.shape != y.shape or x.dtype != y.dtype:
                found.append("%s: shape/dtype %s %s, reference %s %s"
                             % (path, x.shape, x.dtype, y.shape, y.dtype))
            elif not np.array_equal(x, y, equal_nan=x.dtype.kind in "fc"):
                found.append("%s: %d array values differ"
                             % (path, int(np.count_nonzero(x != y))))
        elif _is_path(x) and _is_path(y):
            if not _same_file_name(x, y):
                found.append("%s: %r, reference %r" % (path, x, y))
        elif hasattr(x, "__dict__") and not isinstance(x, type):
            walk(vars(x), vars(y), path)
        else:
            same = x == y
            if isinstance(same, np.ndarray):
                same = bool(same.all())
            if not same and not (isinstance(x, float) and math.isnan(x)
                                 and math.isnan(y)):
                found.append("%s: %r, reference %r" % (path, x, y))

    walk(a, b, where)
    return found


def compare_pickle_pair(run_path, ref_path, rel, out):
    """Tier 3 for ``*.pkl`` outputs (the geometry): compare the contents.

    Only the run's own output and the maintainer's blessed reference are
    unpickled. A file that cannot be unpickled is compared byte for byte.
    """
    try:
        with open(run_path, "rb") as fa, open(ref_path, "rb") as fb:
            a, b = pickle.load(fa), pickle.load(fb)
    except Exception:
        if Path(run_path).read_bytes() != Path(ref_path).read_bytes():
            msg = "%s: bytes differ (could not be unpickled)" % rel
            out.pixels.append({"file": rel, "ext": "-", "where": msg,
                               "n_diff": 1})
            out.pixel_messages.append(msg)
        return
    for d in diff_objects(a, b):
        msg = "%s: %s" % (rel, d)
        out.pixels.append({"file": rel, "ext": "-", "where": msg,
                           "n_diff": 1})
        out.pixel_messages.append(msg)


def _files_under(root):
    root = Path(root)
    if not root.is_dir():
        return {}
    return {str(p.relative_to(root)): p for p in sorted(root.rglob("*"))
            if p.is_file() or (p.is_symlink() and p.exists())}


def png_size(path):
    """(width, height) from the PNG IHDR chunk, or None."""
    with open(path, "rb") as f:
        head = f.read(24)
    if len(head) < 24 or head[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return struct.unpack(">II", head[16:24])


def compare_outputs(run_dir, ref_dir, procfile, tolerances=None,
                    diag_dir=None):
    """Compare ``redux/``, the procfile and ``plots/`` of a run with a reference."""
    run_dir, ref_dir = Path(run_dir), Path(ref_dir)
    out = OutputComparison()
    tolerances = tolerances or {}
    if diag_dir is not None:
        Path(diag_dir).mkdir(parents=True, exist_ok=True)

    # file sets (tier 2): redux/ and the procfile
    run_files = {"redux/" + k: v for k, v in
                 _files_under(run_dir / "redux").items()}
    ref_files = {"redux/" + k: v for k, v in
                 _files_under(ref_dir / "redux").items()}
    run_proc, ref_proc = run_dir / procfile, ref_dir / procfile
    if run_proc.is_file():
        run_files[procfile] = run_proc
    if ref_proc.is_file():
        ref_files[procfile] = ref_proc
    only_run = sorted(set(run_files) - set(ref_files))
    only_ref = sorted(set(ref_files) - set(run_files))
    if only_run:
        out.structure.append("files only in the run: %s" % only_run)
    if only_ref:
        out.structure.append("files missing from the run: %s" % only_ref)

    for rel in sorted(set(run_files) & set(ref_files)):
        out.n_files += 1
        a, b = run_files[rel], ref_files[rel]
        try:
            if rel == procfile:
                compare_proc_tables(a, b, rel, out)
            elif rel.endswith((".fits", ".fit")):
                compare_fits_pair(a, b, rel[len("redux/"):], tolerances, out,
                                  diag_dir)
            elif rel.endswith(".pkl"):
                compare_pickle_pair(a, b, rel, out)
            elif a.read_bytes() != b.read_bytes():
                msg = "%s: bytes differ (not a FITS file)" % rel
                out.pixels.append({"file": rel, "ext": "-", "where": msg,
                                   "n_diff": 1})
                out.pixel_messages.append(msg)
        except Exception as exc:
            out.structure.append("%s: could not be compared: %r" % (rel, exc))

    # plots (tier 6)
    pr, pb = _files_under(run_dir / "plots"), _files_under(ref_dir / "plots")
    out.n_plots = len(pb)
    if set(pr) - set(pb):
        out.plots.append("plots only in the run: %s" % sorted(set(pr) - set(pb)))
    if set(pb) - set(pr):
        out.plots.append("plots missing from the run: %s" % sorted(
            set(pb) - set(pr)))
    for name in sorted(set(pr) & set(pb)):
        if name.lower().endswith(".png"):
            sa, sb = png_size(pr[name]), png_size(pb[name])
            if sa != sb:
                out.plots.append("%s: size %s, reference %s" % (name, sa, sb))
    return out


def structure_result(cmp_):
    n = len(cmp_.structure)
    return TierResult(2, "pass" if not n else "fail",
                      "%d files compared, same structure" % cmp_.n_files
                      if not n else "%d structural difference(s)" % n,
                      cmp_.structure[:MAX_FINDINGS])


def pixels_result(cmp_):
    n = len(cmp_.pixels)
    files = sorted({f["file"] for f in cmp_.pixels})
    return TierResult(
        3, "pass" if not n else "fail",
        "all values equal in %d files" % cmp_.n_files if not n else
        "%d extension(s) differ in %d file(s): %s" % (
            n, len(files), ", ".join(files[:5])),
        cmp_.pixel_messages[:MAX_FINDINGS], {"findings": cmp_.pixels[:50]})


def headers_result(cmp_):
    n = len(cmp_.headers)
    return TierResult(4, "pass" if not n else "fail",
                      "header values equal" if not n else
                      "%d header value difference(s)" % n,
                      cmp_.headers[:MAX_FINDINGS])


def plots_result(cmp_):
    n = len(cmp_.plots)
    return TierResult(6, "pass" if not n else "fail",
                      "%d plot files, same names and PNG sizes" % cmp_.n_plots
                      if not n else "%d plot difference(s)" % n,
                      cmp_.plots[:MAX_FINDINGS])


# -- tier 1: completion -----------------------------------------------------------

def read_prompt_record(path):
    path = Path(path)
    if not path.is_file():
        return []
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            out.append(json.loads(line))
    return out


def scan_logs(run_dir, channel=None):
    """(failure lines, other ERROR lines, number of log files).

    ERROR lines in ``ALLOWED_ERRORS[channel]`` are not reported.
    """
    failures, errors, nlogs = [], [], 0
    log_dir = Path(run_dir) / "logs"
    allowed = [re.compile(p) for p in ALLOWED_ERRORS.get(channel, [])]
    for log in sorted(log_dir.glob("*.log")) if log_dir.is_dir() else []:
        nlogs += 1
        with open(log, errors="replace") as f:
            for line in f:
                line = line.rstrip("\n")
                if any(p in line for p in FAILURE_PATTERNS):
                    failures.append("%s: %s" % (log.name, line[:300]))
                elif ":ERROR:" in line and not any(
                        r.search(line) for r in allowed):
                    errors.append("%s: %s" % (log.name, line[:300]))
    return failures, errors, nlogs


def missing_cubes(run_dir, channel, data=None):
    """Object frames and master arcs in the run's proctab without a cube.

    Every frame the run's own proctab lists as ``OBJECT`` (science, standard
    and sky frames) and every ``MARC`` must have ``redux/<frame>_icube.fits``.
    This holds with or without a reference, so it also protects an initial
    bless.
    """
    run_dir = Path(run_dir)
    proc = run_dir / PROCFILES.get(channel, "")
    if not channel or not proc.is_file():
        return ["no proctab %s: cannot check that every object frame has a "
                "cube" % proc.name]
    try:
        tab = read_proc(proc)
    except Exception as exc:
        return ["proctab %s could not be read: %r" % (proc.name, exc)]
    want = [os.path.basename(str(r["filename"]))
            for r in tab if str(r["TYPE"]).strip() in ("OBJECT", "MARC")]
    missing = []
    for name in want:
        stem = name[:-len(".fits")] if name.endswith(".fits") else name
        if not (run_dir / "redux" / (stem + "_icube.fits")).is_file():
            missing.append(name)
    if data is not None:
        data["cubes_expected"] = len(want)
    if not want:
        return ["the proctab lists no object frame or master arc"]
    return ["no cube for %s (%s_icube.fits missing)" % (
        n, n[:-len(".fits")] if n.endswith(".fits") else n) for n in missing]


def failing_primitive(line):
    m = re.search(r"Failed executing primitive\W*(\w+)", line) or re.search(
        r"while invoking\W*(\w+)", line)
    return m.group(1) if m else None


def tier1_result(run_dir, entry, answers_path=None, ref_dir=None,
                 run_result=None):
    run_dir = Path(run_dir)
    problems = []
    data = {}
    if run_result is not None:
        data["exit_code"] = run_result.exit_code     # recorded, never checked
        if run_result.status in ("timeout", "error"):
            problems.append("the reduction did not finish: %s" %
                            run_result.reason)
    if not (run_dir / "redux").is_dir():
        problems.append("no redux/ directory: the reduction produced nothing")
    failures, errors, nlogs = scan_logs(run_dir, entry.get("channel"))
    if not nlogs:
        problems.append("no logs/*.log: the reduction did not start")
    for line in failures:
        prim = failing_primitive(line)
        problems.append("primitive failure%s: %s" % (
            " (%s)" % prim if prim else "", line))
    data["failed_primitives"] = sorted(
        {p for p in map(failing_primitive, failures) if p})
    for line in errors:
        problems.append("unexpected ERROR line: %s" % line)
    problems.extend(missing_cubes(run_dir, entry.get("channel"), data))

    record = read_prompt_record(run_dir / "prompt_record.jsonl")
    data["n_prompts"] = len(record)
    level = entry["plot_level"]
    if level == 3:
        lines = Path(answers_path).read_text().splitlines()
        data["n_answers"] = len(lines)
        eofs = [i for i, r in enumerate(record) if r.get("eof")]
        if eofs:
            problems.append("EOF at prompt %d (answers file too short): %r"
                            % (eofs[0] + 1, record[eofs[0]]["prompt"][:100]))
        answered = len(record) - len(eofs)
        if answered != len(lines):
            problems.append("prompt count %d != %d lines in the answers file"
                            % (answered, len(lines)))
        for i, (r, ans) in enumerate(zip(record, lines)):
            if "answer" in r and r["answer"] != ans:
                problems.append("prompt %d got answer %r, answers file line "
                                "is %r" % (i + 1, r["answer"], ans))
                break
        ref_record = (read_prompt_record(Path(ref_dir) / "prompt_record.jsonl")
                      if ref_dir else [])
        if ref_record:
            pa = [r["prompt"] for r in record]
            pb = [r["prompt"] for r in ref_record]
            if pa != pb:
                idx = next((i for i, (x, y) in enumerate(zip(pa, pb))
                            if x != y), min(len(pa), len(pb)))
                problems.append(
                    "prompt strings differ from the reference record at "
                    "prompt %d (run %d prompts, reference %d)" % (
                        idx + 1, len(pa), len(pb)))
        else:
            data["prompt_strings"] = "not checked: no reference prompt record"
    elif level in (0, 1):
        if record:
            problems.append("%d unexpected prompt(s), first: %r" % (
                len(record), record[0].get("prompt", "")[:100]))
    summary = ("no failures, %d prompts" % len(record) if not problems
               else "%d problem(s); first: %s" % (len(problems), problems[0]))
    return TierResult(1, "pass" if not problems else "fail", summary,
                      problems[:MAX_FINDINGS], data)


# -- tier 5: science metrics --------------------------------------------------------

METRICS = {}


def register_metric(name):
    """Register ``func(run_dir, ref_dir, entry) -> (ok, message)`` (M.0)."""
    def deco(func):
        METRICS[name] = func
        return func
    return deco


def tier5_result(run_dir, ref_dir, entry):
    if not METRICS:
        return TierResult(5, "pass", "no metrics registered (added by M.0)")
    bad = []
    for name, func in sorted(METRICS.items()):
        ok, message = func(run_dir, ref_dir, entry)
        if not ok:
            bad.append("%s: %s" % (name, message))
    return TierResult(5, "pass" if not bad else "fail",
                      "%d metric(s)" % len(METRICS) if not bad else
                      "%d metric(s) missed" % len(bad), bad)


# -- tier 7: timings ----------------------------------------------------------------

EVENT_RE = re.compile(
    r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d):DRPF:INFO: Event (completed|failed): "
    r"name (\S+), action (\w+)")


def parse_timings(log_path):
    """Seconds per action: the gap since the previous event, summed."""
    from datetime import datetime
    per = {}
    prev = None
    log_path = Path(log_path)
    if not log_path.is_file():
        return per
    with open(log_path, errors="replace") as f:
        for line in f:
            m = EVENT_RE.match(line)
            if not m:
                continue
            t = datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S")
            if prev is not None:
                rec = per.setdefault(m.group(4), {"seconds": 0.0, "calls": 0})
                rec["seconds"] += (t - prev).total_seconds()
                rec["calls"] += 1
            else:
                per.setdefault(m.group(4), {"seconds": 0.0, "calls": 0})
                per[m.group(4)]["calls"] += 1
            prev = t
    return per


def tier7_result(run_dir, wall_time):
    per = parse_timings(Path(run_dir) / "logs" / "DRPF.log")
    top = sorted(per.items(), key=lambda kv: -kv[1]["seconds"])
    details = ["%s: %.0f s in %d call(s)" % (k, v["seconds"], v["calls"])
               for k, v in top if v["seconds"] > 0][:15]
    summary = "wall time %s" % ("%.0f s" % wall_time if wall_time is not None
                                else "n/a")
    return TierResult(7, "info", summary, details,
                      {"wall_time": wall_time, "primitives": per})
