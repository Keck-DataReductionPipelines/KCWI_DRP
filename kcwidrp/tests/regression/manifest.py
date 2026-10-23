"""Regression manifest: load, validate, checksums and the draft builder (T.7).

The data directory (``$KCWI_REGRESSION_DATA``) holds::

    manifest.json
    raw/<dataset>/        lev0/*.fits, files_flatlamp_<channel>.txt, kcwi.sky
    reference/<entry>/    redux/, plots/, <procfile>, prompt_record.jsonl
    answers/<file>        level-3 answer files

``raw/<dataset>`` and ``reference/<entry>`` may be symlinks.
"""
import hashlib
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

SCHEMA_VERSION = 1
CHANNELS = {"blue": ("-b", "kcwib.proc"), "red": ("-r", "kcwir.proc")}
PROFILES = ("quick", "full")
LEVELS = (0, 1, 3)
# S.3 conventions (docs/baseline/README.md): L1 is the quick profile.
LEVEL_PROFILE = {0: "full", 1: "quick", 3: "full"}
DEFAULT_OVERRIDES = {"enable_bokeh": False, "plot_pause": 0}
DEFAULT_ENVIRONMENT_FILE = "docs/baseline/env_baseline.txt"


class ManifestError(Exception):
    """The manifest or the data directory is unusable."""


def data_dir_from_env():
    value = os.environ.get("KCWI_REGRESSION_DATA")
    if not value:
        raise ManifestError("KCWI_REGRESSION_DATA not set")
    return Path(value)


def sha256_file(path, block=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(block)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def manifest_path(data_dir):
    return Path(data_dir) / "manifest.json"


def load_manifest(data_dir):
    path = manifest_path(data_dir)
    try:
        with open(path) as f:
            manifest = json.load(f)
    except FileNotFoundError:
        raise ManifestError("no manifest at %s" % path)
    except json.JSONDecodeError as exc:
        raise ManifestError("%s is not valid JSON: %s" % (path, exc))
    errors = validate_manifest(manifest)
    if errors:
        raise ManifestError("invalid manifest %s:\n  " % path
                            + "\n  ".join(errors))
    return manifest


def save_manifest(data_dir, manifest):
    """Write the manifest atomically."""
    path = manifest_path(data_dir)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=False)
        f.write("\n")
    os.replace(tmp, path)


def validate_manifest(manifest):
    """Return a list of problems (empty when the manifest is well formed)."""
    errors = []
    if not isinstance(manifest, dict):
        return ["top level is not an object"]
    if manifest.get("schema_version") != SCHEMA_VERSION:
        errors.append("schema_version must be %d" % SCHEMA_VERSION)
    datasets = manifest.get("datasets")
    entries = manifest.get("entries")
    if not isinstance(datasets, dict):
        errors.append("'datasets' must be an object keyed by dataset id")
        datasets = {}
    if not isinstance(entries, list):
        errors.append("'entries' must be a list")
        entries = []
    for ds_id, ds in datasets.items():
        if not isinstance(ds, dict) or "path" not in ds or "files" not in ds:
            errors.append("dataset %s needs 'path' and 'files'" % ds_id)
    ids = set()
    for ent in entries:
        eid = ent.get("id", "<no id>") if isinstance(ent, dict) else "<?>"
        pre = "entry %s: " % eid
        if not isinstance(ent, dict):
            errors.append(pre + "not an object")
            continue
        if eid in ids:
            errors.append(pre + "duplicate id")
        ids.add(eid)
        for key in ("id", "dataset", "channel", "file_list", "plot_level",
                    "profile", "config_overrides", "compare", "tolerances",
                    "provenance", "bless_history"):
            if key not in ent:
                errors.append(pre + "missing '%s'" % key)
        if ent.get("dataset") not in datasets:
            errors.append(pre + "unknown dataset %r" % ent.get("dataset"))
        if ent.get("channel") not in CHANNELS:
            errors.append(pre + "channel must be blue or red")
        if ent.get("profile") not in PROFILES:
            errors.append(pre + "profile must be quick or full")
        if ent.get("plot_level") not in (0, 1, 2, 3):
            errors.append(pre + "plot_level must be 0..3")
        elif ent.get("plot_level") == 3 and not ent.get("answers"):
            errors.append(pre + "plot_level 3 needs an 'answers' file")
        ov = ent.get("config_overrides", {})
        if not isinstance(ov, dict):
            errors.append(pre + "config_overrides must be an object")
        elif "plot_level" in ov and ov["plot_level"] != ent.get("plot_level"):
            errors.append(pre + "config_overrides.plot_level disagrees with "
                          "plot_level")
        cmp_ = ent.get("compare", {})
        if cmp_.get("mode") not in ("self", "against"):
            errors.append(pre + "compare.mode must be self or against")
        for pat, tol in (ent.get("tolerances") or {}).items():
            if not isinstance(tol, dict) or not tol.get("reason"):
                errors.append(pre + "tolerance %r needs a written 'reason' "
                              "(docs/testing.md section 6)" % pat)
    for ent in entries:
        cmp_ = ent.get("compare", {}) if isinstance(ent, dict) else {}
        if cmp_.get("mode") == "against" and cmp_.get("entry") not in ids:
            errors.append("entry %s: compare.entry %r is not an entry" % (
                ent.get("id"), cmp_.get("entry")))
    return errors


def get_entry(manifest, entry_id):
    for ent in manifest["entries"]:
        if ent["id"] == entry_id:
            return ent
    raise ManifestError("no entry %r in the manifest" % entry_id)


def select_entries(manifest, profile="quick"):
    """Entries of the quick profile, or all of them for ``full`` (D16)."""
    if profile == "full":
        return list(manifest["entries"])
    return [e for e in manifest["entries"] if e["profile"] == "quick"]


def dataset_dir(data_dir, manifest, dataset_id):
    """Absolute (not symlink-resolved) directory of a raw dataset."""
    return Path(os.path.abspath(
        Path(data_dir) / manifest["datasets"][dataset_id]["path"]))


def reference_dir(data_dir, manifest, entry):
    """Directory of the reference this entry is compared with (may not exist)."""
    cmp_ = entry["compare"]
    ref_id = cmp_["entry"] if cmp_["mode"] == "against" else entry["id"]
    return Path(data_dir) / "reference" / ref_id


def read_file_list(path):
    """Non-blank, non-comment lines of a file list."""
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                out.append(line)
    return out


def entry_raw_files(data_dir, manifest, entry):
    """Relative paths (under the dataset dir) of the FITS an entry reduces."""
    ds_dir = dataset_dir(data_dir, manifest, entry["dataset"])
    return read_file_list(ds_dir / entry["file_list"])


def verify_checksums(data_dir, manifest, entry):
    """Tier 0. Return a list of problems with the entry's raw files."""
    ds = manifest["datasets"][entry["dataset"]]
    ds_dir = dataset_dir(data_dir, manifest, entry["dataset"])
    problems = []
    try:
        listed = entry_raw_files(data_dir, manifest, entry)
    except OSError as exc:
        return ["cannot read file list %s: %s" % (entry["file_list"], exc)]
    for rel in listed:
        want = ds["files"].get(rel)
        path = ds_dir / rel
        if want is None:
            problems.append("%s: no checksum in the manifest" % rel)
        elif not path.is_file():
            problems.append("%s: missing" % rel)
        else:
            got = sha256_file(path)
            if got != want:
                problems.append("%s: SHA-256 %s... != manifest %s..." % (
                    rel, got[:12], want[:12]))
    ans = entry.get("answers")
    if ans and entry.get("answers_sha256"):
        path = Path(data_dir) / "answers" / ans
        if not path.is_file():
            problems.append("answers/%s: missing" % ans)
        elif sha256_file(path) != entry["answers_sha256"]:
            problems.append("answers/%s: SHA-256 differs from manifest" % ans)
    elif ans and not (Path(data_dir) / "answers" / ans).is_file():
        problems.append("answers/%s: missing" % ans)
    return problems


# -- draft ------------------------------------------------------------------

def short_dataset_name(dataset_id):
    """``ds01_F22491_pt2_60266`` -> ``ds01``."""
    return dataset_id.split("_")[0]


def read_reference_provenance(ref_dir):
    """DRP version and git version from the HISTORY of a reference cube."""
    prov = {"drp_version": None, "drp_git_version": None, "date": None}
    ref_dir = Path(ref_dir)
    redux = ref_dir / "redux"
    if not redux.is_dir():
        return prov
    try:
        from astropy.io import fits
    except ImportError:                                  # pragma: no cover
        return prov
    cubes = sorted(redux.glob("*_icubes.fits")) or sorted(
        redux.glob("*_icube*.fits"))
    for cube in cubes[:1]:
        with fits.open(cube) as hdul:
            for line in hdul[0].header.get("HISTORY", []):
                m = re.search(r"kcwidrp version=(\S+)", str(line))
                if m:
                    prov["drp_version"] = m.group(1)
                m = re.search(r"git version=(\S+)", str(line))
                if m:
                    prov["drp_git_version"] = m.group(1)
    procs = sorted(ref_dir.glob("*.proc"))
    if procs:
        mtime = os.stat(procs[0]).st_mtime
        prov["date"] = datetime.fromtimestamp(mtime, timezone.utc).strftime(
            "%Y-%m-%d")
    return prov


def draft_manifest(data_dir):
    """Build a draft manifest from the raw datasets found under ``raw/``.

    Returns ``(manifest, notes)``; ``notes`` lists what could not be determined.
    """
    data_dir = Path(data_dir)
    notes = []
    raw = data_dir / "raw"
    if not raw.is_dir():
        raise ManifestError("no raw/ directory in %s" % data_dir)
    datasets, entries = {}, []
    for ds_path in sorted(p for p in raw.iterdir() if p.is_dir()):
        ds_id = ds_path.name
        files = {}
        for channel in CHANNELS:
            lst = ds_path / ("files_flatlamp_%s.txt" % channel)
            if not lst.is_file():
                continue
            rels = read_file_list(lst)
            for rel in rels:
                if rel not in files:
                    f = ds_path / rel
                    if f.is_file():
                        files[rel] = sha256_file(f)
                    else:
                        notes.append("%s: listed file %s is missing; no "
                                     "checksum" % (ds_id, rel))
            short = short_dataset_name(ds_id)
            for level in LEVELS:
                eid = "%s_%s_L%d" % (short, channel, level)
                ans = None
                if level == 3:
                    ans = "answers_%s_%s_3.txt" % (ds_id, channel)
                ref = data_dir / "reference" / eid
                prov = read_reference_provenance(ref)
                if not ref.exists():
                    notes.append("%s: no reference/%s (bless it with reason "
                                 "'initial' after a first run)" % (eid, eid))
                elif prov["drp_version"] is None:
                    notes.append("%s: DRP version not found in the "
                                 "reference cube HISTORY" % eid)
                elif prov["drp_git_version"] is None:
                    notes.append("%s: no 'git version=' in the reference "
                                 "HISTORY (not an editable install); only "
                                 "version %s is recorded" % (
                                     eid, prov["drp_version"]))
                prov["environment_file"] = DEFAULT_ENVIRONMENT_FILE
                prov["iers"] = None
                ent = {
                    "id": eid,
                    "dataset": ds_id,
                    "channel": channel,
                    "description": "%s %s channel, plot_level %d" % (
                        ds_id, channel, level),
                    "file_list": lst.name,
                    "sky_file": "kcwi.sky" if (ds_path / "kcwi.sky").is_file()
                    else None,
                    "plot_level": level,
                    "profile": LEVEL_PROFILE[level],
                    "answers": ans,
                    "config_overrides": dict(DEFAULT_OVERRIDES,
                                             plot_level=level),
                    "compare": {"mode": "self"},
                    "tolerances": {},
                    "provenance": prov,
                    "bless_history": [],
                }
                if ans:
                    apath = data_dir / "answers" / ans
                    if apath.is_file():
                        ent["answers_sha256"] = sha256_file(apath)
                    else:
                        notes.append("%s: answers/%s not found" % (eid, ans))
                entries.append(ent)
        if not files:
            notes.append("%s: no files_flatlamp_<channel>.txt; dataset "
                         "skipped" % ds_id)
            continue
        datasets[ds_id] = {
            "id": ds_id,
            "path": "raw/%s" % ds_id,
            "provenance": {"source": "KOA, via the maintainer's test suite"},
            "files": files,
        }
    if entries:
        notes.append("IERS table versions (D18) not determined for any "
                     "entry; provenance.iers is null until a bless records it")
    manifest = {"schema_version": SCHEMA_VERSION, "datasets": datasets,
                "entries": entries}
    return manifest, notes
