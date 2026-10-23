"""Bless a run as the new reference (T.7; rules in docs/testing.md section 10).

A bless never happens automatically and only the maintainer runs it. It needs
a reason from the plan's T.7 table, a run whose outputs still exist, and (unless
``--yes``) a typed confirmation.
"""
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from . import manifest as mf
from . import runner
from .report import SessionResults

# The plan's T.7 bless table, plus "initial" for an entry with no reference yet.
ALLOWED_REASONS = (
    "initial", "environment", "T.11 viewer-off", "bokeh port 4.3",
    "prompt text PR3a", "M.2c header", "M.2e initial", "M.3", "M.4", "M.5",
    "matrix default", "release version",
)
BUG_REASON = re.compile(r"P\.5 B\d+[a-z]*$")


class BlessError(Exception):
    """The bless was refused."""


def parse_reason(reason):
    """Split ``"<allowed prefix>: <text>"``; refuse anything else."""
    if not reason or not reason.strip():
        raise BlessError("--reason is required")
    head, sep, text = reason.partition(":")
    head, text = head.strip(), text.strip()
    if head not in ALLOWED_REASONS and not BUG_REASON.match(head):
        raise BlessError(
            "reason %r is not allowed. Use '<reason>: <text>' with one of: %s"
            ", or 'P.5 <bug ID>' (docs/drp_upgrade_plan.md, T.7 bless table)"
            % (head, ", ".join(ALLOWED_REASONS)))
    if not sep or not text:
        raise BlessError("the reason needs a written explanation after the "
                         "colon: '%s: <text>'" % head)
    return head, text


def clone_path(src, dst):
    """Copy a file or tree; APFS clone (``cp -c``) when it works."""
    src, dst = Path(src), Path(dst)
    if sys.platform == "darwin":
        proc = subprocess.run(["cp", "-c", "-R", str(src), str(dst)],
                              capture_output=True)
        if proc.returncode == 0:
            return "clone"
        if dst.exists():                 # a partial clone: start again
            remove_path(dst)
    if src.is_dir():
        shutil.copytree(src, dst, symlinks=True)
    else:
        shutil.copy2(src, dst)
    return "copy"


def remove_path(path):
    path = Path(path)
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.exists():
        shutil.rmtree(path)


def locate_run(entry_id, run=None, out_root=None):
    """(session dir, run dir) for an entry."""
    if run:
        run = Path(run)
        if (run / "results.json").is_file():
            session = run
            run_dir = session / entry_id
        elif (run / "run.json").is_file() and (
                run.parent / "results.json").is_file():
            session, run_dir = run.parent, run
        else:
            raise BlessError("%s is neither a session directory nor a run "
                             "directory" % run)
    else:
        session = runner.latest_session_dir(out_root)
        if session is None:
            raise BlessError("no session found under %s; pass --run" %
                             (out_root or runner.default_out_root()))
        run_dir = session / entry_id
    return session, run_dir


def check_blessable(results, entry_id, run_dir, procfile):
    rec = results.data["entries"].get(entry_id)
    if rec is None or not rec.get("run"):
        raise BlessError("entry %s was not run in session %s" % (
            entry_id, results.session_dir))
    if rec["run"]["status"] != "ok":
        raise BlessError("the run of %s did not complete (status %s)" % (
            entry_id, rec["run"]["status"]))
    for tier in (0, 1):
        if results.tier_status(entry_id, tier) == "fail":
            raise BlessError("tier %d failed for %s; a failed reduction "
                             "cannot become a reference" % (tier, entry_id))
    redux = run_dir / "redux"
    if rec.get("outputs_deleted") or not redux.is_dir() or not any(
            redux.iterdir()):
        raise BlessError(
            "the outputs of %s were deleted (a passing run is torn down); "
            "rerun with --keep-runs" % entry_id)
    if not (run_dir / procfile).is_file():
        raise BlessError("%s has no %s" % (run_dir, procfile))
    return rec


def bless(data_dir, entry_id, reason, run=None, out_root=None,
          assume_yes=False, input_fn=input, out=None, today=None):
    """Replace ``reference/<entry>`` with a run's outputs and record why."""
    out = out or sys.stdout
    head, text = parse_reason(reason)
    mfst = mf.load_manifest(data_dir)
    entry = mf.get_entry(mfst, entry_id)
    procfile = runner.procfile_name(entry)
    session, run_dir = locate_run(entry_id, run, out_root)
    results = SessionResults.load(session)
    rec = check_blessable(results, entry_id, run_dir, procfile)

    data_dir = Path(data_dir)
    ref = data_dir / "reference" / entry_id
    stamp = today or datetime.now().strftime("%Y-%m-%d")
    sources = [("redux", run_dir / "redux"), ("plots", run_dir / "plots"),
               (procfile, run_dir / procfile),
               ("prompt_record.jsonl", run_dir / "prompt_record.jsonl")]
    sources = [(n, p) for n, p in sources if p.exists()]
    old_prov = entry.get("provenance")
    new_prov = rec["run"].get("provenance") or runner.current_provenance()
    new_prov = dict(new_prov)
    new_prov.setdefault("environment_file", (old_prov or {}).get(
        "environment_file"))
    print("Bless %s" % entry_id, file=out)
    print("  reason:       %s: %s" % (head, text), file=out)
    print("  from run:     %s" % run_dir, file=out)
    print("  replaces:     %s%s" % (ref, "" if ref.exists() else
                                    " (does not exist yet)"), file=out)
    print("  copies:       %s" % ", ".join(n for n, _ in sources), file=out)
    print("  old code:     %s %s" % ((old_prov or {}).get("drp_git_version"),
                                    (old_prov or {}).get("drp_version")),
          file=out)
    print("  new code:     %s %s" % (new_prov.get("drp_git_version"),
                                    new_prov.get("drp_version")), file=out)
    if entry["compare"]["mode"] == "against":
        print("  compare mode: 'against %s' becomes 'self'" %
              entry["compare"]["entry"], file=out)
    for k, v in results.diff_summary(entry_id).items():
        print("  tier %s: %s - %s" % (k, v["status"], v["summary"]), file=out)
    if not assume_yes:
        answer = input_fn("Type 'yes' to replace the reference: ")
        if answer.strip().lower() != "yes":
            raise BlessError("not confirmed; nothing changed")

    new = ref.with_name(".%s.new" % entry_id)
    old = ref.with_name(".%s.old" % entry_id)
    for leftover in (new, old):
        remove_path(leftover)
    new.mkdir(parents=True)
    how = None
    try:
        for name, src in sources:
            how = clone_path(src, new / name)
    except BaseException:
        remove_path(new)               # the old reference is untouched
        raise
    if ref.exists() or ref.is_symlink():
        ref.rename(old)
    new.rename(ref)
    remove_path(old)

    history = {
        "date": stamp, "reason": "%s: %s" % (head, text),
        "run": str(run_dir), "old_provenance": old_prov,
        "new_provenance": new_prov, "diff_summary":
            results.diff_summary(entry_id),
        "compare_before": dict(entry["compare"]), "copy_method": how,
    }
    entry.setdefault("bless_history", []).append(history)
    entry["provenance"] = new_prov
    if entry["compare"]["mode"] == "against":
        entry["compare"] = {"mode": "self"}
    mf.save_manifest(data_dir, mfst)
    print("Blessed %s. History entry recorded in manifest.json." % entry_id,
          file=out)
    return history
