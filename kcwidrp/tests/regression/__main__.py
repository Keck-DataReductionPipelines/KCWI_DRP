"""Command line of the regression tooling: ``python -m kcwidrp.tests.regression``.

Commands: draft-manifest, verify-checksums, run, bless, report.
"""
import argparse
import sys
from pathlib import Path

from . import bless as blessmod
from . import manifest as mf
from . import runner
from .compare import TIER_NAMES
from .report import SessionResults
from .session import EntryRun


def _data_dir(args):
    if args.data_dir:
        return Path(args.data_dir)
    return mf.data_dir_from_env()


def cmd_draft(args):
    data_dir = _data_dir(args)
    out_path = Path(args.output) if args.output else mf.manifest_path(data_dir)
    if out_path.exists() and not args.force:
        print("%s exists; use --force to overwrite or --output" % out_path,
              file=sys.stderr)
        return 2
    manifest, notes = mf.draft_manifest(data_dir)
    errors = mf.validate_manifest(manifest)
    if errors:
        print("The draft is not a valid manifest; nothing written:",
              file=sys.stderr)
        for e in errors:
            print("  - " + e, file=sys.stderr)
        for n in notes:
            print("  note: " + n, file=sys.stderr)
        return 2
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(".json.tmp")
    import json
    tmp.write_text(json.dumps(manifest, indent=2) + "\n")
    tmp.replace(out_path)
    nfiles = sum(len(d["files"]) for d in manifest["datasets"].values())
    print("Wrote %s: %d dataset(s), %d entries, %d checksums" % (
        out_path, len(manifest["datasets"]), len(manifest["entries"]),
        nfiles))
    for e in manifest["entries"]:
        print("  %-14s profile %-5s plot_level %d  %s" % (
            e["id"], e["profile"], e["plot_level"],
            e["provenance"].get("drp_git_version")
            or e["provenance"].get("drp_version") or "(no provenance)"))
    if notes:
        print("Could not determine / check by hand:")
        for n in notes:
            print("  - " + n)
    return 0


def cmd_verify(args):
    data_dir = _data_dir(args)
    mfst = mf.load_manifest(data_dir)
    entries = [mf.get_entry(mfst, i) for i in args.entry] if args.entry \
        else mfst["entries"]
    bad = 0
    for e in entries:
        problems = mf.verify_checksums(data_dir, mfst, e)
        print("%-16s %s" % (e["id"], "ok" if not problems else
                            "%d problem(s)" % len(problems)))
        for p in problems[:10]:
            print("    " + p)
        bad += bool(problems)
    return 1 if bad else 0


def cmd_run(args):
    data_dir = _data_dir(args)
    mfst = mf.load_manifest(data_dir)
    entries = [mf.get_entry(mfst, i) for i in args.entry] if args.entry \
        else mf.select_entries(mfst, args.profile)
    session_dir = runner.new_session_dir(args.out)
    results = SessionResults(session_dir, {
        "profile": args.profile, "data_dir": str(data_dir),
        "drp": runner.git_describe()})
    results.write()
    failed = False
    for ent in entries:
        print("== %s" % ent["id"], flush=True)
        run = EntryRun(data_dir, mfst, ent, results, session_dir,
                       timeout=args.timeout, keep_runs=args.keep_runs)
        run.start()
        for t in sorted(TIER_NAMES):
            r = run.tier(t)
            print("   tier %d %-14s %-7s %s" % (t, r.name, r.status,
                                               r.summary[:100]))
            failed |= r.status == "fail"
        run.finish()
    print("report: %s" % (session_dir / "report.md"))
    return 1 if failed else 0


def cmd_bless(args):
    try:
        blessmod.bless(_data_dir(args), args.entry, args.reason, args.run,
                       args.out, args.yes)
    except blessmod.BlessError as exc:
        print("bless refused: %s" % exc, file=sys.stderr)
        return 2
    return 0


def cmd_report(args):
    session = Path(args.run) if args.run else runner.latest_session_dir(
        args.out)
    if session is None:
        print("no session found", file=sys.stderr)
        return 2
    SessionResults.load(session).write()
    print(session / "report.md")
    return 0


def build_parser():
    p = argparse.ArgumentParser(prog="python -m kcwidrp.tests.regression",
                                description=__doc__)
    p.add_argument("--data-dir", help="default: $KCWI_REGRESSION_DATA")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("draft-manifest", help="draft manifest.json from raw/")
    s.add_argument("--output", help="default: <data dir>/manifest.json")
    s.add_argument("--force", action="store_true")
    s.set_defaults(func=cmd_draft)

    s = sub.add_parser("verify-checksums", help="tier 0 for all entries")
    s.add_argument("--entry", action="append")
    s.set_defaults(func=cmd_verify)

    s = sub.add_parser("run", help="run entries without pytest")
    s.add_argument("--entry", action="append")
    s.add_argument("--profile", choices=mf.PROFILES, default="quick")
    s.add_argument("--out", help="default: $KCWI_REGRESSION_OUT or a "
                   "directory under the system temp dir")
    s.add_argument("--keep-runs", action="store_true")
    s.add_argument("--timeout", type=int, default=runner.DEFAULT_TIMEOUT)
    s.set_defaults(func=cmd_run)

    s = sub.add_parser("bless", help="make a run the new reference")
    s.add_argument("--entry", required=True)
    s.add_argument("--reason", help="'<reason>: <text>'; reason from the "
                   "plan's T.7 bless table")
    s.add_argument("--run", help="session or run directory (default: latest "
                   "session under $KCWI_REGRESSION_OUT)")
    s.add_argument("--out", help=argparse.SUPPRESS)
    s.add_argument("--yes", action="store_true", help="do not ask")
    s.set_defaults(func=cmd_bless)

    s = sub.add_parser("report", help="rewrite report.md from results.json")
    s.add_argument("--run", help="session directory (default: latest)")
    s.add_argument("--out", help=argparse.SUPPRESS)
    s.set_defaults(func=cmd_report)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except mf.ManifestError as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
