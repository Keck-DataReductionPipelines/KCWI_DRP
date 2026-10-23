"""Session results (``results.json``), the Markdown report and diff images (T.7)."""
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np

from .compare import TIER_NAMES

# Tiers that must all pass before a run's redux/ is deleted at teardown.
DELETE_REQUIRES = (0, 1, 2, 3, 4, 6)
DETAIL_LINES = 10
STATUS_TEXT = {"pass": "pass", "fail": "**FAIL**", "info": "info",
               "skipped": "skip"}


def write_diff_image(run_plane, ref_plane, path, title):
    """Difference image (run - reference) of one 2-D plane, Agg only."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    diff = np.asarray(run_plane, dtype=np.float64) - np.asarray(
        ref_plane, dtype=np.float64)
    finite = np.isfinite(diff)
    lim = float(np.max(np.abs(diff[finite]))) if finite.any() else 1.0
    lim = lim or 1.0
    fig = Figure(figsize=(7, 5))
    FigureCanvasAgg(fig)
    ax = fig.add_subplot(111)
    im = ax.imshow(diff, origin="lower", aspect="auto", cmap="RdBu_r",
                   vmin=-lim, vmax=lim, interpolation="nearest")
    fig.colorbar(im, ax=ax, label="run - reference")
    ax.set_title(title, fontsize=8)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(path), dpi=100)


class SessionResults:
    """Machine-readable results of one session, written after every entry."""

    def __init__(self, session_dir, meta=None):
        self.session_dir = Path(session_dir)
        self.data = {"session": str(self.session_dir),
                     "created": datetime.now().isoformat(timespec="seconds"),
                     "meta": meta or {}, "entries": {}}

    @classmethod
    def load(cls, session_dir):
        obj = cls(session_dir)
        with open(Path(session_dir) / "results.json") as f:
            obj.data = json.load(f)
        return obj

    def _entry(self, entry_id):
        return self.data["entries"].setdefault(
            entry_id, {"run": None, "tiers": {}, "outputs_deleted": False})

    def set_run(self, entry_id, run_result):
        self._entry(entry_id)["run"] = run_result.to_dict()

    def set_tier(self, entry_id, result):
        self._entry(entry_id)["tiers"][str(result.tier)] = result.to_dict()

    def set_outputs_deleted(self, entry_id, deleted):
        self._entry(entry_id)["outputs_deleted"] = bool(deleted)

    def tier_status(self, entry_id, tier):
        rec = self.data["entries"].get(entry_id, {}).get("tiers", {})
        return rec.get(str(tier), {}).get("status")

    def all_passed(self, entry_id):
        return all(self.tier_status(entry_id, t) == "pass"
                   for t in DELETE_REQUIRES)

    def any_failed(self, entry_id):
        rec = self.data["entries"].get(entry_id, {}).get("tiers", {})
        return any(t["status"] == "fail" for t in rec.values())

    def diff_summary(self, entry_id):
        rec = self.data["entries"].get(entry_id, {}).get("tiers", {})
        return {k: {"status": v["status"], "summary": v["summary"]}
                for k, v in sorted(rec.items())}

    def write(self):
        self.session_dir.mkdir(parents=True, exist_ok=True)
        for name, text in (("results.json",
                            json.dumps(self.data, indent=2, default=str)),
                           ("report.md", render_markdown(self.data))):
            tmp = self.session_dir / (name + ".tmp")
            tmp.write_text(text)
            os.replace(tmp, self.session_dir / name)


def render_markdown(data):
    lines = ["# KCWI DRP regression report", ""]
    meta = data.get("meta", {})
    lines.append("Session `%s`, created %s." % (data["session"],
                                                data.get("created")))
    for key in ("profile", "data_dir", "drp", "python"):
        if key in meta:
            lines.append("- %s: `%s`" % (key, meta[key]))
    lines.append("")
    tiers = sorted(TIER_NAMES)
    lines.append("## Summary")
    lines.append("")
    lines.append("| entry | " + " | ".join(
        "%d %s" % (t, TIER_NAMES[t]) for t in tiers) + " |")
    lines.append("|---|" + "---|" * len(tiers))
    for eid, rec in sorted(data["entries"].items()):
        cells = [STATUS_TEXT.get(rec["tiers"].get(str(t), {}).get("status"),
                                 "-") for t in tiers]
        lines.append("| %s | %s |" % (eid, " | ".join(cells)))
    lines.append("")
    lines.append("pass = equal to the reference; info = reported, not failed; "
                 "skip = not applicable (for example no reference yet); "
                 "- = not run.")
    lines.append("")
    lines.append("## Details")
    for eid, rec in sorted(data["entries"].items()):
        lines.append("")
        lines.append("### %s" % eid)
        run = rec.get("run") or {}
        lines.append("")
        lines.append("- run: status `%s`%s; exit code %s (recorded, never "
                     "checked); wall time %s; run directory `%s`" % (
                         run.get("status"),
                         " (%s)" % run["reason"] if run.get("reason") else "",
                         run.get("exit_code"),
                         "%.0f s" % run["wall_time"]
                         if run.get("wall_time") is not None else "n/a",
                         run.get("run_dir")))
        prov = run.get("provenance") or {}
        if prov:
            lines.append("- code: %s, version %s" % (
                prov.get("drp_git_version"), prov.get("drp_version")))
        lines.append("- outputs: %s" % (
            "deleted after a passing run" if rec.get("outputs_deleted")
            else "kept"))
        for p in run.get("problems", [])[:DETAIL_LINES]:
            lines.append("- raw file problem: %s" % p)
        for t in tiers:
            res = rec["tiers"].get(str(t))
            if not res:
                continue
            lines.append("- tier %d %s: %s. %s" % (
                t, res["name"], STATUS_TEXT.get(res["status"]),
                res["summary"]))
            show = res["details"][:DETAIL_LINES]
            for d in show:
                lines.append("    - %s" % d)
            more = res.get("n_details", len(res["details"])) - len(show)
            if more > 0:
                lines.append("    - ... and %d more" % more)
            for f in res.get("data", {}).get("findings", []):
                if f.get("diff_png"):
                    lines.append("    - difference image: `%s`"
                                 % f["diff_png"])
    lines.append("")
    return "\n".join(lines)
