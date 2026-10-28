"""One entry's run and its checks, shared by the pytest fixture and the CLI."""
import shutil
from pathlib import Path

from . import compare, runner
from . import manifest as mf


class EntryRun:
    """Runs one manifest entry once and serves the tier results lazily.

    ``finish()`` tears the run down: ``redux/`` is deleted when every
    required tier passed (``report.DELETE_REQUIRES``) and ``keep_runs`` is
    false. Logs, the prompt record, ``reduce.out``, ``kcwi.cfg``, ``files.txt``
    and ``entry_report.md`` are always kept.
    """

    def __init__(self, data_dir, mfst, entry, results, session_dir,
                 command=None, timeout=runner.DEFAULT_TIMEOUT,
                 keep_runs=False, extra_env=None, free_fn=None):
        self.data_dir = Path(data_dir)
        self.mfst = mfst
        self.entry = entry
        self.results = results
        self.run_dir = Path(session_dir) / entry["id"]
        self.command = command
        self.timeout = timeout
        self.keep_runs = keep_runs
        self.extra_env = extra_env
        self.free_fn = free_fn
        self.res = None
        self._cmp = None
        self._tiers = {}
        self.redux_deleted = False

    # -- lifecycle ---------------------------------------------------------
    def start(self):
        self.res = runner.run_entry(
            self.data_dir, self.mfst, self.entry, self.run_dir,
            command=self.command, timeout=self.timeout,
            extra_env=self.extra_env, free_fn=self.free_fn)
        self.results.set_run(self.entry["id"], self.res)
        self.results.write()
        return self

    @property
    def ref_dir(self):
        return mf.reference_dir(self.data_dir, self.mfst, self.entry)

    @property
    def has_reference(self):
        return (self.ref_dir / "redux").is_dir()

    @property
    def ran(self):
        return self.res.status in ("ok", "timeout", "error")

    def finish(self):
        eid = self.entry["id"]
        if (self.ran and self.results.all_passed(eid) and not self.keep_runs
                and not self.results.any_failed(eid)):
            shutil.rmtree(self.run_dir / "redux", ignore_errors=True)
            self.redux_deleted = True
        self.results.set_outputs_deleted(eid, self.redux_deleted)
        self.results.write()
        if self.run_dir.is_dir():
            from .report import render_markdown
            sub = dict(self.results.data)
            sub["entries"] = {eid: self.results.data["entries"].get(eid, {})}
            (self.run_dir / "entry_report.md").write_text(
                render_markdown(sub))

    # -- tiers -------------------------------------------------------------
    def comparison(self):
        if self._cmp is None:
            self._cmp = compare.compare_outputs(
                self.run_dir, self.ref_dir, runner.procfile_name(self.entry),
                self.entry["tolerances"], self.run_dir / "diff_images")
        return self._cmp

    def answers_path(self):
        ans = self.entry.get("answers")
        return self.data_dir / "answers" / ans if ans else None

    def tier(self, n):
        if n not in self._tiers:
            result = self._compute(n)
            compare.apply_informational(result, self.entry)
            self._tiers[n] = result
            self.results.set_tier(self.entry["id"], result)
            self.results.write()
        return self._tiers[n]

    def _skipped(self, n, why):
        return compare.TierResult(n, "skipped", why)

    def _compute(self, n):
        st = self.res.status
        if n == 0:
            if st == "checksum_failed":
                return compare.TierResult(0, "fail", self.res.reason,
                                          self.res.problems)
            if st == "skipped":
                return self._skipped(0, self.res.reason)
            return compare.TierResult(0, "pass", "raw checksums match")
        if st in ("checksum_failed", "skipped"):
            return self._skipped(n, "not run: tier 0 %s" % (
                "failed" if st == "checksum_failed" else "skipped (%s)" %
                self.res.reason))
        if st == "build_failed":
            if n == 1:
                return compare.TierResult(1, "fail", self.res.reason,
                                          [self.res.reason])
            return self._skipped(n, "run directory could not be built")
        if n == 1:
            return compare.tier1_result(
                self.run_dir, self.entry, self.answers_path(),
                self.ref_dir, self.res)
        if n == 5:
            return compare.tier5_result(self.run_dir, self.ref_dir,
                                        self.entry)
        if n == 7:
            return compare.tier7_result(self.run_dir, self.res.wall_time)
        if not self.has_reference:
            return self._skipped(
                n, "no reference at %s (bless with reason 'initial')"
                % self.ref_dir)
        if self.redux_deleted:
            return self._skipped(n, "outputs already deleted")
        cmp_ = self.comparison()
        return {2: compare.structure_result, 3: compare.pixels_result,
                4: compare.headers_result, 6: compare.plots_result}[n](cmp_)
