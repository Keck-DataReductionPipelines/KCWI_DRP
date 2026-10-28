"""Build a run directory and reduce one manifest entry in a subprocess (T.7).

The reduction command is injectable: ``run_entry(..., command=callable)`` where
the callable maps ``(entry, run_dir)`` to an argument list. The default is
``default_command``; unit tests substitute a fake reduction.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

from . import manifest as mf

SITECUSTOMIZE_DIR = Path(__file__).resolve().parent / "sitecustomize_dir"
DEFAULT_TIMEOUT = 2 * 3600
MIN_FREE_NO_REFERENCE = 5 * 1024 ** 3
FREE_FACTOR = 1.5


class RunError(Exception):
    """The run directory could not be built."""


def default_out_root():
    """``$KCWI_REGRESSION_OUT`` or a fixed directory under the system temp dir."""
    env = os.environ.get("KCWI_REGRESSION_OUT")
    if env:
        return Path(env)
    return Path(tempfile.gettempdir()) / "kcwi_regression_out"


def new_session_dir(out_root=None):
    root = Path(out_root) if out_root else default_out_root()
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    path = root / stamp
    n = 1
    while path.exists():
        n += 1
        path = root / ("%s-%d" % (stamp, n))
    path.mkdir(parents=True)
    return path


def latest_session_dir(out_root=None):
    root = Path(out_root) if out_root else default_out_root()
    if not root.is_dir():
        return None
    sessions = sorted(p for p in root.iterdir()
                      if p.is_dir() and (p / "results.json").is_file())
    return sessions[-1] if sessions else None


def default_command(entry, run_dir):
    """The reduction exactly as S.3 ran it (docs/baseline/README.md)."""
    flag = mf.CHANNELS[entry["channel"]][0]
    return [sys.executable, "-c",
            "from kcwidrp.scripts.reduce_kcwi import main; main()",
            flag, "-c", "kcwi.cfg", "-l", "files.txt", "-g"]


def procfile_name(entry):
    return mf.CHANNELS[entry["channel"]][1]


# -- run directory ------------------------------------------------------------

# The viewer switch was renamed in PR1 (bokeh plan 5.1).  A manifest drafted
# on either side of the rename must run on code from either side, so an
# override of one name applies to whichever of the two the config has.
KEY_ALIASES = {"enable_bokeh": "enable_plot_viewer",
               "enable_plot_viewer": "enable_bokeh"}


def _key_pattern(key):
    return re.compile(r"^(\s*)%s\s*=.*$" % re.escape(key), re.M)


def apply_config_overrides(text, overrides):
    """Rewrite only the ``key = value`` lines named in ``overrides``."""
    for key, value in overrides.items():
        pattern = _key_pattern(key)
        if not pattern.search(text) and key in KEY_ALIASES:
            alias = KEY_ALIASES[key]
            if _key_pattern(alias).search(text):
                key, pattern = alias, _key_pattern(alias)
        if not pattern.search(text):
            raise RunError("config key %r not found in kcwi.cfg" % key)
        text = pattern.sub(lambda m: "%s%s = %s" % (m.group(1), key, value),
                           text, count=1)
    return text


def packaged_config_text():
    """The shipped default ``kcwi.cfg`` (what ``--write_config`` copies)."""
    from kcwidrp.core.kcwi_pkg_resources import get_resource_path
    return Path(get_resource_path("kcwidrp", "configs/kcwi.cfg")).read_text()


def entry_overrides(entry):
    ov = {"plot_level": entry["plot_level"]}
    ov.update(entry["config_overrides"])
    return ov


def check_path_for_channel_bug(paths):
    """B8: ``kr``/``kb`` in a directory part triggers extra prompts."""
    bad = []
    for p in paths:
        directory = os.path.dirname(p).lower()
        if "kr" in directory or "kb" in directory:
            bad.append(p)
    return bad


def build_run_dir(data_dir, mfst, entry, run_dir, config_text=None):
    """Write files.txt, kcwi.sky and kcwi.cfg into ``run_dir``."""
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    ds_dir = mf.dataset_dir(data_dir, mfst, entry["dataset"])
    lines = mf.read_file_list(ds_dir / entry["file_list"])
    absolute = []
    for rel in lines:
        absolute.append(rel if os.path.isabs(rel) else str(ds_dir / rel))
    bad = check_path_for_channel_bug(absolute)
    if bad:
        raise RunError(
            "raw path %s has 'kr' or 'kb' in a directory name; bug B8 would "
            "add prompts that consume answer lines. Move the data directory."
            % bad[0])
    (run_dir / "files.txt").write_text("\n".join(absolute) + "\n")
    if entry.get("sky_file"):
        shutil.copyfile(ds_dir / entry["sky_file"], run_dir / "kcwi.sky")
    text = config_text if config_text is not None else packaged_config_text()
    (run_dir / "kcwi.cfg").write_text(
        apply_config_overrides(text, entry_overrides(entry)))
    return run_dir


# -- sizes and space ---------------------------------------------------------------

def tree_size(path):
    total = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            try:
                total += os.stat(os.path.join(root, name)).st_size
            except OSError:
                pass
    return total


def required_free_bytes(data_dir, mfst, entry):
    ref = mf.reference_dir(data_dir, mfst, entry)
    redux = ref / "redux"
    if redux.is_dir():
        return int(FREE_FACTOR * tree_size(redux))
    return MIN_FREE_NO_REFERENCE


def free_bytes(path):
    return shutil.disk_usage(path).free


# -- provenance --------------------------------------------------------------------

def git_describe():
    """``git describe`` of the repository the package is imported from."""
    import kcwidrp
    repo = Path(kcwidrp.__file__).resolve().parent.parent
    try:
        out = subprocess.run(
            ["git", "describe", "--always", "--dirty", "--tags"], cwd=repo,
            capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() or None


def current_provenance():
    """DRP commit and version plus the numeric library versions (D18)."""
    prov = {"drp_git_version": git_describe(), "drp_version": None,
            "date": datetime.now().strftime("%Y-%m-%d"),
            "python": sys.version.split()[0]}
    try:
        from importlib.metadata import version
        prov["drp_version"] = version("kcwidrp")
    except Exception:
        pass
    for name in ("numpy", "scipy", "astropy", "scikit-image", "astroscrappy"):
        try:
            from importlib.metadata import version
            prov[name] = version(name)
        except Exception:
            prov[name] = None
    try:
        from astropy.utils import iers
        table = iers.IERS_B.open()
        prov["iers"] = {"IERS_B_last_mjd": float(table["MJD"][-1]),
                        "auto_download": bool(iers.conf.auto_download)}
    except Exception:
        prov["iers"] = None
    return prov


# -- running -----------------------------------------------------------------------

def reduction_env(run_dir, base_env=None, extra_env=None):
    env = dict(os.environ if base_env is None else base_env)
    env.pop("DISPLAY", None)
    old = env.get("PYTHONPATH")
    env["PYTHONPATH"] = str(SITECUSTOMIZE_DIR) + (
        os.pathsep + old if old else "")
    env["KCWI_PROMPT_RECORD"] = str(Path(run_dir) / "prompt_record.jsonl")
    env["KCWI_STRICT_INPUT"] = "1"
    if extra_env:
        env.update(extra_env)
    return env


class RunResult:
    """What happened when one entry was run."""

    def __init__(self, entry_id, run_dir):
        self.entry_id = entry_id
        self.run_dir = Path(run_dir)
        # ok | checksum_failed | skipped | build_failed | timeout | error
        self.status = "ok"
        self.reason = ""
        self.problems = []
        self.exit_code = None
        self.wall_time = None
        self.command = None
        self.provenance = None

    def to_dict(self):
        return {"status": self.status, "reason": self.reason,
                "problems": self.problems[:50], "exit_code": self.exit_code,
                "wall_time": self.wall_time, "command": self.command,
                "run_dir": str(self.run_dir), "provenance": self.provenance}


def run_entry(data_dir, mfst, entry, run_dir, command=None,
              timeout=DEFAULT_TIMEOUT, extra_env=None, free_fn=None,
              config_text=None):
    """Verify checksums, check space, build the run directory, reduce.

    Never raises for an expected failure; the outcome is in ``status``.
    """
    command = command or default_command
    free_fn = free_fn or free_bytes
    run_dir = Path(run_dir)
    res = RunResult(entry["id"], run_dir)
    problems = mf.verify_checksums(data_dir, mfst, entry)
    if problems:
        res.status = "checksum_failed"
        res.reason = "%d raw file problem(s)" % len(problems)
        res.problems = problems
        return res
    run_dir.mkdir(parents=True, exist_ok=True)
    need = required_free_bytes(data_dir, mfst, entry)
    have = free_fn(run_dir)
    if have < need:
        res.status = "skipped"
        res.reason = "free space %.1f GB < required %.1f GB" % (
            have / 1e9, need / 1e9)
        return res
    try:
        build_run_dir(data_dir, mfst, entry, run_dir, config_text)
    except (RunError, OSError) as exc:
        res.status = "build_failed"
        res.reason = str(exc)
        return res
    res.provenance = current_provenance()
    argv = list(command(entry, run_dir))
    res.command = argv
    answers = entry.get("answers") if entry["plot_level"] == 3 else None
    stdin = (open(Path(data_dir) / "answers" / answers, "rb") if answers
             else open(os.devnull, "rb"))
    start = time.monotonic()
    try:
        with stdin, open(run_dir / "reduce.out", "wb") as out:
            proc = subprocess.run(
                argv, cwd=run_dir, stdin=stdin, stdout=out,
                stderr=subprocess.STDOUT, timeout=timeout,
                env=reduction_env(run_dir, extra_env=extra_env))
        res.exit_code = proc.returncode
    except subprocess.TimeoutExpired:
        res.status = "timeout"
        res.reason = "no exit after %d s" % timeout
    except OSError as exc:
        res.status = "error"
        res.reason = "could not start the reduction: %s" % exc
    res.wall_time = time.monotonic() - start
    (run_dir / "run.json").write_text(json.dumps(res.to_dict(), indent=2))
    return res
