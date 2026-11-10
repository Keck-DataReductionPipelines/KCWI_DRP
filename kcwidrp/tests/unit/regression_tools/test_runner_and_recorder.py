"""The runner, the input recorder (sitecustomize) and tier 1 end to end.

Each test starts a small subprocess; the "reduction" is the fake script in
conftest.py. The recorder under test is the real sitecustomize.py.
"""
import json
import subprocess
import sys

import pytest

from kcwidrp.tests.regression import compare, runner
from kcwidrp.tests.regression import manifest as mf

RECORDER_DIR = str(runner.SITECUSTOMIZE_DIR)


def run(world, command, entry_id="ds9_blue_L1", where="run", **kw):
    entry = mf.get_entry(world.mfst, entry_id)
    return entry, runner.run_entry(world.data, world.mfst, entry,
                                   world.root / where / entry_id,
                                   command=command, **kw)


# -- run directory -----------------------------------------------------------

def test_run_directory_files_and_config(world):
    entry, res = run(world, world.command, where="rundir")
    d = res.run_dir
    ds = world.data / "raw" / "ds9_fake"
    # source: analytic -- files.txt is the list with lev0/ rewritten against
    # the dataset dir (what `sed "s|^lev0/|$DATA/$DS/lev0/|"` did in S.3)
    assert (d / "files.txt").read_text().splitlines() == [
        str(ds / "lev0/KB.1.fits"), str(ds / "lev0/KB.2.fits")]
    assert (d / "kcwi.sky").read_text() == (ds / "kcwi.sky").read_text()
    cfg = (d / "kcwi.cfg").read_text().splitlines()
    shipped = runner.packaged_config_text().splitlines()
    changed = {a.split("=")[0].strip() for a, b in zip(shipped, cfg) if a != b}
    # the packaged default has plot_level = 1 already, so only the other two
    # named keys differ from it
    assert changed == {"enable_plot_viewer", "plot_pause"}
    assert "enable_plot_viewer = False" in cfg and "plot_pause = 0" in cfg
    assert "plot_level = 1" in cfg
    assert len(cfg) == len(shipped)
    assert res.exit_code == 0 and res.wall_time > 0
    assert (d / "reduce.out").is_file()


def test_unknown_config_key_is_an_error():
    with pytest.raises(runner.RunError, match="no_such_key"):
        runner.apply_config_overrides("a = 1\n", {"no_such_key": 2})


def test_override_rewrites_only_the_named_line():
    text = "# plot_level = 9\nplot_level = 1\nother = 5\n"
    assert runner.apply_config_overrides(text, {"plot_level": 3}) == \
        "# plot_level = 9\nplot_level = 3\nother = 5\n"


def test_subprocess_environment(world, tmp_path, monkeypatch):
    monkeypatch.setenv("DISPLAY", ":0")
    monkeypatch.setenv("PYTHONPATH", "/existing")
    script = tmp_path / "envdump.py"
    script.write_text(
        "import json, os, sys\n"
        "json.dump({k: os.environ.get(k) for k in ('DISPLAY', 'PYTHONPATH',"
        " 'KCWI_PROMPT_RECORD', 'KCWI_STRICT_INPUT')}, open('env.json','w'))\n"
        "sys.stdout.write(sys.stdin.read())\n")
    entry, res = run(world, lambda e, d: [sys.executable, str(script)],
                     where="envrun")
    env = json.loads((res.run_dir / "env.json").read_text())
    assert env["DISPLAY"] is None
    assert env["PYTHONPATH"] == RECORDER_DIR + ":/existing"
    assert env["KCWI_PROMPT_RECORD"] == str(res.run_dir / "prompt_record.jsonl")
    assert env["KCWI_STRICT_INPUT"] == "1"
    # level 1: stdin is /dev/null, so nothing was read
    assert (res.run_dir / "reduce.out").read_bytes() == b""


def test_level3_stdin_is_the_answers_file(world, tmp_path):
    script = tmp_path / "cat.py"
    script.write_text("import sys; sys.stdout.write(sys.stdin.read())\n")
    entry, res = run(world, lambda e, d: [sys.executable, str(script)],
                     "ds9_blue_L3", where="stdin")
    assert (res.run_dir / "reduce.out").read_text() == "a\nb\nc\n"


# -- checksums, space, timeout -------------------------------------------------

def test_changed_raw_byte_fails_checksum_before_the_reduction(make_world):
    w = make_world()
    raw = w.data / "raw" / "ds9_fake" / "lev0" / "KB.2.fits"
    data = bytearray(raw.read_bytes())
    data[10] ^= 0x01                      # change exactly one byte
    raw.write_bytes(bytes(data))
    calls = []
    entry = mf.get_entry(w.mfst, "ds9_blue_L1")
    res = runner.run_entry(w.data, w.mfst, entry, w.root / "r",
                           command=lambda e, d: calls.append(1) or ["false"])
    assert res.status == "checksum_failed"
    assert "lev0/KB.2.fits" in res.problems[0]
    assert calls == []                     # the reduction was never started
    assert not (w.root / "r" / "files.txt").exists()


def test_missing_raw_file_fails_checksum(make_world):
    w = make_world()
    (w.data / "raw" / "ds9_fake" / "lev0" / "KB.1.fits").unlink()
    assert mf.verify_checksums(w.data, w.mfst, mf.get_entry(
        w.mfst, "ds9_blue_L1"))[0].endswith("missing")


def test_low_free_space_skips_the_entry(make_world):
    w = make_world()
    calls = []
    entry = mf.get_entry(w.mfst, "ds9_blue_L1")
    res = runner.run_entry(w.data, w.mfst, entry, w.root / "r",
                           command=lambda e, d: calls.append(1) or ["false"],
                           free_fn=lambda p: 1024)
    # no reference: 5 GiB (= 5.4e9 bytes) is required
    assert res.status == "skipped" and "required 5.4 GB" in res.reason
    assert not calls


def test_required_space_is_one_and_a_half_times_the_reference(world):
    entry = mf.get_entry(world.mfst, "ds9_blue_L1")
    size = runner.tree_size(world.data / "reference" / "ds9_blue_L1" / "redux")
    assert runner.required_free_bytes(world.data, world.mfst, entry) == \
        int(1.5 * size)


def test_timeout_is_reported(world, tmp_path):
    script = tmp_path / "sleep.py"
    script.write_text("import time; time.sleep(30)\n")
    entry, res = run(world, lambda e, d: [sys.executable, str(script)],
                     where="to", timeout=0.5)
    assert res.status == "timeout"


def test_channel_bug_path_check():
    # bug B8: 'kr' or 'kb' in a directory part triggers extra prompts
    assert runner.check_path_for_channel_bug(["/data/kbx/lev0/KB.1.fits"])
    assert not runner.check_path_for_channel_bug(["/data/ok/lev0/KB.1.fits"])


# -- the recorder ----------------------------------------------------------------

def py(code, record, stdin, extra_env=None, pythonpath=RECORDER_DIR):
    import os
    env = dict(os.environ, PYTHONPATH=pythonpath)
    env.pop("KCWI_PROMPT_RECORD", None)
    if record:
        env["KCWI_PROMPT_RECORD"] = str(record)
    env.update(extra_env or {})
    return subprocess.run([sys.executable, "-c", code], input=stdin, text=True,
                          capture_output=True, env=env)


def test_recorder_writes_prompt_and_answer_then_eof(tmp_path):
    rec = tmp_path / "rec.jsonl"
    code = ("a = input('first? ')\n"
            "try:\n    input('second? ')\nexcept EOFError:\n"
            "    print('EOF raised')\n")
    out = py(code, rec, "yes\n")
    assert "EOF raised" in out.stdout            # the EOFError is re-raised
    lines = [json.loads(x) for x in rec.read_text().splitlines()]
    # source: analytic -- one answered prompt, then end of stdin
    assert lines == [{"prompt": "first? ", "answer": "yes"},
                     {"prompt": "second? ", "eof": True}]


def test_recorder_does_nothing_when_unset(tmp_path):
    out = py("import builtins; print(hasattr(builtins.input,"
             " '_kcwi_recorder'))", None, "")
    assert out.stdout.strip() == "False"


def test_recorder_is_installed_when_set(tmp_path):
    out = py("import builtins; print(hasattr(builtins.input,"
             " '_kcwi_recorder'))", tmp_path / "r.jsonl", "")
    assert out.stdout.strip() == "True"


def test_recorder_chains_to_another_sitecustomize(tmp_path):
    other = tmp_path / "other"
    other.mkdir()
    (other / "sitecustomize.py").write_text(
        "import builtins; builtins.CHAINED = 'yes'\n")
    out = py("import builtins; print(getattr(builtins, 'CHAINED', 'no'))",
             tmp_path / "r.jsonl", "", pythonpath=RECORDER_DIR + ":" + str(other))
    assert out.stdout.strip() == "yes"


def test_iers_auto_download_is_switched_off_lazily(tmp_path):
    code = ("import sys\n"
            "assert 'astropy' not in sys.modules\n"      # not imported early
            "from astropy.utils import iers\n"
            "print(iers.conf.auto_download)\n")
    on = py(code, tmp_path / "r.jsonl", "")
    off = py("from astropy.utils import iers; print(iers.conf.auto_download)",
             None, "")
    assert on.stdout.strip() == "False", on.stderr
    assert off.stdout.strip() == "True"        # untouched without the record


# -- tier 1 end to end through the recorder ---------------------------------------

def level3(world, answers_lines, n_inputs, command_factory, tmp_path):
    entry = dict(mf.get_entry(world.mfst, "ds9_blue_L3"))
    entry.pop("answers_sha256", None)
    answers = world.data / "answers" / "variant.txt"
    answers.write_text("".join(x + "\n" for x in answers_lines))
    entry["answers"] = "variant.txt"
    cmd = command_factory(world.script, {3: n_inputs})
    res = runner.run_entry(world.data, world.mfst, entry,
                           tmp_path / ("r%d_%d" % (len(answers_lines), n_inputs)),
                           command=cmd)
    assert res.status == "ok"
    ref = world.data / "reference" / "ds9_blue_L3"
    return compare.tier1_result(res.run_dir, entry, answers, ref, res), res


def test_level3_exact_answers_pass_and_record_matches(world, command_factory,
                                                      tmp_path):
    r, res = level3(world, ["a", "b", "c"], 3, command_factory, tmp_path)
    assert r.status == "pass", r.details
    rec = compare.read_prompt_record(res.run_dir / "prompt_record.jsonl")
    assert [x["answer"] for x in rec] == ["a", "b", "c"]
    assert [x["prompt"] for x in rec] == ["prompt %d: " % i for i in range(3)]


def test_answers_file_one_line_short_fails_tier1_on_eof(world, command_factory,
                                                        tmp_path):
    # the fake asks 3 questions, the answers file has 2 lines
    r, _ = level3(world, ["a", "b"], 3, command_factory, tmp_path)
    assert r.status == "fail" and "EOF at prompt 3" in r.details[0]


def test_answers_file_one_line_extra_fails_tier1_on_count(world, command_factory,
                                                          tmp_path):
    r, _ = level3(world, ["a", "b", "c", "d"], 3, command_factory, tmp_path)
    assert r.status == "fail"
    assert "prompt count 3 != 4 lines" in r.details[0]


def test_changed_prompt_text_fails_against_the_reference_record(
        world, command_factory, tmp_path):
    ref = world.data / "reference" / "ds9_blue_L3" / "prompt_record.jsonl"
    saved = ref.read_text()
    try:
        ref.write_text(saved.replace("prompt 1", "prompt one"))
        r, _ = level3(world, ["a", "b", "c"], 3, command_factory, tmp_path)
    finally:
        ref.write_text(saved)
    assert r.status == "fail" and "prompt strings differ" in r.details[0]


def test_primitive_failure_in_log_fails_tier1_though_exit_code_is_zero(
        world, command_factory, tmp_path):
    entry = mf.get_entry(world.mfst, "ds9_blue_L1")
    res = runner.run_entry(world.data, world.mfst, entry, tmp_path / "f",
                           command=command_factory(
                               world.script, fail_ids=("ds9_blue_L1",)))
    r = compare.tier1_result(res.run_dir, entry, run_result=res)
    assert res.exit_code == 0
    assert r.status == "fail" and r.data["failed_primitives"] == ["MakeCube"]
    assert r.data["exit_code"] == 0
