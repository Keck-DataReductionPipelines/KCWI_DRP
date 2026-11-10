"""Tests for kcwidrp.core.kcwi_interactive (plan task P.2).

Written from plan section 8.2 ("Interface fixed for parallel work"), not from
the implementation.  Input is scripted by patching ``builtins.input``; the
module is imported inside a fixture so that the file collects without it.
"""
import logging
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

# source: analytic - every expected return value below follows from the
# contract in plan 8.2 (quit words, help words, range current+1 .. n_items-1,
# defaults at end of input); there is no numerical reference data.

REPO_ROOT = Path(__file__).resolve().parents[4]


class Config:
    """Stand-in for ``self.config.instrument``."""

    def __init__(self, plot_level=2, plot_pause=0.25):
        self.plot_level = plot_level
        self.plot_pause = plot_pause


class Script:
    """Scripted ``input``: answers in order, EOFError when exhausted."""

    def __init__(self, answers):
        self.answers = list(answers)
        self.prompts = []

    def __call__(self, prompt=""):
        self.prompts.append(prompt)
        if not self.answers:
            raise EOFError
        return self.answers.pop(0)

    @property
    def consumed_all(self):
        return not self.answers


@pytest.fixture
def ki(monkeypatch):
    """The module under test, with a clean end-of-input state."""
    monkeypatch.delenv("KCWI_STRICT_INPUT", raising=False)
    import kcwidrp.core.kcwi_interactive as mod

    mod.reset_end_of_input()
    yield mod
    mod.reset_end_of_input()


@pytest.fixture
def script(monkeypatch):
    def install(*answers):
        s = Script(answers)
        monkeypatch.setattr("builtins.input", s)
        return s

    return install


@pytest.fixture
def log(caplog):
    caplog.set_level(logging.DEBUG, logger="KCWI")
    return logging.getLogger("KCWI")


def recs(caplog, level=None):
    out = [r for r in caplog.records if r.name == "KCWI"]
    if level is not None:
        out = [r for r in out if r.levelno == level]
    return out


def no_input(prompt=""):
    pytest.fail("input() called unexpectedly with prompt %r" % (prompt,))


def fparse(text):
    return float(text)


def positive(v):
    return None if v > 0 else "value must be positive"


# ---------------------------------------------------------------- pause

def test_pause_low_plot_level_sleeps_and_never_asks(ki, monkeypatch, log):
    monkeypatch.setattr("builtins.input", no_input)
    calls = []
    monkeypatch.setattr(time, "sleep", lambda s: calls.append(s))
    result = ki.pause(Config(plot_level=1, plot_pause=0.37), log, what="W")
    assert result is None
    # source: analytic - plot_pause is passed through unchanged
    assert calls == [0.37]


def test_pause_plot_level_zero_sleeps(ki, monkeypatch, log):
    monkeypatch.setattr("builtins.input", no_input)
    calls = []
    monkeypatch.setattr(time, "sleep", lambda s: calls.append(s))
    ki.pause(Config(plot_level=0, plot_pause=2.5), log)
    assert calls == [2.5]


def test_pause_interactive_prints_what_and_asks(ki, monkeypatch, script, capsys, log):
    s = script("")
    monkeypatch.setattr(time, "sleep", lambda s_: pytest.fail("slept"))
    assert ki.pause(Config(plot_level=2), log, what="Arc fit residuals") is None
    assert "Arc fit residuals" in capsys.readouterr().out
    assert s.prompts == ["Next? <cr> (? for help): "]


def test_pause_high_plot_level_asks(ki, monkeypatch, script, log):
    s = script("")
    monkeypatch.setattr(time, "sleep", lambda s_: pytest.fail("slept"))
    ki.pause(Config(plot_level=3), log)
    assert len(s.prompts) == 1


def test_pause_help_then_cr(ki, script, capsys, log):
    s = script("?", "help", "")
    ki.pause(Config(plot_level=2), log, what="W", help="LONG HELP TEXT")
    out = capsys.readouterr().out
    assert out.count("LONG HELP TEXT") == 2
    assert len(s.prompts) == 3 and s.consumed_all


def test_help_is_set_off_by_blank_lines(ki, script, capsys, log):
    script("?", "")
    ki.pause(Config(plot_level=2), log, help="LONG HELP TEXT")
    # the typed answer ends the prompt line, so one "\n" is a blank line
    assert capsys.readouterr().out == "\nLONG HELP TEXT\n\n"


def _terminal(monkeypatch, **env):
    monkeypatch.setattr(sys.stdout, "isatty", lambda: True)
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "xterm-256color")
    for k, v in env.items():
        monkeypatch.setenv(k, v)


def test_terminal_bolds_what_and_help_headers(ki, script, capsys, log,
                                               monkeypatch):
    _terminal(monkeypatch)
    s = script("?", "")
    ki.pause(Config(plot_level=2), log, what="W",
             help="What you see:\n  Blue: data.\nDefault:\n  <cr>.")
    out = capsys.readouterr().out
    assert "\033[1mW\033[0m\n" in out
    assert "\033[1mWhat you see:\033[0m\n  Blue: data.\n" in out
    assert "\033[1mDefault:\033[0m\n  <cr>." in out
    # the question passed to input() stays plain (readline width)
    assert all("\033" not in p for p in s.prompts)


@pytest.mark.parametrize("env", [{"NO_COLOR": "1"}, {"TERM": "dumb"}])
def test_no_bold_with_no_color_or_dumb_terminal(ki, script, capsys, log,
                                                monkeypatch, env):
    _terminal(monkeypatch, **env)
    script("?", "")
    ki.pause(Config(plot_level=2), log, what="W", help="What you see:")
    assert "\033" not in capsys.readouterr().out


def test_pause_does_not_log_info(ki, script, caplog, log):
    script("")
    ki.pause(Config(plot_level=2), log, what="W", who="Prim")
    assert recs(caplog, logging.INFO) == []


# ------------------------------------------------------------- ask_next

def test_ask_next_cr_is_next(ki, script, log):
    s = script("")
    assert ki.ask_next(log, "what", "help") == ki.NEXT
    assert s.prompts == ["Next? <cr>, q to quit (? for help): "]


def test_ask_next_prompt_with_n_items(ki, script, log):
    s = script("")
    ki.ask_next(log, "what", "help", n_items=10)
    assert s.prompts == ["Next? <cr>, <int> to jump, q to quit (? for help): "]


def test_ask_next_prints_what(ki, script, capsys, log):
    script("")
    ki.ask_next(log, "THE WHAT LINE", "help")
    assert "THE WHAT LINE" in capsys.readouterr().out


@pytest.mark.parametrize("answer", ["q", "Q", "quit", "QUIT", "Quit", "  q  ", "\tquit\n", " Q "])
def test_ask_next_quit_words(ki, script, log, answer):
    script(answer)
    assert ki.ask_next(log, "w", "h") == ki.QUIT


@pytest.mark.parametrize("answer", ["quiet", "sq", "Q1", "qq", "quits", "aqb"])
@pytest.mark.parametrize("n_items", [None, 10])
def test_ask_next_text_containing_q_is_not_quit(ki, script, capsys, log,
                                                answer, n_items):
    # source: analytic - plan 8.2 rule: B5, any input containing q used to
    # quit; now only q/quit quits and "other input containing a q
    # re-prompts".  The second answer, q, shows the question was asked again
    # (an exhausted script would also give NEXT, by the end-of-input default).
    s = script(answer, "q")
    assert ki.ask_next(log, "w", "h", n_items=n_items) == ki.QUIT
    assert len(s.prompts) == 2 and s.consumed_all
    assert capsys.readouterr().out.strip() != ""


def test_ask_next_other_text_is_next(ki, script, log):
    script("whatever")
    assert ki.ask_next(log, "w", "h") == ki.NEXT


def test_ask_next_help_then_valid_not_counted(ki, script, capsys, caplog, log):
    s = script("?", "help", "")
    assert ki.ask_next(log, "w", "LONG HELP", who="Prim") == ki.NEXT
    assert capsys.readouterr().out.count("LONG HELP") == 2
    assert len(s.prompts) == 3
    # help presses are not answers: only the single "next" is logged (DEBUG)
    assert len(recs(caplog)) == 1


def test_ask_next_next_logged_at_debug_not_info(ki, script, caplog, log):
    script("")
    ki.ask_next(log, "w", "h", who="ArcOffsets")
    assert recs(caplog, logging.INFO) == []
    assert len(recs(caplog, logging.DEBUG)) == 1


def test_ask_next_jump_returns_int(ki, script, log):
    script("7")
    result = ki.ask_next(log, "w", "h", n_items=10, current=3)
    assert result == 7 and isinstance(result, int)
    assert result != ki.NEXT and result != ki.QUIT


@pytest.mark.parametrize("number", [4, 9])
def test_ask_next_range_ends_are_valid(ki, script, log, number):
    # source: analytic - valid range is current+1 .. n_items-1 = 4 .. 9
    script(str(number))
    assert ki.ask_next(log, "w", "h", n_items=10, current=3) == number


@pytest.mark.parametrize("bad", ["3", "2", "0", "-1", "-5", "10", "11", "100"])
def test_ask_next_out_of_range_reprompts_with_range(ki, script, capsys, log, bad):
    # B4: <= current, negative, >= n_items used to end the review silently.
    s = script(bad, "5")
    assert ki.ask_next(log, "w", "h", n_items=10, current=3) == 5
    assert len(s.prompts) == 2
    out = capsys.readouterr().out
    # source: analytic - the message must give the valid range 4 .. 9
    assert "4" in out and "9" in out


def test_ask_next_out_of_range_not_logged_as_answer(ki, script, caplog, log):
    script("99", "5")
    ki.ask_next(log, "w", "h", n_items=10, current=3, who="P")
    assert len(recs(caplog, logging.INFO)) <= 1  # only the accepted jump


def test_ask_next_current_none_reads_as_minus_one(ki, script, log):
    # source: analytic - current=None -> valid range 0 .. n_items-1
    script("0")
    assert ki.ask_next(log, "w", "h", n_items=5, current=None) == 0
    script("4")
    assert ki.ask_next(log, "w", "h", n_items=5) == 4


def test_ask_next_current_none_rejects_negative_and_n_items(ki, script, log):
    s = script("-1", "5", "2")
    assert ki.ask_next(log, "w", "h", n_items=5, current=None) == 2
    assert len(s.prompts) == 3


def test_ask_next_jump_current_itself_rejected(ki, script, log):
    s = script("0", "1")
    assert ki.ask_next(log, "w", "h", n_items=5, current=0) == 1
    assert len(s.prompts) == 2


def test_ask_next_q_after_out_of_range(ki, script, log):
    script("50", "q")
    assert ki.ask_next(log, "w", "h", n_items=10, current=3) == ki.QUIT


# ------------------------------------------------------------ ask_value

def test_ask_value_valid(ki, script, caplog, log):
    s = script("12.5")
    v = ki.ask_value(log, "Offset", "h", fparse, 10.0, "offset set to {value}", who="ReadAtlas")
    assert v == 12.5
    assert s.prompts == ["Offset (? for help): "]
    # source: analytic - INFO line is who prefix + effect.format(value=...)
    infos = recs(caplog, logging.INFO)
    assert [r.getMessage() for r in infos] == ["ReadAtlas: offset set to 12.5"]


def test_ask_value_prompt_ending_colon_space_unchanged(ki, script, log):
    s = script("")
    ki.ask_value(log, "Enter value: ", "h", fparse, 1.0, "e {value}")
    assert s.prompts == ["Enter value: "]


def test_ask_value_prompt_without_colon_gets_help_suffix(ki, script, log):
    s = script("")
    ki.ask_value(log, "Enter value", "h", fparse, 1.0, "e {value}")
    assert s.prompts == ["Enter value (? for help): "]


def test_ask_value_empty_returns_default_without_info(ki, script, caplog, log):
    script("")
    v = ki.ask_value(log, "P: ", "h", fparse, 3.5, "set {value}", who="X")
    assert v == 3.5
    assert recs(caplog, logging.INFO) == []


def test_ask_value_no_who_has_no_prefix(ki, script, caplog, log):
    script("2")
    ki.ask_value(log, "P: ", "h", fparse, 1.0, "set to {value}")
    assert [r.getMessage() for r in recs(caplog, logging.INFO)] == ["set to 2.0"]


def test_ask_value_parse_failure_reprompts_logs_once(ki, script, capsys, caplog, log):
    s = script("abc", "4")
    v = ki.ask_value(log, "P: ", "h", fparse, 1.0, "set {value}", who="X")
    assert v == 4.0
    assert len(s.prompts) == 2
    assert capsys.readouterr().out.strip() != ""  # a reason was printed
    assert [r.getMessage() for r in recs(caplog, logging.INFO)] == ["X: set 4.0"]


def test_ask_value_validate_failure_prints_reason_and_reprompts(ki, script, capsys, caplog, log):
    s = script("-3", "0", "6")
    v = ki.ask_value(log, "P: ", "h", fparse, 1.0, "set {value}", validate=positive, who="X")
    assert v == 6.0
    assert len(s.prompts) == 3
    assert capsys.readouterr().out.count("value must be positive") == 2
    assert [r.getMessage() for r in recs(caplog, logging.INFO)] == ["X: set 6.0"]


def test_ask_value_parse_error_reason_is_printed(ki, script, capsys, log):
    def parse(text):
        if text != "ok":
            raise ValueError("needs exactly two numbers")
        return 1

    script("bad", "ok")
    ki.ask_value(log, "P: ", "h", parse, 0, "e {value}")
    assert "needs exactly two numbers" in capsys.readouterr().out


def test_ask_value_help_then_valid_not_counted(ki, script, capsys, caplog, log):
    s = script("?", "help", "8")
    v = ki.ask_value(log, "P: ", "LONG HELP", fparse, 1.0, "set {value}", who="X")
    assert v == 8.0
    assert capsys.readouterr().out.count("LONG HELP") == 2
    assert len(s.prompts) == 3
    assert len(recs(caplog, logging.INFO)) == 1


def test_ask_value_validate_none_accepts_anything_parsed(ki, script, log):
    script("-7")
    assert ki.ask_value(log, "P: ", "h", fparse, 1.0, "e {value}") == -7.0


def test_ask_value_validate_called_with_parsed_value(ki, script, log):
    seen = []

    def validate(v):
        seen.append(v)
        return None

    script("2.5")
    ki.ask_value(log, "P: ", "h", fparse, 1.0, "e {value}", validate=validate)
    assert seen == [2.5]


def test_ask_value_default_returned_as_is(ki, script, log):
    # An empty answer returns the default untouched, even if it would not
    # pass validation.
    script("")
    v = ki.ask_value(log, "P: ", "h", fparse, -1.0, "e", validate=positive)
    assert v == -1.0


# --------------------------------------------------------- ask_yes_quit

def test_yes_quit_cr_proceeds(ki, script, log):
    script("")
    assert ki.ask_yes_quit(log, "Proceed? <cr>=yes or Q=quit: ", "h") is True


@pytest.mark.parametrize("answer", ["q", "Q", "quit", "QuIt", "  q ", " quit "])
def test_yes_quit_quit_words(ki, script, log, answer):
    script(answer)
    assert ki.ask_yes_quit(log, "Proceed? ", "h") is False


@pytest.mark.parametrize("answer", ["quiet", "sq", "Q1"])
def test_yes_quit_text_containing_q_reprompts(ki, script, log, answer):
    # source: analytic - plan 8.2 rule: B5, text containing q must not quit;
    # it re-prompts.  The second answer, q, shows it was asked again.
    s = script(answer, "q")
    assert ki.ask_yes_quit(log, "Proceed? ", "h") is False
    assert len(s.prompts) == 2 and s.consumed_all


@pytest.mark.parametrize("answer", ["yes", "no", "y", "zzz"])
def test_yes_quit_other_text_proceeds_as_before(ki, script, log, answer):
    # source: analytic - D12: today any text without a q proceeds
    # (reduce_kcwi.py and kcwi_rti.py 'Proceed?' prompts), so it still does.
    s = script(answer)
    assert ki.ask_yes_quit(log, "Proceed? ", "h") is True
    assert len(s.prompts) == 1 and s.consumed_all


def test_yes_quit_help_then_quit(ki, script, capsys, log):
    s = script("?", "help", "q")
    assert ki.ask_yes_quit(log, "Proceed? ", "LONG HELP") is False
    assert capsys.readouterr().out.count("LONG HELP") == 2
    assert len(s.prompts) == 3


def test_yes_quit_invalid_then_valid_prints_reason(ki, script, capsys, log):
    script("quiet", "")
    ki.ask_yes_quit(log, "Proceed? ", "h")
    assert capsys.readouterr().out.strip() != ""


# ------------------------------------------------------- end of input

def test_eof_ask_next_warns_once_returns_next(ki, script, caplog, log):
    script()
    assert ki.ask_next(log, "w", "h", who="ArcOffsets") == ki.NEXT
    assert len(recs(caplog, logging.WARNING)) == 1


def test_eof_warning_names_the_prompt_ask_value(ki, script, caplog, log):
    script()
    v = ki.ask_value(log, "Enter atlas offset", "h", fparse, 9.0, "e {value}", who="ReadAtlas")
    assert v == 9.0
    warns = recs(caplog, logging.WARNING)
    assert len(warns) == 1
    # The warning must identify the prompt: its text or the primitive.
    msg = warns[0].getMessage()
    assert "atlas offset" in msg or "ReadAtlas" in msg


def test_eof_warning_names_the_prompt_yes_quit(ki, script, caplog, log):
    script()
    assert ki.ask_yes_quit(log, "Proceed? <cr>=yes or Q=quit", "h", who="Reduce") is True
    msg = recs(caplog, logging.WARNING)[0].getMessage()
    assert "Proceed" in msg or "Reduce" in msg


def test_eof_pause_proceeds(ki, script, caplog, log):
    script()
    assert ki.pause(Config(plot_level=2), log, what="W", who="P") is None
    assert len(recs(caplog, logging.WARNING)) == 1


def test_eof_then_all_later_calls_return_defaults_silently(ki, script, monkeypatch, caplog, log):
    script()
    assert ki.ask_next(log, "w", "h", who="A") == ki.NEXT
    assert len(recs(caplog, logging.WARNING)) == 1
    caplog.clear()

    # From here on input must not be called at all.
    monkeypatch.setattr("builtins.input", no_input)
    assert ki.ask_next(log, "w", "h", n_items=5, current=0) == ki.NEXT
    assert ki.ask_value(log, "P: ", "h", fparse, 2.5, "e {value}") == 2.5
    assert ki.ask_yes_quit(log, "Proceed? ", "h") is True
    assert ki.pause(Config(plot_level=2), log) is None
    assert recs(caplog, logging.WARNING) == []
    assert recs(caplog, logging.INFO) == []


@pytest.mark.parametrize("first", ["next", "value", "yes", "pause"])
def test_eof_in_any_function_sets_sticky_state(ki, script, monkeypatch, caplog, log, first):
    script()
    if first == "next":
        ki.ask_next(log, "w", "h")
    elif first == "value":
        ki.ask_value(log, "P: ", "h", fparse, 1.0, "e")
    elif first == "yes":
        ki.ask_yes_quit(log, "Proceed? ", "h")
    else:
        ki.pause(Config(plot_level=2), log)
    caplog.clear()
    monkeypatch.setattr("builtins.input", no_input)
    assert ki.ask_next(log, "w", "h") == ki.NEXT
    assert ki.ask_value(log, "P: ", "h", fparse, 7.0, "e") == 7.0
    assert recs(caplog, logging.WARNING) == []


def test_eof_after_some_answers_only_warns_at_exhaustion(ki, script, caplog, log):
    script("", "")
    ki.ask_next(log, "w", "h")
    ki.ask_next(log, "w", "h")
    assert recs(caplog, logging.WARNING) == []
    ki.ask_next(log, "w", "h")
    assert len(recs(caplog, logging.WARNING)) == 1


def test_reset_end_of_input_clears_state(ki, script, caplog, log):
    script()
    ki.ask_next(log, "w", "h")
    ki.reset_end_of_input()
    caplog.clear()
    s = script("q")
    # input is read again, so the scripted answer is honored
    assert ki.ask_next(log, "w", "h") == ki.QUIT
    assert len(s.prompts) == 1
    # and a new EOF warns again
    ki.reset_end_of_input()
    script()
    ki.ask_next(log, "w", "h")
    assert len(recs(caplog, logging.WARNING)) == 1


def test_eof_strict_raises_ask_next(ki, script, monkeypatch, log):
    monkeypatch.setenv("KCWI_STRICT_INPUT", "1")
    script()
    with pytest.raises(EOFError):
        ki.ask_next(log, "w", "h")


def test_eof_strict_raises_every_function(ki, script, monkeypatch, log):
    monkeypatch.setenv("KCWI_STRICT_INPUT", "1")
    script()
    with pytest.raises(EOFError):
        ki.ask_value(log, "P: ", "h", fparse, 1.0, "e")
    with pytest.raises(EOFError):
        ki.ask_yes_quit(log, "Proceed? ", "h")
    with pytest.raises(EOFError):
        ki.pause(Config(plot_level=2), log)


def test_eof_strict_raises_every_time_no_sticky_default(ki, script, monkeypatch, caplog, log):
    monkeypatch.setenv("KCWI_STRICT_INPUT", "1")
    script()
    for _ in range(2):
        with pytest.raises(EOFError):
            ki.ask_next(log, "w", "h")
    assert recs(caplog, logging.WARNING) == []


def test_strict_is_read_at_call_time(ki, script, monkeypatch, log):
    # The module is already imported (fixture); set the variable afterwards.
    script()
    assert ki.ask_next(log, "w", "h") == ki.NEXT  # default mode first
    ki.reset_end_of_input()
    monkeypatch.setenv("KCWI_STRICT_INPUT", "1")
    script()
    with pytest.raises(EOFError):
        ki.ask_next(log, "w", "h")
    monkeypatch.delenv("KCWI_STRICT_INPUT")
    ki.reset_end_of_input()
    script()
    assert ki.ask_next(log, "w", "h") == ki.NEXT


def test_strict_mode_does_not_affect_scripted_answers(ki, script, monkeypatch, log):
    monkeypatch.setenv("KCWI_STRICT_INPUT", "1")
    script("q")
    assert ki.ask_next(log, "w", "h") == ki.QUIT


def test_pause_low_plot_level_unaffected_by_eof_state(ki, script, monkeypatch, log):
    script()
    ki.ask_next(log, "w", "h")  # puts the process in end-of-input state
    calls = []
    monkeypatch.setattr(time, "sleep", lambda s: calls.append(s))
    ki.pause(Config(plot_level=1, plot_pause=0.5), log)
    # source: analytic - the sleep branch does not depend on input state
    assert calls == [0.5]


# ----------------------------------------------------------- misc rules

def test_input_looked_up_at_call_time(ki, monkeypatch, log):
    # The module is already imported; patching afterwards must take effect.
    seen = []

    def late(prompt=""):
        seen.append(prompt)
        return "quit"

    monkeypatch.setattr("builtins.input", late)
    assert ki.ask_next(log, "w", "h") == ki.QUIT
    assert len(seen) == 1

    seen2 = []
    monkeypatch.setattr("builtins.input", lambda p="": seen2.append(p) or "")
    assert ki.ask_next(log, "w", "h") == ki.NEXT
    assert len(seen2) == 1


def test_who_none_gives_no_prefix(ki, script, caplog, log):
    script("")
    ki.ask_next(log, "w", "h")
    for r in recs(caplog):
        assert not r.getMessage().startswith(": ")


def test_import_does_not_load_matplotlib_or_bokeh():
    code = (
        "import sys, kcwidrp.core.kcwi_interactive\n"
        "bad = sorted(m for m in sys.modules if m.startswith(('matplotlib', 'bokeh')))\n"
        "assert not bad, bad\n"
    )
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    proc = subprocess.run([sys.executable, "-c", code], cwd=str(REPO_ROOT),
                          env=env, capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, proc.stdout + proc.stderr
