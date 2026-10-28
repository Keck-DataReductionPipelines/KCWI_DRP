"""``what=`` of ``ask_value`` and ``ask_yes_quit`` (plan 8.2, "Interface
fixed for parallel work": "``what``, where given, is printed once before the
first question, as for ``pause`` and ``ask_next``").

Input is scripted by replacing ``builtins.input``; the fake prints a marker
line when it is asked, so the order of the ``what`` text and the questions
can be read from the captured standard output.  Nothing here reads the
helper's source.
"""
import builtins
import logging

import pytest

from kcwidrp.core import kcwi_interactive as ki

WHAT = "WHAT-TEXT: look at the plot and decide"
HELP = "HELP-TEXT: long explanation"
ASKED = "<<asked>>"
LOG = logging.getLogger("test_helper_what")


@pytest.fixture(autouse=True)
def fresh_state(monkeypatch):
    monkeypatch.delenv("KCWI_STRICT_INPUT", raising=False)
    ki.reset_end_of_input()
    yield
    ki.reset_end_of_input()


@pytest.fixture
def script(monkeypatch):
    def install(*answers):
        left = list(answers)

        def fake_input(prompt=""):
            print(ASKED)
            if not left:
                raise AssertionError("asked more often than scripted")
            return left.pop(0)

        monkeypatch.setattr(builtins, "input", fake_input)
        return left
    return install


def ask_value(**kw):
    return ki.ask_value(LOG, "Offset? ", HELP, parse=int, default=None,
                        effect="set to {value}", who="Prim",
                        validate=lambda v: "negative" if v < 0 else None,
                        **kw)


def ask_yes_quit(**kw):
    return ki.ask_yes_quit(LOG, "Proceed? ", HELP, who="Prim", **kw)


def assert_once_before_first_question(out):
    # source: analytic - plan 8.2: printed once, before the first question
    assert out.count(WHAT) == 1
    assert out.index(WHAT) < out.index(ASKED)


# ------------------------------------------------------------- ask_value
def test_ask_value_prints_what_once_before_the_question(script, capsys):
    script("12")
    assert ask_value(what=WHAT) == 12
    assert_once_before_first_question(capsys.readouterr().out)


def test_ask_value_prints_what_once_when_a_question_is_repeated(script,
                                                                 capsys):
    # "?" (help), a parse failure and a refusal each ask again; what stays
    # printed once
    answers = script("?", "abc", "1.5", "-7", "")
    assert ask_value(what=WHAT,) is None       # <cr>: the default
    out = capsys.readouterr().out
    assert_once_before_first_question(out)
    assert out.count(ASKED) == 5 and not answers
    assert HELP in out


def test_ask_value_what_also_printed_when_first_answer_is_default(script,
                                                                   capsys):
    script("")
    ask_value(what=WHAT)
    assert_once_before_first_question(capsys.readouterr().out)


def test_ask_value_without_what_prints_no_what(script, capsys):
    script("?", "4")
    assert ask_value(what=None) == 4
    out = capsys.readouterr().out
    # source: analytic - only the help text and the markers are printed
    assert WHAT not in out and "None" not in out
    assert out.count(ASKED) == 2


def test_ask_value_default_for_what_is_none(script, capsys):
    script("4")
    assert ki.ask_value(LOG, "Offset? ", HELP, parse=int, default=None,
                        effect="set to {value}") == 4
    assert WHAT not in capsys.readouterr().out


def test_ask_value_what_not_logged(script, caplog):
    # source: plan 8.2 - explanations go to stdout and are not logged
    script("3")
    with caplog.at_level(logging.DEBUG, logger=LOG.name):
        ask_value(what=WHAT)
    assert WHAT not in caplog.text


# ----------------------------------------------------------- ask_yes_quit
def test_ask_yes_quit_prints_what_once_before_the_question(script, capsys):
    script("")
    assert ask_yes_quit(what=WHAT) is True
    assert_once_before_first_question(capsys.readouterr().out)


def test_ask_yes_quit_prints_what_once_when_a_question_is_repeated(script,
                                                                    capsys):
    # "?" and "eq" (contains a q: asks again, B5) repeat the question
    answers = script("?", "eq", "q")
    assert ask_yes_quit(what=WHAT) is False
    out = capsys.readouterr().out
    assert_once_before_first_question(out)
    assert out.count(ASKED) == 3 and not answers
    assert HELP in out


def test_ask_yes_quit_without_what_prints_no_what(script, capsys):
    script("?", "")
    assert ask_yes_quit(what=None) is True
    out = capsys.readouterr().out
    assert WHAT not in out and "None" not in out
    assert out.count(ASKED) == 2


def test_ask_yes_quit_default_for_what_is_none(script, capsys):
    script("")
    assert ki.ask_yes_quit(LOG, "Proceed? ", HELP) is True
    assert capsys.readouterr().out.count(ASKED) == 1


# -------------------------------------------------------- end of input
def test_what_is_printed_once_when_input_ends(monkeypatch, capsys):
    # EOF takes the default (D4); what was still shown before the question
    def eof(prompt=""):
        print(ASKED)
        raise EOFError

    monkeypatch.setattr(builtins, "input", eof)
    assert ask_value(what=WHAT) is None
    assert_once_before_first_question(capsys.readouterr().out)
