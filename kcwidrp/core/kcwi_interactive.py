"""Interactive prompts for the primitives and launch scripts (Track P).

Every prompt the pipeline asks goes through one of four functions:

* ``pause(config, logger, what=None, help=None, who=None)``: the
  "look at this plot" pause.  At ``plot_level >= 2`` it asks
  ``Next? <cr> (? for help): ``; below that it sleeps ``plot_pause``
  seconds and reads no input.
* ``ask_next(logger, what, help, n_items=None, current=None, who=None)``:
  review loops.  Returns ``NEXT``, ``QUIT`` or the item number to jump to.
* ``ask_value(logger, prompt, help, parse, default, effect, validate=None,
  who=None)``: a value that changes the reduction.
* ``ask_yes_quit(logger, prompt, help, who=None)``: proceed or quit.

Answers common to all of them (whitespace is stripped first):

* **Quit** is exactly ``q`` or ``quit``, in any case.  In ``ask_next`` and
  ``ask_yes_quit``, other text that contains a ``q`` is not a quit: it asks
  again (before Track P it quit, bug B5).  Text without a ``q`` keeps the
  meaning it had before (D12): "next" or "proceed".
* **Help** is ``?`` or ``help``: the long ``help`` text is printed and the
  question is asked again.  It never counts as an answer.

End of input (Ctrl-D, closed or exhausted standard input):

* By default one WARNING naming the prompt is logged and the default answer
  is used.  From then on every prompt in the process returns its default
  without reading input and without warning again (``pause`` proceeds,
  ``ask_next`` returns ``NEXT``, ``ask_value`` returns its ``default``,
  ``ask_yes_quit`` returns ``True``).  ``reset_end_of_input()`` clears this
  state (for tests).
* With the environment variable ``KCWI_STRICT_INPUT=1`` (read at each call)
  the ``EOFError`` is re-raised, as before Track P.  Replay checks use this.

Explanations (``what``, ``help``, the reason an answer was refused) are
printed to standard output and not logged.  Answers that change something
are logged at INFO as ``"<who>: <effect>"``; plain "next" answers at DEBUG.

Input is read with ``builtins.input`` looked up at every call, so a wrapper
installed later (the regression prompt recorder) sees every prompt.  This
module imports neither Matplotlib nor Bokeh nor any kcwidrp plotting module.
"""

import builtins
import logging
import os
import time

NEXT = "next"
QUIT = "quit"

_QUIT_WORDS = ("q", "quit")
_HELP_WORDS = ("?", "help")
_HELP_HINT = " (? for help): "
_NOT_QUIT = "To quit, type q or quit; press <cr> to go on."

_end_of_input = False


def reset_end_of_input():
    """Forget an earlier end of input, so prompts read input again."""
    global _end_of_input
    _end_of_input = False


def _strict():
    return os.environ.get("KCWI_STRICT_INPUT", "").strip() == "1"


def _log(logger):
    return logger if logger is not None else logging.getLogger(__name__)


def _prefix(who):
    return "%s: " % who if who else ""


def _show_help(help):
    # blank line before and after, to set the help off from the prompts
    text = help if help else "No further help for this prompt."
    print("\n" + text + "\n")


def _read(prompt, logger, who, default_text):
    """Read one answer, stripped; ``None`` means "use the default".

    Handles end of input as described in the module docstring:
    ``default_text`` names the default in the one warning.
    """
    global _end_of_input
    strict = _strict()
    if _end_of_input and not strict:
        return None
    try:
        return builtins.input(prompt).strip()
    except EOFError:
        if strict:
            raise
        _end_of_input = True
        _log(logger).warning(
            "%send of input at prompt %r: using the default (%s) here and "
            "for all later prompts" % (_prefix(who), prompt.strip(),
                                       default_text))
        return None


def _with_hint(prompt):
    return prompt if prompt.endswith(": ") else prompt + _HELP_HINT


def pause(config, logger, what=None, help=None, who=None):
    """Pause after a plot: ask at ``plot_level >= 2``, else sleep.

    ``config`` is the instrument config (``self.config.instrument``).  At
    ``plot_level >= 2`` prints ``what`` (if given) and asks
    ``Next? <cr> (? for help): ``; ``?``/``help`` prints ``help`` and asks
    again, any other answer goes on.  Below that, sleeps
    ``config.plot_pause`` seconds and reads no input.

    End of input: by default warns once and goes on (and later prompts go
    on without asking); with ``KCWI_STRICT_INPUT=1`` raises ``EOFError``.
    Returns ``None``.
    """
    if int(config.plot_level) < 2:
        time.sleep(float(config.plot_pause))
        return None
    if what:
        print(what)
    while True:
        answer = _read("Next? <cr>" + _HELP_HINT, logger, who, "next")
        if answer is not None and answer.lower() in _HELP_WORDS:
            _show_help(help)
            continue
        _log(logger).debug("%snext" % _prefix(who))
        return None


def ask_next(logger, what, help, n_items=None, current=None, who=None):
    """Ask whether to go on with a review loop, jump ahead or quit.

    Always asks (the caller keeps its own ``plot_level`` guard).  Prompt:
    ``Next? <cr>, q to quit (? for help): ``, or with ``n_items`` given
    ``Next? <cr>, <int> to jump, q to quit (? for help): ``.

    Returns ``QUIT`` for ``q``/``quit`` (exact, any case), the item number
    for an integer in ``current+1 ... n_items-1`` (``current=None`` is read
    as -1), and ``NEXT`` for an empty answer or other text without a ``q``.
    Text containing a ``q`` that is not a quit word asks again (B5).  With
    ``n_items`` given, an integer outside that range prints the valid range
    and asks again.  ``?``/``help`` prints ``help`` and asks again.

    End of input: by default warns once and returns ``NEXT`` (later prompts
    return ``NEXT`` without asking); with ``KCWI_STRICT_INPUT=1`` raises
    ``EOFError``.
    """
    log = _log(logger)
    if n_items is None:
        prompt = "Next? <cr>, q to quit" + _HELP_HINT
    else:
        prompt = "Next? <cr>, <int> to jump, q to quit" + _HELP_HINT
    first = (current if current is not None else -1) + 1
    if what:
        print(what)
    while True:
        answer = _read(prompt, logger, who, "next")
        if answer is None or answer == "":
            break
        word = answer.lower()
        if word in _HELP_WORDS:
            _show_help(help)
            continue
        if word in _QUIT_WORDS:
            log.info("%suser quit the review%s" % (
                _prefix(who),
                "" if current is None else " at item %d" % current))
            return QUIT
        if "q" in word:
            print(_NOT_QUIT)
            continue
        if n_items is None:
            break
        try:
            item = int(answer)
        except ValueError:
            break
        if first <= item <= n_items - 1:
            log.info("%suser jumped to item %d" % (_prefix(who), item))
            return item
        if first > n_items - 1:
            print("There is no later item to jump to; press <cr> to go on "
                  "or q to quit.")
        else:
            print("Item number must be from %d to %d." % (first, n_items - 1))
    log.debug("%snext" % _prefix(who))
    return NEXT


def ask_value(logger, prompt, help, parse, default, effect, validate=None,
              who=None, what=None):
    """Ask for a value that changes the reduction.

    The prompt shown is ``prompt`` plus `` (? for help): `` unless it already
    ends with ``": "``.  An empty answer returns ``default``.  Otherwise
    ``parse(text)`` makes the value (raising ``ValueError`` on bad input) and
    ``validate(value)``, if given, returns ``None`` or the reason it is
    refused; a refusal prints the reason and asks again.  ``?``/``help``
    prints ``help`` and asks again.  An accepted typed value is logged at
    INFO as ``"<who>: " + effect.format(value=value)``; with ``effect=None``
    nothing is logged, for a caller that logs the effect itself once it has
    been applied.

    End of input: by default warns once and returns ``default`` (later
    prompts return their defaults without asking); with
    ``KCWI_STRICT_INPUT=1`` raises ``EOFError``.
    """
    log = _log(logger)
    shown = _with_hint(prompt)
    if what:
        print(what)
    while True:
        answer = _read(shown, logger, who, "%r" % (default,))
        if answer is None or answer == "":
            log.debug("%sdefault kept: %r" % (_prefix(who), default))
            return default
        if answer.lower() in _HELP_WORDS:
            _show_help(help)
            continue
        try:
            value = parse(answer)
        except ValueError as err:
            print("Not accepted (%s), try again." % (err or "bad format"))
            continue
        reason = validate(value) if validate is not None else None
        if reason:
            print("Not accepted (%s), try again." % reason)
            continue
        if effect is not None:
            log.info(_prefix(who) + effect.format(value=value))
        return value


def ask_yes_quit(logger, prompt, help, who=None, what=None):
    """Ask whether to proceed: ``True`` for ``<cr>``, ``False`` for quit.

    The prompt shown is ``prompt`` plus `` (? for help): `` unless it already
    ends with ``": "``.  ``q``/``quit`` (exact, any case) returns ``False``;
    ``?``/``help`` prints ``help`` and asks again; other text containing a
    ``q`` asks again (B5); any other text proceeds, as before (D12).

    End of input: by default warns once and returns ``True`` (later prompts
    return their defaults without asking); with ``KCWI_STRICT_INPUT=1``
    raises ``EOFError``.
    """
    log = _log(logger)
    shown = _with_hint(prompt)
    if what:
        print(what)
    while True:
        answer = _read(shown, logger, who, "proceed")
        if answer is None or answer == "":
            log.debug("%sproceed" % _prefix(who))
            return True
        word = answer.lower()
        if word in _QUIT_WORDS:
            log.info("%suser chose to quit" % _prefix(who))
            return False
        if word in _HELP_WORDS:
            _show_help(help)
        elif "q" in word:
            print(_NOT_QUIT)
        else:
            log.debug("%sproceed" % _prefix(who))
            return True
