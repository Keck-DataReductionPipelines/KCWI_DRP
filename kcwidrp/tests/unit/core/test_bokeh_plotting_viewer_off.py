"""Viewer-off crashes B9 and B10 (plan T.11): bokeh_plotting helpers.

No Bokeh server, browser, network or display is used.
"""


def test_bokeh_clear_none_returns_without_raising():
    """B10: bokeh_clear(None) returns at once, as bokeh_plot(None) does."""
    from kcwidrp.core.bokeh_plotting import bokeh_clear

    result = bokeh_clear(None)

    # source: specification (plan T.11 / §13 B9, B10)
    assert result is None


def test_init_plot_session_sets_none_when_attribute_missing():
    """B9: a context without bokeh_session gets bokeh_session = None."""
    from kcwidrp.core.bokeh_plotting import init_plot_session

    class Context:
        pass

    context = Context()
    assert not hasattr(context, "bokeh_session")

    init_plot_session(context)

    # source: specification (plan T.11 / §13 B9, B10)
    assert hasattr(context, "bokeh_session")
    # source: specification (plan T.11 / §13 B9, B10)
    assert context.bokeh_session is None


def test_init_plot_session_keeps_existing_session():
    """B9: an existing bokeh_session is left as the same object."""
    from kcwidrp.core.bokeh_plotting import init_plot_session

    class Context:
        pass

    sentinel = object()
    context = Context()
    context.bokeh_session = sentinel

    init_plot_session(context)

    # source: specification (plan T.11 / §13 B9, B10)
    assert context.bokeh_session is sentinel
