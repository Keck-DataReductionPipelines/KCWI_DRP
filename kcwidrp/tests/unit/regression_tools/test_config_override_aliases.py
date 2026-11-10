"""apply_config_overrides treats enable_bokeh / enable_plot_viewer as aliases.

PR1 renamed the viewer switch; a manifest drafted on either side of the rename
must run on code from either side.
"""
import pytest

from kcwidrp.tests.regression import runner

NEW_ONLY = "# header\nenable_plot_viewer = True\nplot_level = 1\n"
OLD_ONLY = "# header\nenable_bokeh = True\nplot_level = 1\n"


def test_old_name_override_rewrites_new_name_line():
    out = runner.apply_config_overrides(NEW_ONLY, {"enable_bokeh": False})
    # source: analytic - only the enable_plot_viewer line changes
    assert out == "# header\nenable_plot_viewer = False\nplot_level = 1\n"
    assert "enable_bokeh" not in out    # source: analytic - no old-name line added


def test_new_name_override_rewrites_old_name_line():
    out = runner.apply_config_overrides(OLD_ONLY, {"enable_plot_viewer": False})
    # source: analytic - only the enable_bokeh line changes
    assert out == "# header\nenable_bokeh = False\nplot_level = 1\n"
    assert "enable_plot_viewer" not in out    # source: analytic - none added


def test_name_present_itself_is_rewritten():
    out = runner.apply_config_overrides(OLD_ONLY, {"enable_bokeh": False})
    assert out == "# header\nenable_bokeh = False\nplot_level = 1\n"  # source: analytic
    out = runner.apply_config_overrides(NEW_ONLY,
                                        {"enable_plot_viewer": False})
    assert out == "# header\nenable_plot_viewer = False\nplot_level = 1\n"  # source: analytic


def test_packaged_config_with_manifest_overrides():
    text = runner.packaged_config_text()
    out = runner.apply_config_overrides(
        text, {"enable_bokeh": False, "plot_pause": 0, "plot_level": 3})
    lines = {ln.split("=")[0].strip(): ln.split("=")[1].strip()
             for ln in out.splitlines()
             if "=" in ln and not ln.lstrip().startswith("#")}
    # source: analytic - the values the overrides ask for, rendered with str()
    assert lines["enable_plot_viewer"] == "False"
    assert lines["plot_pause"] == "0"
    assert lines["plot_level"] == "3"
    assert "enable_bokeh" not in lines    # source: analytic - no line added


def test_neither_name_present_raises():
    with pytest.raises(runner.RunError):    # source: analytic - key in neither form
        runner.apply_config_overrides("plot_level = 1\n",
                                      {"enable_bokeh": False})
    with pytest.raises(runner.RunError):    # source: analytic - key in neither form
        runner.apply_config_overrides("plot_level = 1\n",
                                      {"enable_plot_viewer": False})


def test_unrelated_missing_key_raises():
    with pytest.raises(runner.RunError):    # source: analytic - key absent
        runner.apply_config_overrides(NEW_ONLY, {"no_such_key": 1})


def test_only_one_line_changes_and_indent_kept():
    text = "a = 1\n  enable_plot_viewer = True\nb = 2\n"
    out = runner.apply_config_overrides(text, {"enable_bokeh": False})
    # source: analytic - the indented line changes, indentation preserved
    assert out == "a = 1\n  enable_plot_viewer = False\nb = 2\n"
