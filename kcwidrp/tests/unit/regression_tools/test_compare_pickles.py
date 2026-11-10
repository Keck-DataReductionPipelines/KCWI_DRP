"""Pickle comparison (diff_objects, compare_pickle_pair, compare_outputs).

Written from the specification only. The motivating case: a geometry pickle
whose file-path strings were read via another directory must compare equal,
while any number, array value, transform parameter or file name change must
not.
"""
import pickle

import numpy as np
from skimage.transform import PolynomialTransform

from kcwidrp.tests.regression import compare
from kcwidrp.tests.regression.compare import (OutputComparison, compare_outputs,
                                              compare_pickle_pair,
                                              diff_objects)

PROC = "kcwib.proc"


def tform(scale=1.0):
    params = np.array([[0.0, scale, 0.0, 0.0, 0.0, 0.0],
                       [0.0, 0.0, scale, 0.0, 0.0, 0.0]])
    return PolynomialTransform(params)


def geom(arcfl="/data/a/kb_arc.fits", cbarsfl="/data/a/kb_cbars.fits",
         xl=3.5, scale=1.0):
    return {"arcfl": arcfl, "cbarsfl": cbarsfl, "xl0": xl, "nasmask": 4,
            "tform": [tform(scale), tform(scale)],
            "invtf": [tform(scale)]}


def test_equal_objects_have_no_differences():
    assert diff_objects(geom(), geom()) == []     # source: analytic (a == a)


def test_path_in_another_directory_with_same_file_name_is_equal():
    other = geom(arcfl="/elsewhere/b/kb_arc.fits",
                 cbarsfl="/other/c/kb_cbars.fits")
    assert diff_objects(geom(), other) == []      # source: spec, path rule


def test_different_file_name_is_a_difference_naming_the_key():
    diffs = diff_objects(geom(), geom(arcfl="/data/a/kb_other.fits"))
    assert len(diffs) == 1                        # source: spec, one per diff
    assert "arcfl" in diffs[0]


def test_float_change_of_1e_minus_9_is_a_difference():
    diffs = diff_objects({"xl0": 3.5}, {"xl0": 3.5 + 1e-9})
    assert len(diffs) == 1 and "xl0" in diffs[0]  # source: spec, exact


def test_nan_equals_nan():
    nan = float("nan")
    assert diff_objects({"v": nan}, {"v": float("nan")}) == []
    # source: spec, NaN rule for floats


def test_nan_in_float_and_complex_arrays_is_equal():
    f = np.array([1.0, np.nan])
    c = np.array([1 + 1j, complex(np.nan, 0)])
    assert diff_objects(f, f.copy()) == []        # source: spec
    assert diff_objects(c, c.copy()) == []        # source: spec


def test_array_value_change_is_a_difference():
    a = np.arange(6, dtype=float).reshape(2, 3)
    b = a.copy()
    b[1, 2] += 1e-12
    assert len(diff_objects({"arr": a}, {"arr": b})) == 1  # source: spec
    assert "arr" in diff_objects({"arr": a}, {"arr": b})[0]


def test_array_shape_change_is_a_difference():
    a = np.zeros((2, 3))
    assert len(diff_objects(a, np.zeros((3, 2)))) >= 1     # source: spec


def test_array_dtype_change_is_a_difference():
    a = np.zeros(4, dtype="f8")
    assert len(diff_objects(a, a.astype("f4"))) >= 1       # source: spec


def test_transform_parameter_change_is_a_difference():
    diffs = diff_objects(geom(), geom(scale=1.0 + 1e-9))
    assert len(diffs) >= 1                        # source: spec, __dict__ rule
    assert diff_objects(tform(), tform()) == []   # source: spec


def test_type_change_int_vs_float_is_a_difference():
    diffs = diff_objects({"n": 4}, {"n": 4.0})
    assert len(diffs) == 1 and "n" in diffs[0]    # source: spec, type rule


def test_missing_key_is_a_difference_on_either_side():
    d1 = diff_objects({"a": 1, "b": 2}, {"a": 1})
    d2 = diff_objects({"a": 1}, {"a": 1, "b": 2})
    assert len(d1) == 1 and "b" in d1[0]          # source: spec
    assert len(d2) == 1 and "b" in d2[0]          # source: spec


def test_list_and_tuple_length_and_elements():
    assert len(diff_objects([1, 2, 3], [1, 2])) >= 1       # source: spec
    assert len(diff_objects((1, 2), (1, 2, 3))) >= 1       # source: spec
    assert len(diff_objects([1, 2, 3], [1, 9, 3])) == 1    # source: spec
    assert diff_objects([1, (2, 3)], [1, (2, 3)]) == []


def test_limit_caps_the_number_of_differences():
    a = {k: 0 for k in range(50)}
    b = {k: 1 for k in range(50)}
    assert len(diff_objects(a, b)) == 20          # source: spec, default limit
    assert len(diff_objects(a, b, limit=3)) == 3  # source: spec
    assert len(diff_objects(a, b, limit=100)) == 50


def test_where_prefix_appears_in_messages():
    diffs = diff_objects({"k": 1}, {"k": 2}, where="top")
    assert len(diffs) == 1 and "top" in diffs[0] and "k" in diffs[0]


def dump(path, obj):
    with open(path, "wb") as f:
        pickle.dump(obj, f)


def test_pickle_pair_equal_objects_add_no_findings(tmp_path):
    dump(tmp_path / "a.pkl", geom())
    dump(tmp_path / "b.pkl", geom(arcfl="/elsewhere/kb_arc.fits"))
    out = OutputComparison()
    compare_pickle_pair(tmp_path / "a.pkl", tmp_path / "b.pkl",
                        "redux/x_geom.pkl", out)
    assert out.pixels == [] and out.pixel_messages == []


def test_pickle_pair_differences_start_with_the_relative_name(tmp_path):
    dump(tmp_path / "a.pkl", geom())
    dump(tmp_path / "b.pkl", geom(xl=3.6))
    out = OutputComparison()
    compare_pickle_pair(tmp_path / "a.pkl", tmp_path / "b.pkl",
                        "redux/x_geom.pkl", out)
    assert len(out.pixels) == 1                   # source: spec, one per diff
    assert len(out.pixel_messages) == 1
    assert out.pixel_messages[0].startswith("redux/x_geom.pkl")
    assert "xl0" in out.pixel_messages[0]


def test_unpicklable_files_fall_back_to_bytes(tmp_path):
    (tmp_path / "a.pkl").write_bytes(b"not a pickle")
    (tmp_path / "same.pkl").write_bytes(b"not a pickle")
    (tmp_path / "diff.pkl").write_bytes(b"not a pickle either")
    out = OutputComparison()
    compare_pickle_pair(tmp_path / "a.pkl", tmp_path / "same.pkl", "r.pkl", out)
    assert out.pixels == [] and out.pixel_messages == []
    compare_pickle_pair(tmp_path / "a.pkl", tmp_path / "diff.pkl", "r.pkl", out)
    assert len(out.pixels) == 1 and len(out.pixel_messages) == 1
    assert out.pixel_messages[0].startswith("r.pkl")


def make_tree(root, obj):
    (root / "redux").mkdir(parents=True)
    (root / "plots").mkdir()
    dump(root / "redux" / "x_geom.pkl", obj)
    (root / PROC).write_text(
        "| FRAMENO | TYPE | filename |\n"
        "|       1 | BIAS | KB.1.fits |\n")


def test_compare_outputs_ignores_pickle_path_directory(tmp_path):
    make_tree(tmp_path / "run", geom(arcfl="/run/dir/kb_arc.fits"))
    make_tree(tmp_path / "ref", geom(arcfl="/ref/dir/kb_arc.fits"))
    out = compare_outputs(tmp_path / "run", tmp_path / "ref", PROC)
    assert out.pixels == [] and out.pixel_messages == []
    assert out.structure == []


def test_compare_outputs_flags_pickle_number_difference(tmp_path):
    make_tree(tmp_path / "run", geom(xl=3.5 + 1e-9))
    make_tree(tmp_path / "ref", geom())
    out = compare_outputs(tmp_path / "run", tmp_path / "ref", PROC)
    assert len(out.pixel_messages) == 1           # source: spec
    assert out.pixel_messages[0].startswith("redux/x_geom.pkl")
