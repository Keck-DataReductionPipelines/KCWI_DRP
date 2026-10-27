"""bless with a symlinked reference writes at the link's target.

Written from the specification: ``reference/<entry>`` may be a directory or a
symlink (e.g. to external storage); for a symlink the new reference is built
and swapped in at the target and the link itself is left alone.
"""
import io
import os
import stat

import pytest

from kcwidrp.tests.regression import bless as blessmod
from kcwidrp.tests.regression import manifest as mf
from kcwidrp.tests.regression.report import SessionResults
from kcwidrp.tests.regression.session import EntryRun

EID = "ds9_blue_L1"


@pytest.fixture
def world(make_world, tmp_path):
    """A kept run of an entry that has no reference yet."""
    w = make_world()
    session_dir = tmp_path / "out" / "20260101-000000"
    results = SessionResults(session_dir)
    entry = mf.get_entry(w.mfst, EID)
    run = EntryRun(w.data, w.mfst, entry, results, session_dir,
                   command=w.command, keep_runs=True).start()
    for t in range(8):
        run.tier(t)
    run.finish()
    w.session = session_dir
    return w


def do_bless(w, out=None, reason="environment: x"):
    return blessmod.bless(w.data, EID, reason, run=w.session,
                          assume_yes=True, out=out)


def make_target(parent, name="ext_ref"):
    """An existing external reference directory holding an old file."""
    target = parent / name
    (target / "redux").mkdir(parents=True)
    (target / "redux" / "old.txt").write_text("old reference")
    return target


def link_reference(w, target):
    (w.data / "reference").mkdir(exist_ok=True)
    link = w.data / "reference" / EID
    link.symlink_to(target)
    return link


@pytest.fixture
def readonly_parent(tmp_path):
    parent = tmp_path / "ro"
    parent.mkdir()
    yield parent
    parent.chmod(stat.S_IRWXU)


def test_bless_writes_the_new_reference_at_the_symlink_target(world, tmp_path):
    w = world
    target = make_target(tmp_path)
    link = link_reference(w, target)
    link_to = os.readlink(link)
    do_bless(w)
    # the link is unchanged
    assert link.is_symlink()
    assert os.readlink(link) == link_to
    # the target holds the new run's outputs, not the old reference
    assert (target / "redux" / "ds_icubes.fits").is_file()  # source: redux file written by the fake reduction (test_bless_requires_pass)
    assert (target / "kcwib.proc").is_file()  # source: procfile name (test_bless_records_its_reason...)
    assert not (target / "redux" / "old.txt").exists()
    # no leftovers next to the target, nothing but the link in reference/
    assert not (tmp_path / ".ext_ref.new").exists()
    assert not (tmp_path / ".ext_ref.old").exists()
    assert os.listdir(w.data / "reference") == [EID]
    hist = mf.get_entry(mf.load_manifest(w.data), EID)["bless_history"]
    assert len(hist) == 1


def test_bless_symlink_does_not_copy_outputs_into_the_data_dir(world,
                                                               tmp_path):
    w = world
    target = make_target(tmp_path)
    link_reference(w, target)
    do_bless(w)
    # only the link itself lives under reference/; no real directory
    entries = list((w.data / "reference").iterdir())
    assert [p.name for p in entries] == [EID]
    assert entries[0].is_symlink()
    assert not (w.data / "reference" / ".ds9_blue_L1.new").exists()


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores directory modes")
def test_unwritable_target_parent_is_refused_untouched(world, tmp_path,
                                                       readonly_parent):
    w = world
    target = make_target(readonly_parent)
    link = link_reference(w, target)
    readonly_parent.chmod(stat.S_IRUSR | stat.S_IXUSR)
    with pytest.raises(blessmod.BlessError) as exc:
        do_bless(w)
    msg = str(exc.value)
    assert str(target) in msg
    assert str(link) in msg
    assert (target / "redux" / "old.txt").read_text() == "old reference"
    assert sorted(p.name for p in target.iterdir()) == ["redux"]
    assert sorted(p.name for p in readonly_parent.iterdir()) == ["ext_ref"]
    assert link.is_symlink() and os.readlink(link) == str(target)
    assert mf.get_entry(mf.load_manifest(w.data), EID)["bless_history"] == []


def test_missing_target_parent_is_refused(world, tmp_path):
    w = world
    target = tmp_path / "no_such_dir" / "ref"
    link = link_reference(w, target)
    with pytest.raises(blessmod.BlessError) as exc:
        do_bless(w)
    assert str(target) in str(exc.value)
    assert str(link) in str(exc.value)
    assert not (tmp_path / "no_such_dir").exists()
    assert link.is_symlink() and os.readlink(link) == str(target)
    assert mf.get_entry(mf.load_manifest(w.data), EID)["bless_history"] == []


def test_plan_names_the_symlink_target(world, tmp_path):
    w = world
    target = make_target(tmp_path)
    link_reference(w, target)
    buf = io.StringIO()
    do_bless(w, out=buf)
    text = buf.getvalue()
    assert "symlink to" in text
    assert str(target) in text


def test_plan_of_a_plain_directory_has_no_symlink_line(world):
    w = world
    buf = io.StringIO()
    do_bless(w, out=buf)
    assert "symlink to" not in buf.getvalue()


def test_reference_destination_of_a_directory_is_the_directory(tmp_path):
    ref = tmp_path / "ref"
    ref.mkdir()
    assert blessmod.reference_destination(ref) == ref


def test_reference_destination_of_a_symlink_is_its_resolved_target(tmp_path):
    target = make_target(tmp_path)
    link = tmp_path / "link"
    link.symlink_to(target)
    assert blessmod.reference_destination(link) == target.resolve()


def test_reference_destination_of_a_symlink_with_missing_parent_raises(
        tmp_path):
    link = tmp_path / "link"
    link.symlink_to(tmp_path / "nope" / "ref")
    with pytest.raises(blessmod.BlessError):
        blessmod.reference_destination(link)


def test_plain_directory_reference_is_still_replaced(world):
    w = world
    ref = w.data / "reference" / EID
    (ref / "redux").mkdir(parents=True)
    (ref / "redux" / "stale.txt").write_text("old")
    do_bless(w)
    assert not ref.is_symlink()
    assert not (ref / "redux" / "stale.txt").exists()
    assert (ref / "redux" / "ds_icubes.fits").is_file()
    assert sorted(os.listdir(w.data / "reference")) == [EID]
    assert len(mf.get_entry(mf.load_manifest(w.data),
                            EID)["bless_history"]) == 1
