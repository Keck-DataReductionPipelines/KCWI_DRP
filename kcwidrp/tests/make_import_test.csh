# to run this command:
# csh make_import_test.csh


if ( -f  test_import.py ) then
  rm test_import.py

echo "import pytest" > test_import.py

echo "" | awk '{print "\n\n"}' >> test_import.py

ls -1 ../primitives/*.py | grep -v init | sed s/.py//g | cut -d / -f 3 | awk '{print "def test_import_"$1"():\n    import kcwidrp.primitives."$1"\n\n"'} >> test_import.py
