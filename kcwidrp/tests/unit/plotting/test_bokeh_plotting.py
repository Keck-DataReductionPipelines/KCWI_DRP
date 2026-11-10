import multiprocessing
from kcwidrp.core.bokeh_plotting import check_running_process
import time


def test_check_running_process():

    # check that a python process is running
    assert check_running_process('python') is True
    # check that a non-running process is not running
    assert check_running_process('abcdef') is False


