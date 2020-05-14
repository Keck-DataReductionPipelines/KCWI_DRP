"""
Created on Aug 8, 2019

@author: lrizzi
"""

# from bokeh.io import curdoc
from bokeh.plotting import output_file, save
from bokeh.plotting.figure import Figure
from bokeh.models import Column
import psutil
# import subprocess
import os
# import time


def bokeh_plot(plot, session):

    figure = session.document.select_one(selector=dict(type=Figure))
    layout = session.document.select_one(selector=dict(type=Column))

    layout.children.remove(figure)
    layout.children.insert(0, plot)
    session.push()


def bokeh_add_plot(glyph, session):

    figure = session.document.select_one(selector=dict(type=Figure))
    new_figure = figure.glyph
    c.childen.remove(p)
    c.children.insert(0, new_p)
    session.push()
    # session.show(c)


def check_bokeh_server():
    """
    Check if there is any running process that contains
    the given name processName.
    """
    # Iterate over the all the running process
    for proc in psutil.process_iter(attrs=["name", "cmdline"]):
        try:
            # Check if process name contains the given name string.
            for command in proc.cmdline():
                if 'bokeh' in command:
                    return True
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass
    return False


def bokeh_save(plot):
    cwd = os.getcwd()
    filename = os.path.join(cwd, 'plots', 'plot.html')
    # try:
    #    os.remove(filename)
    # except OSError:
    #    pass
    output_file(filename)
    save(plot)

    # time.sleep(1)
    # subprocess.Popen("open %s" % filename, shell=True)
