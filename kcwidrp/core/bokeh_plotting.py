"""
Created on Aug 8, 2019

@author: lrizzi
"""

from bokeh.io import curdoc
from bokeh.plotting.figure import Figure
from bokeh.models import Column


def bokeh_plot(plot):

    p = curdoc().select_one(selector=dict(type=Figure))
    c = curdoc().select_one(selector=dict(type=Column))
    c.children.remove(p)
    c.children.insert(0, plot)


def bokeh_add_plot(glyph):

    p = curdoc().select_one(selector=dict(type=Figure))
    new_p = p.glyph
    c.childen.remove(p)
    c.children.insert(0, new_p)
