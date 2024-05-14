from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import read_table, kcwi_fits_reader
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.primitives.kcwi_file_primitives import strip_fname, plotlabel

import numpy as np
from numpy.polynomial import polynomial as poly
from scipy.signal import savgol_filter
import os
from skimage import transform as tf
from bokeh.plotting import figure


class ExtractArcs(BasePrimitive):
    """
    Use derived traces to extract arc spectra along continuum bars.

    Reads in traces from continuum bars and then uses them to extract spectra
    along each trace.  Also performs a background subtraction for each extracted
    spectrum.  Adds a HISTORY record to the FITS header.

    Uses the following configuration parameter:

        * saveintims: if set to ``True`` write out a warped version of arc image in \*_warped.fits.  Defaults to ``False``.

    Adds extracted arcs to returned arguments.

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        self.action.args.reference_bar_separation = None
        self.action.args.contbar_image_number = None
        self.action.args.contbar_image = None
        self.action.args.arc_number = None
        self.action.args.arc_image = None
        self.action.args.source_control_points = None
        self.action.args.destination_control_points = None
        self.action.args.bar_id = None
        self.action.args.slice_id = None

    def _pre_condition(self):
        marcs_in_proctable = self.context.proctab.search_proctab(
            frame=self.action.args.ccddata, target_type='MARC',
            nearest=True
        )
        self.logger.info("%d master arc frames found" %
                         len(marcs_in_proctable))
        contbars_in_proctable = self.context.proctab.search_proctab(
            frame=self.action.args.ccddata, target_type='MCBARS',
            nearest=True)
        self.logger.info("%d master continuum bars frames found" %
                         len(contbars_in_proctable))
        if len(contbars_in_proctable) > 0 and len(marcs_in_proctable) > 0:
            self.action.args.original_filename = strip_fname(
                contbars_in_proctable['filename'][0]) + "_trace.fits"
            return True
        else:
            self.action.args.original_filename = None
            return False

    def _perform(self):
        do_plot = (self.config.instrument.plot_level >= 3)
        plab = plotlabel(self.action.args)
        self.logger.info("Extracting arc spectra")
        # Double check
        if not self.action.args.original_filename:
            self.logger.error("No traces found")
            return self.action.args
        # All is ready
        original_filename = self.action.args.original_filename
        self.logger.info("Trace table found: %s" % original_filename)
        # trace = read_table(tab=tab, indir='redux', suffix='trace')
        # Find  and read control points from continuum bars
        trace = read_table(
            input_dir=os.path.join(self.config.instrument.cwd,
                                   self.config.instrument.output_directory),
            file_name=original_filename)
        self.context.trace = {}
        for key in trace.meta.keys():
            self.context.trace[key] = trace.meta[key]
        middle_row = self.context.trace['MIDROW']
        window = self.context.trace['WINDOW']
        self.action.args.reference_bar_separation = self.context.trace[
            'REFDELX']
        self.action.args.contbar_image_number = self.context.trace['CBARSNO']
        self.action.args.contbar_image = self.context.trace['CBARSFL']
        self.action.args.arc_number = self.action.args.ccddata.header['FRAMENO']
        self.action.args.arc_image = self.action.args.name
        # self.action.args.arc_image = self.action.args.ccddata.header['OFNAME']

        self.action.args.source_control_points = trace['src']
        self.action.args.destination_control_points = trace['dst']
        self.action.args.bar_id = trace['barid']
        self.action.args.slice_id = trace['slid']

        self.logger.info("Fitting spatial control points")

        # This is strictly for AIT data!
        if self.config.instrument.NBARS < 60:
            # get trimmed bars image so geometry is consistent
            barsfn = self.action.args.contbar_image.split(
                '.fits')[0] + '_int.fits'
            barim, bartab = kcwi_fits_reader(
                os.path.join(self.config.instrument.output_directory, barsfn))
            arcs = []
            bars = []
            for ib in range(self.config.instrument.NBARS):
                bind = np.where(trace['barid'] == ib)
                xp = trace['dst'][bind][:, 0]
                yp = trace['dst'][bind][:, 1]
                sind = np.argsort(yp)
                yps = yp[sind]
                xps = xp[sind]
                c, stats = poly.polyfit(yps, xps, deg=3, full=True)
                arc = []
                bar = []
                bkg = []
                for iy in range(self.action.args.ccddata.header['NAXIS2']):
                    xv = int(poly.polyval(iy, c))
                    yy = self.action.args.ccddata.data[
                         iy, (xv - window):(xv + window + 1)]
                    byy = barim.data[iy, (xv - window):(xv + window + 1)]
                    nyy = len(yy)
                    if nyy > 0:
                        bkg.append(np.nanmin(yy) * float(nyy))
                        arc.append(np.nansum(yy))
                        bar.append(np.nansum(byy))
                    else:
                        bkg.append(0.)
                        arc.append(0.)
                        bar.append(0.)
                arc = np.asarray(arc, dtype=float)
                bar = np.asarray(bar, dtype=float)
                bkg = np.asarray(bkg, dtype=float)
                bkgf = savgol_filter(bkg, 51, 5)
                if do_plot:
                    xp = np.arange(len(arc))
                    p = figure(title=plab + "ARC # %d" % len(arcs),
                               x_axis_label="Y CCD Pixel",
                               y_axis_label="Flux",
                               plot_width=self.config.instrument.plot_width,
                               plot_height=self.config.instrument.plot_height)
                    p.line(xp, arc, legend_label='Arc', color='blue')
                    p.line(xp, bkg, legend_label='Bkg', color='red')
                    p.line(xp, bkgf, legend_label='BkgFit', color='green')
                    bokeh_plot(p, self.context.bokeh_session)
                    q = input("Next? <cr>, q to quit: ")
                    if 'Q' in q.upper():
                        do_plot = False
                arc -= bkgf
                arcs.append(arc)
                bars.append(bar)
        else:
            # fit transform
            # NOTE: we do not need an asymmetric polynomial for this
            # global fit (only for per slice fitting: see SolveGeom.py)
            tform = tf.estimate_transform(
                'polynomial', self.action.args.source_control_points,
                self.action.args.destination_control_points, order=3)

            self.logger.info("Transforming arc image")
            warped_image = tf.warp(self.action.args.ccddata.data, tform)
            # Write warped arcs if requested
            if self.config.instrument.saveintims:
                from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer
                # write out warped image
                self.action.args.ccddata.data = warped_image
                kcwi_fits_writer(
                    self.action.args.ccddata, table=self.action.args.table,
                    output_file=self.action.args.name,
                    output_dir=self.config.instrument.output_directory,
                    suffix="warped")
                self.logger.info("Transformed arcs produced")
            # extract arcs
            self.logger.info("Extracting arcs")
            arcs = []
            bars = []
            # sectors for bkgnd subtraction
            sectors = 16
            for xyi, xy in enumerate(self.action.args.source_control_points):
                if xy[1] == middle_row:
                    xi = int(xy[0]+0.5)
                    arc = np.nanmedian(
                        warped_image[:, (xi - window):(xi + window + 1)],
                        axis=1)
                    # divide spectrum into sectors
                    div = int((len(arc)-100) / sectors)
                    # get minimum for each sector
                    xv = []
                    yv = []
                    for i in range(sectors):
                        try:
                            mi = np.nanargmin(arc[50+i*div:50+(i+1)*div])
                            mn = np.nanmin(arc[50+i*div:50+(i+1)*div])
                            xv.append(mi+50+i*div)
                            yv.append(mn)
                        except ValueError:
                            self.logger.warning("Bad sector %d" % i)
                    # fit minima to model background
                    res = np.polyfit(xv, yv, 3)
                    xp = np.arange(len(arc))
                    bkg = np.polyval(res, xp)   # resulting model
                    # plot if requested
                    if do_plot:
                        p = figure(
                            title=plab + "ARC # %d" % len(arcs),
                            x_axis_label="Y CCD Pixel", y_axis_label="Flux",
                            plot_width=self.config.instrument.plot_width,
                            plot_height=self.config.instrument.plot_height)
                        p.line(xp, arc, legend_label='Arc', color='blue')
                        p.line(xp, bkg, legend_label='Bkg', color='red')
                        bokeh_plot(p, self.context.bokeh_session)
                        q = input("Next? <cr>, q to quit: ")
                        if 'Q' in q.upper():
                            do_plot = False
                    # subtract model background
                    arc -= bkg
                    # add to arcs list
                    arcs.append(arc)
        # Did we get the correct number of arcs?
        if len(arcs) == self.config.instrument.NBARS:
            self.logger.info("Extracted %d arcs" % len(arcs))
            self.context.arcs = arcs
            if self.config.instrument.NBARS < 60:
                self.context.bars = bars
        else:
            self.logger.error("Did not extract %d arcs, extracted %d" %
                              (self.config.instrument.NBARS, len(arcs)))

        log_string = ExtractArcs.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
# END: class ExtractArcs()
