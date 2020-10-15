from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import read_table
from kcwidrp.core.bokeh_plotting import bokeh_plot

import numpy as np
import os
from skimage import transform as tf
from bokeh.plotting import figure


class ExtractArcs(BasePrimitive):
    """Use derived traces to extract arc spectra along bars"""

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
        contbars_in_proctable = self.context.proctab.n_proctab(
            frame=self.action.args.ccddata, target_type='CONTBARS',
            nearest=True)
        self.logger.info("%d continuum bars frames found" %
                         len(contbars_in_proctable))
        if len(contbars_in_proctable) > 0:
            self.action.args.original_filename = contbars_in_proctable[
                'OFNAME'][0].split('.')[0] + "_trace.fits"
            return True
        else:
            self.action.args.original_filename = None
            return False

    def _perform(self):
        do_plot = (self.config.instrument.plot_level >= 3)
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
        if hasattr(self.context, 'trace'):
            trace = self.context.trace
        else:
            trace = read_table(
                input_dir=os.path.join(os.path.dirname(self.action.args.name),
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
        self.action.args.arc_image = self.action.args.ccddata.header['OFNAME']

        self.action.args.source_control_points = trace['src']
        self.action.args.destination_control_points = trace['dst']
        self.action.args.bar_id = trace['barid']
        self.action.args.slice_id = trace['slid']

        self.logger.info("Fitting spatial control points")
        transformation = tf.estimate_transform(
            'polynomial', self.action.args.source_control_points,
            self.action.args.destination_control_points, order=3)

        self.logger.info("Transforming arc image")
        warped_image = tf.warp(self.action.args.ccddata.data, transformation)
        # Write warped arcs if requested
        if self.config.instrument.saveintims:
            from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer
            # write out warped image
            self.action.args.ccddata.data = warped_image
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name,
                             output_dir=self.config.instrument.output_directory,
                             suffix="warped")
            self.logger.info("Transformed arcs produced")
        # extract arcs
        self.logger.info("Extracting arcs")
        arcs = []
        # sectors for bkgnd subtraction
        sectors = 16
        for xyi, xy in enumerate(self.action.args.source_control_points):
            if xy[1] == middle_row:
                xi = int(xy[0]+0.5)
                arc = np.median(
                    warped_image[:, (xi - window):(xi + window + 1)], axis=1)
                # divide spectrum into sectors
                div = int((len(arc)-100) / sectors)
                # get minimum for each sector
                xv = []
                yv = []
                for i in range(sectors):
                    mi = np.nanargmin(arc[50+i*div:50+(i+1)*div])
                    mn = np.nanmin(arc[50+i*div:50+(i+1)*div])
                    xv.append(mi+50+i*div)
                    yv.append(mn)
                # fit minima to model background
                res = np.polyfit(xv, yv, 3)
                xp = np.arange(len(arc))
                bkg = np.polyval(res, xp)   # resulting model
                # plot if requested
                if do_plot:
                    p = figure(title=self.action.args.plotlabel + "ARC # %d" %
                               len(arcs),
                               x_axis_label="Y CCD Pixel",
                               y_axis_label="Flux",
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
        else:
            self.logger.error("Did not extract %d arcs, extracted %d" %
                              (self.config.instrument.NBARS, len(arcs)))

        log_string = ExtractArcs.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
# END: class ExtractArcs()
