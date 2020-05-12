from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import read_table

import numpy as np
import os
from skimage import transform as tf


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

    def _perform(self):
        self.logger.info("Extracting arc spectra")
        contbars_in_proctable = self.context.proctab.n_proctab(
            frame=self.action.args.ccddata, target_type='CONTBARS',
            nearest=True)
        self.logger.info("%d continuum bars frames found" %
                         len(contbars_in_proctable))
        # ofname = tab['OFNAME'][0]
        original_filename = contbars_in_proctable['OFNAME'][0].split('.')[0] + \
            "_trace.fits"
        self.logger.info("Reading Trace Table: %s" % original_filename)
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
        for xyi, xy in enumerate(self.action.args.source_control_points):
            if xy[1] == middle_row:
                xi = int(xy[0]+0.5)
                arc = np.median(
                    warped_image[:, (xi - window):(xi + window + 1)], axis=1)
                arc = arc - np.nanmin(arc[100:-100])    # avoid ends
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
