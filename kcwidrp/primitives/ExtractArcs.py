from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from kcwidrp.primitives.kcwi_file_primitives import read_table

from skimage import transform as tf
import numpy as np
import os


class ExtractArcs(BasePrimitive):
    """Use derived traces to extract arc spectra along bars"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        self.action.args.refdelx = None
        self.action.args.cbarsno = None
        self.action.args.cbarsfl = None
        self.action.args.arcno = None
        self.action.args.arcfl = None
        self.action.args.src = None
        self.action.args.dst = None
        self.action.args.barid = None
        self.action.args.slid = None

    def _perform(self):
        self.logger.info("Extracting arc spectra")
        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type='CONTBARS',
                                             nearest=True)
        self.logger.info("%d continuum bars frames found" % len(tab))
        # ofname = tab['OFNAME'][0]
        ofname = tab['OFNAME'][0].split('.')[0] + "_trace.fits"
        print("*************** READING TABLE: %s" % ofname)
        # trace = read_table(tab=tab, indir='redux', suffix='trace')
        # Find  and read control points from continuum bars
        if hasattr(self.context, 'trace'):
            trace = self.context.trace
        else:
            trace = read_table(input_dir=os.path.dirname(self.action.args.name),
                               file_name=ofname)
            self.context.trace = {}
            for key in trace.meta.keys():
                self.context.trace[key] = trace.meta[key]
        midrow = self.context.trace['MIDROW']
        win = self.context.trace['WINDOW']
        self.action.args.refdelx = self.context.trace['REFDELX']
        self.action.args.cbarsno = self.context.trace['CBARSNO']
        self.action.args.cbarsfl = self.context.trace['CBARSFL']
        self.action.args.arcno = self.action.args.ccddata.header['FRAMENO']
        self.action.args.arcfl = self.action.args.ccddata.header['OFNAME']

        self.action.args.src = trace['src']  # source control points
        self.action.args.dst = trace['dst']  # destination control points
        self.action.args.barid = trace['barid']
        self.action.args.slid = trace['slid']

        self.logger.info("Fitting spatial control points")
        tform = tf.estimate_transform('polynomial',
                                      self.action.args.src,
                                      self.action.args.dst, order=3)

        self.logger.info("Transforming arc image")
        warped = tf.warp(self.action.args.ccddata.data, tform)
        # Write warped arcs if requested
        if self.config.instrument.saveintims:
            # write out warped image
            self.action.args.ccddata.data = warped
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name,
                             suffix="warped")
            self.logger.info("Transformed arcs produced")
        # extract arcs
        self.logger.info("Extracting arcs")
        arcs = []
        for xyi, xy in enumerate(self.action.args.src):
            if xy[1] == midrow:
                xi = int(xy[0]+0.5)
                arc = np.median(
                    warped[:, (xi - win):(xi + win + 1)], axis=1)
                arc = arc - np.nanmin(arc[100:-100])    # avoid ends
                arcs.append(arc)
        # Did we get the correct number of arcs?
        if len(arcs) == self.config.instrument.NBARS:
            self.logger.info("Extracted %d arcs" % len(arcs))
            self.context.arcs = arcs
        else:
            self.logger.error("Did not extract %d arcs, extracted %d" %
                              (self.config.instrument.NBARS, len(arcs)))

        logstr = ExtractArcs.__module__ + "." + ExtractArcs.__qualname__
        self.action.args.ccddata.header['HISTORY'] = logstr
        self.logger.info(logstr)

        return self.action.args
# END: class ExtractArcs()
