from keckdrpframework.models.arguments import Arguments
from keckdrpframework.primitives.base_img import BaseImg
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer, parse_imsec
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import save_plot

from bokeh.plotting import figure
import ccdproc
import numpy as np
from scipy.stats import sigmaclip
import time


class MakeMasterBias(BaseImg):
    """Generate a master bias image from individual bias frames"""

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can build a stacked frame based on the processing table
        :return:
        """
        # Add to proctab
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix='RAW')
        self.context.proctab.write_proctab()
        # Get bias count
        self.logger.info("Checking precondition for MakeMasterBias")
        self.combine_list = self.context.proctab.n_proctab(
            frame=self.action.args.ccddata, target_type='BIAS',
            target_group=self.action.args.groupid)
        self.logger.info(f"pre condition got {len(self.combine_list)},"
                         f" expecting {self.action.args.min_files}")
        # Did we meet our pre-condition?
        if len(self.combine_list) >= self.action.args.min_files:
            return True
        else:
            return False

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation
        """
        method = 'average'
        suffix = self.action.args.new_type.lower()

        combine_list = list(self.combine_list['OFNAME'])
        # get master bias output name
        mbname = combine_list[0].split('.fits')[0] + '_' + suffix + '.fits'
        stack = []
        stackf = []
        for bias in combine_list:
            stackf.append(bias)
            # using [0] drops the table
            stack.append(kcwi_fits_reader(bias)[0])

        stacked = ccdproc.combine(stack, method=method, sigma_clip=True,
                                  sigma_clip_low_thresh=None,
                                  sigma_clip_high_thresh=2.0)
        stacked.header['IMTYPE'] = self.action.args.new_type
        stacked.header['NSTACK'] = (len(combine_list),
                                    'number of images stacked')
        stacked.header['STCKMETH'] = (method, 'method used for stacking')
        for ii, fname in enumerate(stackf):
            stacked.header['STACKF%d' % (ii + 1)] = (fname, "stack input file")

        # for readnoise stats use 2nd and 3rd bias
        diff = stack[1].data.astype(np.float32) - \
            stack[2].data.astype(np.float32)
        namps = stack[1].header['NVIDINP']
        for ia in range(namps):
            # get gain
            gain = stacked.header['GAIN%d' % (ia + 1)]
            # get amp section
            sec, rfor = parse_imsec(stacked.header['DSEC%d' % (ia + 1)])
            noise = diff[sec[0]:(sec[1]+1), sec[2]:(sec[3]+1)]
            noise = np.reshape(noise, noise.shape[0]*noise.shape[1]) * \
                gain / 1.414
            # get stats on noise
            c, low, upp = sigmaclip(noise, low=3.5, high=3.5)
            bias_rn = c.std()
            self.logger.info("Amp%d read noise from bias in e-: %.3f" %
                             ((ia + 1), bias_rn))
            stacked.header['BIASRN%d' % (ia + 1)] = \
                (float("%.3f" % bias_rn), "RN in e- from bias")
            if self.config.instrument.plot_level >= 1:
                # output filename stub
                biasfnam = "bias_%05d_amp%d_rdnoise" % \
                          (self.action.args.ccddata.header['FRAMENO'], ia+1)
                plabel = '[ Img # %d' % self.action.args.ccddata.header[
                    'FRAMENO']
                plabel += ' (Bias)'
                plabel += ' %s' % self.action.args.ccddata.header['BINNING']
                plabel += ' %s' % self.action.args.ccddata.header['AMPMODE']
                plabel += ' %d' % self.action.args.ccddata.header['GAINMUL']
                plabel += ' %s' % ('fast' if
                                   self.action.args.ccddata.header['CCDMODE']
                                   else 'slow')
                plabel += ' ] '
                hist, edges = np.histogram(noise, range=(low, upp),
                                           density=False, bins=50)
                x = np.linspace(low, upp, 500)
                pdf = np.max(hist)*np.exp(-x**2/(2.*bias_rn**2))
                p = figure(title=plabel+'BIAS NOISE amp %d = %.3f' %
                           (ia+1, bias_rn),
                           x_axis_label='e-', y_axis_label='N',
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
                       fill_color="navy", line_color="white", alpha=0.5)
                p.line(x, pdf, line_color="#ff8888", line_width=4, alpha=0.7,
                       legend_label="PDF")
                p.line([-bias_rn, -bias_rn], [0, np.max(hist)], color='red',
                       legend_label="Sigma")
                p.line([bias_rn, bias_rn], [0, np.max(hist)], color='red')
                p.y_range.start = 0
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
                save_plot(p, filename=biasfnam+".png")

        log_string = MakeMasterBias.__module__
        stacked.header['HISTORY'] = log_string
        self.logger.info(log_string)

        kcwi_fits_writer(stacked, output_file=mbname,
                         output_dir=self.config.instrument.output_directory)
        self.context.proctab.update_proctab(frame=stacked, suffix=suffix,
                                            newtype=self.action.args.new_type)
        self.context.proctab.write_proctab()
        return Arguments(name=mbname)
    # END: class ProcessBias()
