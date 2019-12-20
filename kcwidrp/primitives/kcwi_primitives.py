
from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from keckdrpframework.primitives.base_img import BaseImg
from .kcwi_file_primitives import *
from keckdrpframework.core.bokeh_plotting import bokeh_plot
import ccdproc
from astropy.io import fits as pf
import os

import matplotlib.pyplot as pl
import numpy as np
from bokeh.plotting import figure
from bokeh.models import Range1d
from bokeh.models.markers import X
from bokeh.io import curdoc
from bokeh.util.logconfig import basicConfig, bokeh_logger as bl
import logging

from scipy.signal import find_peaks
from skimage import transform as tf
from scipy.interpolate import interpolate
from scipy.signal import find_peaks
from scipy import signal
from scipy.ndimage import gaussian_filter1d
import time
from multiprocessing import Pool
import pkg_resources


def pascal_shift(coef=None, x0=None):
    """Shift coefficients to a new reference value (X0)

    This should probably go somewhere else, but will be needed here.
    """
    if not coef:
        print("Error, no coefficients for pascal_shift.")
        return None
    if not x0:
        print("Warning, no reference value (x0) supplied")
        return coef
    if len(coef) == 7:
        usecoeff = list(reversed(coef))
        fincoeff = [0.] * 7
    else:
        if len(coef) > 7:
            print("Warning - this routine only handles up to 7 coefficients.")
            usecoeff = list(reversed(coef[0:7]))
            fincoeff = [0.] * len(coef)
        else:
            usecoeff = [0.] * 7
            fincoeff = usecoeff
            for ic, c in enumerate(coef):
                usecoeff[len(coef)-(ic+1)] = coef[ic]
    # get reference values
    x01 = x0
    x02 = x0**2
    x03 = x0**3
    x04 = x0**4
    x05 = x0**5
    x06 = x0**6
    # use Pascal's Triangle to shift coefficients
    fincoeff[0] = usecoeff[0] - usecoeff[1] * x01 + usecoeff[2] * x02 \
        - usecoeff[3] * x03 + usecoeff[4] * x04 - usecoeff[5] * x05 \
        + usecoeff[6] * x06

    fincoeff[1] = usecoeff[1] - 2.0 * usecoeff[2] * x01 \
        + 3.0 * usecoeff[3] * x02 - 4.0 * usecoeff[4] * x03 \
        + 5.0 * usecoeff[5] * x04 - 6.0 * usecoeff[6] * x05

    fincoeff[2] = usecoeff[2] - 3.0 * usecoeff[3] * x01 \
        + 6.0 * usecoeff[4] * x02 - 10.0 * usecoeff[5] * x03 \
        + 15.0 * usecoeff[6] * x04

    fincoeff[3] = usecoeff[3] - 4.0 * usecoeff[4] * x01 \
        + 10.0 * usecoeff[5] * x02 - 20.0 * usecoeff[6] * x03

    fincoeff[4] = usecoeff[4] - 5.0 * usecoeff[5] * x01 \
        + 15.0 * usecoeff[6] * x02

    fincoeff[5] = usecoeff[5] - 6.0 * usecoeff[6] * x01

    fincoeff[6] = usecoeff[6]
    # Trim if needed
    if len(coef) < 7:
        fincoeff = fincoeff[0:len(coef)]
    # Reverse for python
    return list(reversed(fincoeff))

class subtract_overscan(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)


    def _perform(self):
        # image sections for each amp
        bsec, dsec, tsec, direc = self.action.args.map_ccd
        namps = len(bsec)
        # polynomial fit order
        if namps == 4:
            porder = 2
        else:
            porder = 7
        # header keyword to update
        key = 'OSCANSUB'
        # is it performed?
        performed = False
        # loop over amps
        #session = self.context.bokeh_session

        for ia in range(namps):
            # get gain
            gain = self.context.data_set.get_info_column(self.action.args.name, 'GAIN%d' % (ia + 1))
            # check if we have enough data to fit
            if (bsec[ia][3] - bsec[ia][2]) > self.context.config.instrument.minoscanpix:
                # pull out an overscan vector
                x0 = bsec[ia][2] + self.context.config.instrument.oscanbuf
                x1 = bsec[ia][3] - self.context.config.instrument.oscanbuf
                y0 = bsec[ia][0]
                y1 = bsec[ia][1] + 1
                osvec = np.nanmedian(self.action.args.ccddata.data[y0:y1, x0:x1], axis=1)
                nsam = x1 - x0
                xx = np.arange(len(osvec), dtype=np.float)
                # fit it, avoiding first 50 px
                if direc[ia]:
                    # forward read skips first 50 px
                    oscoef = np.polyfit(xx[50:], osvec[50:], porder)
                else:
                    # reverse read skips last 50 px
                    oscoef = np.polyfit(xx[:-50], osvec[:-50], porder)
                # generate fitted overscan vector for full range
                osfit = np.polyval(oscoef, xx)
                # calculate residuals
                resid = (osvec - osfit) * math.sqrt(nsam) * gain / 1.414
                sdrs = float("%.3f" % np.std(resid))
                self.logger.info("Amp%d Read noise from oscan in e-: %.3f" %
                          ((ia + 1), sdrs))
                self.action.args.ccddata.header['OSCNRN%d' % (ia + 1)] = (sdrs, "amp%d RN in e- from oscan" % (ia + 1))

                if self.context.config.plot_level>=1:
                    x=np.arange(len(osvec))
                    p = figure(title='Overscan amp %d' % (ia+1), x_axis_label='x', y_axis_label='counts',
                               plot_width=self.context.config.instrument.plot_width, plot_height=self.context.config.instrument.plot_height)
                    p.line(x,osvec)
                    p.line(x,osfit)
                    bokeh_plot(p)
                    if self.context.config.plot_level >= 2:
                        input("Next? <cr>: ")
                    else:
                        time.sleep(self.context.config.instrument.plot_pause)
                # subtract it
                for ix in range(dsec[ia][2], dsec[ia][3] + 1):
                    self.action.args.ccddata.data[y0:y1, ix] = \
                        self.action.args.ccddata.data[y0:y1, ix] - osfit
            performed = True
            #else:
            #    self.log.info("not enough overscan px to fit amp %d")

        if performed:
            self.action.args.ccddata.header[key] = (True, 'Overscan subtracted')
        else:
            self.action.args.ccddata.header[key] = (False, 'Overscan subtracted')
        logstr = self.__module__ + "." + self.__class__.__name__
        self.action.args.ccddata.header['HISTORY'] = logstr
        #self.log.info(self.subtract_oscan.__qualname__)

        #new_name = "%s_ovsc" % self.action.args.name
        kcwi_fits_writer(self.action.args.ccddata, table=self.action.args.table, output_file=self.action.args.name, suffix="oscan")
        return self.action.args


class trim_overscan(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)

    def _perform(self):

        # parameters
        # image sections for each amp
        bsec, dsec, tsec, direc = self.action.args.map_ccd
        namps = len(bsec)
        # header keyword to update
        key = 'OSCANTRM'
        # get output image dimensions
        max_sec = max(tsec)
        # create new blank image
        new = np.zeros((max_sec[1]+1, max_sec[3]+1), dtype=np.float64)
        # loop over amps
        for ia in range(namps):
            # input range indices
            yi0 = dsec[ia][0]
            yi1 = dsec[ia][1] + 1
            xi0 = dsec[ia][2]
            xi1 = dsec[ia][3] + 1
            # output range indices
            yo0 = tsec[ia][0]
            yo1 = tsec[ia][1] + 1
            xo0 = tsec[ia][2]
            xo1 = tsec[ia][3] + 1
            # transfer to new image
            new[yo0:yo1, xo0:xo1] = self.action.args.ccddata.data[yi0:yi1, xi0:xi1]
            # update amp section
            sec = "[%d:" % (xo0+1)
            sec += "%d," % xo1
            sec += "%d:" % (yo0+1)
            sec += "%d]" % yo1
            self.logger.info("ADDING ATSEC%d" % (ia + 1))
            self.action.args.ccddata.header['ATSEC%d' % (ia+1)] = sec
            # remove obsolete sections
            self.action.args.ccddata.header.pop('ASEC%d' % (ia + 1))
            self.action.args.ccddata.header.pop('BSEC%d' % (ia + 1))
            self.action.args.ccddata.header.pop('DSEC%d' % (ia + 1))
            self.action.args.ccddata.header.pop('CSEC%d' % (ia + 1))
        # update with new image
        self.action.args.ccddata.data = new
        self.action.args.ccddata.header['NAXIS1'] = max_sec[3] + 1
        self.action.args.ccddata.header['NAXIS2'] = max_sec[1] + 1
        self.action.args.ccddata.header[key] = (True, "Overscan trimmed")

        logstr = self.__module__ + "." + \
                 self.__class__.__name__
        self.action.args.ccddata.header['HISTORY'] = logstr
        #new_name = "%s_trim" % self.action.args.name

        kcwi_fits_writer(self.action.args.ccddata, table=self.action.args.table, output_file=self.action.args.name, suffix="trim")
        return self.action.args
        #self.log.info(self.trim_oscan.__qualname__)


class correct_gain(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)

    def _perform(self):
        #print(self.action.args.ccddata.header)
        namps = self.action.args.namps
        for ia in range(namps):
            # get amp section
            section = self.action.args.ccddata.header['ATSEC%d' % (ia +1)]
            sec, rfor = parse_imsec(section)
            # get gain for this amp
            gain = self.context.data_set.get_info_column(self.action.args.name,'GAIN%d' % (ia + 1))
            self.logger.info("Applying gain correction of %.3f in section %s" %(gain, self.action.args.ccddata.header['ATSEC%d' % (ia + 1)]))
            self.action.args.ccddata.data[sec[0]:(sec[1]+1), sec[2]:(sec[3]+1)] *= gain
            #sliced_ccddata.header=self.action.args.header

            #multiplied_ccddata = sliced_ccddata.multiply(gain)
            #mu
            #self.action.args.ccddata=multiplied_ccddata

        self.action.args.ccddata.header['GAINCOR'] = (True, "Gain corrected")
        self.action.args.ccddata.header['BUNIT'] = ('electron','Units set to electrons')
        self.action.args.ccddata.unit = 'electron'

        #logstr = self.correct_gain.__module__ + "." + \
        #         self.correct_gain.__qualname__
        #self.frame.header['HISTORY'] = logstr
        #self.log.info(self.correct_gain.__qualname__)
        #self.logger.info("Writing test file")
        #self.action.args.ccddata.write('test.fits')
        kcwi_fits_writer(self.action.args.ccddata, table=self.action.args.table, output_file=self.action.args.name, suffix="int")
        return self.action.args

class process_bias(BaseImg):

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)

    def hello(self):
        """
        Checks is we can build a stacked  frame
        Expected arguments:
            want_type: ie. BIAS
            min_files: ie 10
            new_type: ie MASTER_BIAS
            new_file_name: master_bias.fits

        """
        try:
            args = self.action.args
            df = self.context.data_set.data_table
            files = df[(df.IMTYPE == args.want_type) & (df.GROUPID == args.groupid)]
            nfiles = len(files)


            self.logger.info(f"pre condition got {nfiles}, expecting {args.min_files}")
            if nfiles < 1 or nfiles < args.min_files:
                return False
            return True
        except Exception as e:
            self.logger.error(f"Exception in base_ccd_primitive: {e}")
            return False

    def _pre_condition(self):
        """
        Checks if we can build a stacked frame based on the processing table
        :return:
        """
        # get current group id
        self.logger.info("Checking precondition for process_bias")
        self.combine_list = self.context.proctab.n_proctab(frame = self.action.args.ccddata, target_type='BIAS', target_group=self.action.args.groupid)
        self.logger.info(f"pre condition got {len(self.combine_list)}, expecting {self.action.args.min_files}")
        # create master bias
        if len(self.combine_list) >= self.action.args.min_files:
            return True
        else:
            return False

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation.
        """
        args = self.action.args

        #df = self.context.data_set.data_table
        #files = list(df[(df.IMTYPE == args.want_type) & (df.GROUPID == args.groupid)].index)
        #print("BEFORE:", files)
        #print("NOW:", self.combine_list)
        combine_list = list(self.combine_list['OFNAME'])
        stack = []
        for bias in combine_list:
            # using [0] drops the table
            stack.append(kcwi_fits_reader(bias)[0])

        stacked = ccdproc.combine(stack)
        stacked.header.IMTYPE=args.new_type

        kcwi_fits_writer(stacked, output_file = args.new_file_name)
        self.context.proctab.update_proctab(frame=stacked, suffix='master_bias', newtype='MBIAS')
        return Arguments(name=args.new_file_name)

class process_contbars(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)

    def _perform(self):
        return self.action.args

class find_bars(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        basicConfig(level=logging.ERROR)

    def _perform(self):
        self.logger.info("Finding continuum bars")
        # Do we plot?
        if self.context.config.plot_level >= 1:
            do_plot = True
        else:
            do_plot = False
        # initialize
        midcntr = []
        # get image dimensions
        nx = self.action.args.ccddata.data.shape[1]
        ny = self.action.args.ccddata.data.shape[0]
        # get binning
        ybin = self.action.args.ybinsize
        win = int(10 / ybin)
        # select from center rows of image
        midy = int(ny / 2)
        midvec = np.median(self.action.args.ccddata.data[(midy-win):(midy+win+1), :], axis=0)
        # set threshold for peak finding
        midavg = np.average(midvec)
        self.logger.info("peak threshold = %f" % midavg)
        # find peaks above threshold
        midpeaks, _ = find_peaks(midvec, height=midavg)
        # do we have the requisite number?
        if len(midpeaks) != self.context.config.instrument.NBARS:
            self.logger.error("Did not find %d peaks: n peaks = %d"
                           % (self.context.config.instrument.NBARS, len(midpeaks)))
        else:
            self.logger.info("found %d bars" % len(midpeaks))
            if do_plot:
                # plot the peak positions
                x = np.arange(len(midvec))
                #pl.plot(midvec, '-')
                p = figure(title="Img %d, Thresh = %.2f" % (self.action.args.ccddata.header['FRAMENO'], midavg), \
                    x_axis_label='CCD X (px)', y_axis_label='e-',
                           plot_width=self.context.config.instrument.plot_width, plot_height=self.context.config.instrument.plot_height)
                p.line(x, midvec, color='blue')
                p.scatter(midpeaks, midvec[midpeaks], marker='x', color='red')
                p.line([0, nx], [midavg, midavg], color='grey', line_dash='dashed')
                bokeh_plot(p)
                time.sleep(self.context.config.instrument.plot_pause)
            # calculate the bar centroids
            plotting_vector_x=[]
            plotting_vector_y=[]
            for peak in midpeaks:
                xs = list(range(peak-win, peak+win+1))
                ys = midvec[xs] - np.nanmin(midvec[xs])
                xc = np.sum(xs*ys) / np.sum(ys)
                midcntr.append(xc)
                plotting_vector_x.append(xc)
                plotting_vector_y.append(midavg)
                plotting_vector_x.append(xc)
                plotting_vector_y.append(midvec[peak])
                plotting_vector_x.append(xc)
                plotting_vector_y.append(midavg)

                #p=figure()
            if do_plot:
                p.line(plotting_vector_x, plotting_vector_y, color='grey')
                #p.line([xc, xc], [midavg, midvec[peak]], color='grey')

            if do_plot:
                #p = figure()
                p.scatter(midcntr, midvec[midpeaks], marker='x', color='green')
                bokeh_plot(p)
                if self.context.config.plot_level >= 2:
                    input("next: ")
                else:
                    time.sleep(self.context.config.instrument.plot_pause)
            self.logger.info("Found middle centroids for continuum bars")
        self.action.args.midcntr = midcntr
        self.action.args.midrow = midy
        self.action.args.win = win
        return self.action.args


class trace_bars(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)

    def _perform(self):
        self.logger.info("Tracing continuum bars")
        if self.context.config.plot_level >= 1:
            do_plot = True
            pl.ion()
        else:
            do_plot = False
        if len(self.action.args.midcntr) < 1:
            self.logger.error("No bars found")
        else:
            # initialize
            samp = int(80 / self.action.args.ybinsize)
            win = self.action.args.win
            xi = []     # x input
            xo = []     # x output
            yi = []     # y input (and output)
            barid = []  # bar id number
            slid = []   # slice id number
            # loop over bars
            for barn, barx in enumerate(self.action.args.midcntr):
                # nearest pixel to bar center
                barxi = int(barx + 0.5)
                self.logger.info("bar number %d is at %.3f" % (barn, barx))
                # middle row data
                xi.append(barx)
                xo.append(barx)
                yi.append(self.action.args.midrow)
                barid.append(barn)
                slid.append(int(barn/5))
                # trace up
                samy = self.action.args.midrow + samp
                done = False
                while samy < (self.action.args.ccddata.data.shape[0] - win) and not done:
                    ys = np.median(
                        self.action.args.ccddata.data[(samy - win):(samy + win + 1),
                                        (barxi - win):(barxi + win + 1)],
                        axis=0)
                    ys = ys - np.nanmin(ys)
                    xs = list(range(barxi - win, barxi + win + 1))
                    xc = np.sum(xs * ys) / np.sum(ys)
                    if np.nanmax(ys) > 255:
                        xi.append(xc)
                        xo.append(barx)
                        yi.append(samy)
                        barid.append(barn)
                        slid.append(int(barn/5))
                    else:
                        done = True
                    samy += samp
                # trace down
                samy = self.action.args.midrow - samp
                done = False
                while samy >= win and not done:
                    ys = np.median(
                        self.action.args.ccddata.data[(samy - win):(samy + win + 1),
                                        (barxi - win):(barxi + win + 1)],
                        axis=0)
                    ys = ys - np.nanmin(ys)
                    xs = list(range(barxi - win, barxi + win + 1))
                    xc = np.sum(xs * ys) / np.sum(ys)
                    if np.nanmax(ys) > 255:
                        xi.append(xc)
                        xo.append(barx)
                        yi.append(samy)
                        barid.append(barn)
                        slid.append(int(barn / 5))
                    else:
                        done = True
                    # disable for now
                    samy -= samp
            # end loop over bars
            # create source and destination coords
            yo = yi
            dst = np.column_stack((xi, yi))
            src = np.column_stack((xo, yo))

            if do_plot:
                # plot them
                # pl.ioff()
                p = figure(title="Img %d" % self.action.args.ccddata.header['FRAMENO'],
                    x_axis_label="CCD X (px)", y_axis_label="CCD Y (px)",
                           plot_width=self.context.config.instrument.plot_width, plot_height=self.context.config.instrument.plot_height)
                p.scatter(xi, yi, marker='x', size=2, color='blue')
                #pl.plot(xi, yi, 'x', ms=0.5)
                p.scatter(self.action.args.midcntr, [self.action.args.midrow]*120, color='red')
                #pl.plot(self.action.args.midcntr, [self.action.args.midrow]*120, 'x', color='red')
                bokeh_plot(p)
                if self.context.config.plot_level >= 2:
                    input("next: ")
                else:
                    time.sleep(self.context.config.instrument.plot_pause)
            trace = {
                'src': src,
                'dst': dst,
                'barid': barid,
                'slid': slid,
                'MIDROW': self.action.args.midrow,
                'WINDOW': self.action.args.win}

            # in this line we pass the trace information to an argument instead of writing it to a table
            self.context.trace = trace
            write_table(table = [src, dst, barid, slid],
                        names = ('src', 'dst', 'barid', 'slid'),
                        #suffix='trace',
                        output_dir = os.path.dirname(self.action.args.name),
                        output_name = 'trace_table.fits',
                        comment = ['Source and destination fiducial points',
                                      'Derived from KCWI continuum bars images',
                                      'For defining spatial transformation'],
                        keywords={'MIDROW': self.action.args.midrow,
                                    'WINDOW': self.action.args.win})

            #self.write_table(table=[src, dst, barid, slid],
            #                 names=('src', 'dst', 'barid', 'slid'),
            #                 suffix='trace',
            #                 comment=['Source and destination fiducial points',
            #                          'Derived from KCWI continuum bars images',
            #                          'For defining spatial transformation'],
            #                 keywords={'MIDROW': self.midrow,
            #                           'WINDOW': self.win})
            if True:
                # fit transform
                self.logger.info("Fitting spatial control points")
                tform = tf.estimate_transform('polynomial', src, dst, order=3)
                self.logger.info("Transforming bars image")
                warped = tf.warp(self.action.args.ccddata.data, tform)
                # write out warped image
                self.action.args.ccddata.data = warped
                kcwi_fits_writer(self.action.args.ccddata, self.action.args.table, output_file=self.action.args.name, suffix='warped')
                self.logger.info("Transformed bars produced")
            return self.action.args



class extract_arcs(BasePrimitive):

    def __init__(self, action, context):
            BasePrimitive.__init__(self, action, context)

    def _perform(self):
        self.logger.info("Extracting arc spectra")
        # Find  and read control points from continuum bars
        if hasattr(self.context, 'trace'):
            trace = self.context.trace
            midrow = trace['MIDROW']
            win = trace['WINDOW']
        else:
            trace = read_table(input_dir = os.path.dirname(self.action.args.name), file_name = "trace_table.fits")
            midrow = trace.meta['MIDROW']
            win = trace.meta['WINDOW']
            self.context.trace = trace
        src = trace['src']  # source control points
        dst = trace['dst']  # destination control points
        barid = trace['barid']
        slid = trace['slid']



        self.logger.info("Fitting spatial control points")
        tform = tf.estimate_transform('polynomial', src, dst, order=3)

        self.logger.info("Transforming arc image")
        warped = tf.warp(self.action.args.ccddata.data, tform)
        # Write warped arcs if requested
        #if self.frame.saveintims():
        #    # write out warped image
        #    self.frame.data = warped
        #    self.write_image(suffix='warped')
        #    self.log.info("Transformed arcs produced")
        # extract arcs
        self.logger.info("Extracting arcs")
        arcs = []
        for xyi, xy in enumerate(src):
            if xy[1] == midrow:
                xi = int(xy[0]+0.5)
                arc = np.median(
                    warped[:, (xi - win):(xi + win + 1)], axis=1)
                arc = arc - np.nanmin(arc[100:-100])    # avoid ends
                arcs.append(arc)
        # Did we get the correct number of arcs?
        if len(arcs) == self.context.config.instrument.NBARS:
            self.logger.info("Extracted %d arcs" % len(arcs))
            self.context.arcs = arcs
        else:
            self.logger.error("Did not extract %d arcs, extracted %d" %
                           (self.context.config.instrument.NBARS, len(arcs)))
        return self.action.args

class arc_offsets(BasePrimitive):

    def __init__(self, action, context):
            BasePrimitive.__init__(self, action, context)

    def _perform(self):
        self.logger.info("Finding inter-bar offsets")
        arcs = self.context.arcs
        if arcs is not None:
            # Do we plot?
            if self.context.config.plot_level >= 2:
                do_plot = True
            else:
                do_plot = False
            # Compare with reference arc
            refarc = arcs[self.context.config.instrument.REFBAR][:]
            # number of cross-correlation samples (avoiding ends)
            nsamp = len(refarc[10:-10])
            # possible offsets
            offar = np.arange(1-nsamp, nsamp)
            # Collect offsets
            offsets = []
            for na, arc in enumerate(arcs):
                # Cross-correlate, avoiding junk on the ends
                xcorr = np.correlate(refarc[10:-10], arc[10:-10], mode='full')
                # Calculate offset
                offset = offar[xcorr.argmax()]
                offsets.append(offset)
                self.logger.info("Arc %d Slice %d XCorr shift = %d" %
                              (na, int(na/5), offset))
                # display if requested
                if do_plot:
                    p=figure(title="Arc %d Slice %d XCorr, Shift = %d" % (na, int(na/5), offset),
                             x_axis_label="CCD y (px)", y_axis_label="e-",
                             plot_width=self.context.config.instrument.plot_width, plot_height=self.context.config.instrument.plot_height)
                    x = range(len(refarc))
                    p.line(x,refarc, color='green')
                    p.line(x, np.roll(arc, offset), color='red')
                    #pl.ylim(bottom=0.)
                    #pl.xlabel("CCD y (px)")
                    #pl.ylabel("e-")
                    #pl.title("Arc %d Slice %d XCorr, Shift = %d" %
                    #         (na, int(na/5), offset))
                    #pl.show()
                    bokeh_plot(p)
                    q = input("<cr> - Next, q to quit: ")
                    if 'Q' in q.upper():
                        do_plot = False
            self.context.baroffs = offsets
        else:
            self.log.error("No extracted arcs found")
        return self.action.args


class calc_prelim_disp(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)

    def _perform(self):
        # get binning
        ybin = self.action.args.ybinsize
        # 0 - compute alpha
        prelim_alpha = self.action.args.grangle - 13.0 - self.action.args.adjang
        # 1 - compute preliminary angle of diffraction
        prelim_beta = self.action.args.camangle - prelim_alpha
        # 2 - compute preliminary dispersion
        prelim_disp = math.cos(prelim_beta/math.degrees(1.)) / \
            self.action.args.rho / self.context.config.instrument.FCAM * (self.context.config.instrument.PIX*ybin) * 1.e4
        prelim_disp *= math.cos(self.context.config.instrument.GAMMA/math.degrees(1.))
        self.logger.info("Initial alpha, beta (deg): %.3f, %.3f" %
                      (prelim_alpha, prelim_beta))
        self.logger.info("Initial calculated dispersion (A/binned pix): %.3f" %
                      prelim_disp)
        self.context.prelim_disp = prelim_disp
        return self.action.args


class read_atlas(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)

    def _perform(self):
        # What lamp are we using?
        lamp = self.action.args.illum
        # rez factor
        if 'fear' in lamp.lower():
            rezfact = 0.5
        else:
            rezfact = 1.0
        #atpath = os.path.join("../data", "%s.fits" % lamp.lower())
        # Does the atlas file exist?
        path = "data/%s.fits" % lamp.lower()  # always use slash
        pkg = __name__.split('.')[0]
        atpath = pkg_resources.resource_filename(pkg, path)
        if os.path.exists(atpath):
            self.logger.info("Reading atlas spectrum in: %s" % atpath)
        else:
            self.logger.error("Atlas spectrum not found for %s" % atpath)
        # Read the atlas
        ff = pf.open(atpath)
        reflux = ff[0].data
        refdisp = ff[0].header['CDELT1']
        refwav = np.arange(0, len(reflux)) * refdisp + ff[0].header['CRVAL1']
        ff.close()
        # Convolve with appropriate Gaussian
        reflux = gaussian_filter1d(reflux, self.action.args.atres*rezfact)
        # Observed arc spectrum
        obsarc = self.context.arcs[self.context.config.instrument.REFBAR]
        # Preliminary wavelength solution
        xvals = np.arange(0, len(obsarc)) - int(len(obsarc)/2)
        obswav = xvals * self.context.prelim_disp + self.action.args.cwave
        # Get central third
        minow = int(len(obsarc)/3)
        maxow = int(2.*len(obsarc)/3)
        # Unless we are low dispersion, then get central 3 5ths
        if 'BL' in self.action.args.grating or 'RL' in self.action.args.grating:
            minow = int(len(obsarc)/5)
            maxow = int(4.*len(obsarc)/5)
        minwav = obswav[minow]
        maxwav = obswav[maxow]
        # Get corresponding ref range
        minrw = [i for i, v in enumerate(refwav) if v >= minwav][0]
        maxrw = [i for i, v in enumerate(refwav) if v <= maxwav][-1]
        # Subsample for cross-correlation
        cc_obsarc = obsarc[minow:maxow]
        cc_obswav = obswav[minow:maxow]
        cc_reflux = reflux[minrw:maxrw]
        cc_refwav = refwav[minrw:maxrw]
        # Resample onto reference wavelength scale
        obsint = interpolate.interp1d(cc_obswav, cc_obsarc, kind='cubic',
                                      bounds_error=False,
                                      fill_value='extrapolate'
                                      )
        cc_obsarc = obsint(cc_refwav)
        # Apply cosign bell taper to both
        cc_obsarc *= signal.windows.tukey(len(cc_obsarc),
                                          alpha=self.context.config.instrument.TAPERFRAC)
        cc_reflux *= signal.windows.tukey(len(cc_reflux),
                                          alpha=self.context.config.instrument.TAPERFRAC)
        nsamp = len(cc_refwav)
        offar = np.arange(1 - nsamp, nsamp)
        # Cross-correlate
        xcorr = np.correlate(cc_obsarc, cc_reflux, mode='full')
        # Get central region
        x0c = int(len(xcorr)/3)
        x1c = int(2*(len(xcorr)/3))
        xcorr_central = xcorr[x0c:x1c]
        offar_central = offar[x0c:x1c]
        # Calculate offset
        offset_pix = offar_central[xcorr_central.argmax()]
        offset_wav = offset_pix * refdisp
        self.logger.info("Initial arc-atlas offset (px, Ang): %d, %.1f" %
                      (offset_pix, offset_wav))
        if self.context.config.plot_level >= 1:
            #if self.frame.inter() >= 2:
            #    pl.ion()
            #else:
            #    pl.ioff()
            # Plot
            p = figure(title="Img # %d (%s), Offset = %d px" % (self.action.args.ccddata.header['FRAMENO'], lamp, offset_pix),
                       x_axis_label="Offset(px)", y_axis_label="X-corr",
                       plot_width=self.context.config.instrument.plot_width, plot_height=self.context.config.instrument.plot_height)

            p.line(offar_central, xcorr_central)
            ylim_min = min(xcorr_central)
            ylim_max = max(xcorr_central)
            #ylim = pl.gca().get_ylim()
            p.line([offset_pix, offset_pix], [ylim_min, ylim_max], color='green')
            bokeh_plot(p)
            if self.context.config.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.context.config.instrument.plot_pause)
            # Get central wavelength
            cwave = self.action.args.cwave
            # Set up offset tweaking
            q = 'test'
            while q:
                # Plot the two spectra
                p = figure(title="Img # %d (%s), Offset = %.1f Ang (%d px)" % (self.action.args.ccddata.header['FRAMENO'], lamp, offset_wav, offset_pix),
                           x_axis_label="Wave(A)", y_axis_label="Rel. Flux",
                           plot_width=self.context.config.instrument.plot_width,
                           plot_height=self.context.config.instrument.plot_height)
                p.line(obswav[minow:maxow] - offset_wav,
                        obsarc[minow:maxow]/np.nanmax(obsarc[minow:maxow]), legend="ref bar (%d)" % self.context.config.instrument.REFBAR)
                p.line(refwav[minrw:maxrw],
                        reflux[minrw:maxrw]/np.nanmax(reflux[minrw:maxrw]), color="red", legend="Atlas")
                p.x_range=Range1d(np.nanmin(obswav[minow:maxow]),np.nanmax(obswav[minow:maxow]))
                ylim_min = min(obsarc[minow:maxow]/np.nanmax(obsarc[minow:maxow]))
                ylim_max = max(obsarc[minow:maxow]/np.nanmax(obsarc[minow:maxow]))
                #ylim = pl.gca().get_ylim()
                p.line([cwave, cwave], [ylim_min, ylim_max], color="green", legend="CWAVE")
                bokeh_plot(p)

                if self.context.config.plot_level >= 2:
                    q = input("Enter: <cr> - next, new offset (px): ")
                    if q:
                        offset_pix = int(q)
                        offset_wav = offset_pix * refdisp
                else:
                    time.sleep(self.context.config.instrument.plot_pause)
                    q = None
            self.logger.info("Final   arc-atlas offset (px, Ang): %d, %.1f" %
                          (offset_pix, offset_wav))
        # Store atlas spectrum
        self.action.args.reflux = reflux
        self.action.args.refwave = refwav
        # Store offsets
        self.action.args.offset_pix = offset_pix
        self.action.args.offset_wave = offset_wav
        # Store reference dispersion
        self.action.args.refdisp = refdisp
        # Store central limits
        self.action.args.minrow = minow
        self.action.args.maxrow = maxow
        # Store x values
        self.action.args.xvals = xvals
        self.action.args.x0 = int(len(obsarc)/2)
        return self.action.args

def myhelper(argument):

    b = argument['b']
    bs = argument['bs']
    # wavelength coefficients
    coeff = [0., 0., 0., 0., 0.]
    # container for maxima, shifts
    maxima = []
    shifts = []
    # get sub spectrum for this bar
    subspec = bs[argument['minrow']:argument['maxrow']]
    # now loop over dispersions
    for di, disp in enumerate(argument['disps']):
        # populate the coefficients
        coeff[4] = argument['p0'][b]
        coeff[3] = disp
        cosbeta = disp / (argument['PIX'] * argument['ybin']) * argument['rho'] * \
                  argument['FCAM'] * 1.e-4
        if cosbeta > 1.:
            cosbeta = 1.
        beta = math.acos(cosbeta)
        coeff[2] = -(argument['PIX'] * argument['ybin'] / argument['FCAM']) ** 2 * \
                   math.sin(beta) / 2. / argument['rho'] * 1.e4
        coeff[1] = -(argument['PIX']* argument['ybin'] / argument['FCAM']) ** 3 * \
                   math.cos(beta) / 6. / argument['rho'] * 1.e4
        coeff[0] = (argument['PIX'] * argument['ybin'] / argument['FCAM']) ** 4 * \
                   math.sin(beta) / 24. / argument['rho'] * 1.e4
        # what are the min and max wavelengths to consider?
        wl0 = np.polyval(coeff, argument['xvals'][argument['minrow']])
        wl1 = np.polyval(coeff, argument['xvals'][argument['maxrow']])
        minwvl = np.nanmin([wl0, wl1])
        maxwvl = np.nanmax([wl0, wl1])
        # where will we need to interpolate to cross-correlate?
        minrw = [i for i, v in enumerate(argument['refwave'])
                 if v >= minwvl][0]
        maxrw = [i for i, v in enumerate(argument['refwave'])
                 if v <= maxwvl][-1]
        subrefwvl = argument['refwave'][minrw:maxrw]
        subrefspec = argument['reflux'][minrw:maxrw]
        # get bell cosine taper to avoid nasty edge effects
        tkwgt = signal.windows.tukey(len(subrefspec),
                                     alpha=argument['taperfrac'])
        # apply taper to atlas spectrum
        subrefspec *= tkwgt
        # adjust wavelengths
        waves = np.polyval(coeff, argument['subxvals'])
        # interpolate the bar spectrum
        obsint = interpolate.interp1d(waves, subspec, kind='cubic',
                                      bounds_error=False,
                                      fill_value='extrapolate')
        intspec = obsint(subrefwvl)
        # apply taper to bar spectrum
        intspec *= tkwgt
        # get a label
        # cross correlate the interpolated spectrum with the atlas spec
        nsamp = len(subrefwvl)
        offar = np.arange(1 - nsamp, nsamp)
        # Cross-correlate
        xcorr = np.correlate(intspec, subrefspec, mode='full')
        # Get central region
        x0c = int(len(xcorr) / 3)
        x1c = int(2 * (len(xcorr) / 3))
        xcorr_central = xcorr[x0c:x1c]
        offar_central = offar[x0c:x1c]
        # Calculate offset
        maxima.append(xcorr_central[xcorr_central.argmax()])
        shifts.append(offar_central[xcorr_central.argmax()])
    # Get interpolations
    int_max = interpolate.interp1d(argument['disps'], maxima, kind='cubic',
                                   bounds_error=False,
                                   fill_value='extrapolate')
    int_shift = interpolate.interp1d(argument['disps'], shifts, kind='cubic',
                                     bounds_error=False,
                                     fill_value='extrapolate')
    xdisps = np.linspace(min(argument['disps']), max(argument['disps']), num=argument['nn'] * 100)
    # get peak values
    maxima_res = int_max(xdisps)
    shifts_res = int_shift(xdisps) * argument['refdisp']
    bardisp = xdisps[maxima_res.argmax()]
    barshift = shifts_res[maxima_res.argmax()]
    # update coeffs
    coeff[4] = argument['p0'][b] - barshift
    coeff[3] = bardisp
    cosbeta = coeff[3] / (argument['PIX'] * argument['ybin']) * argument['rho'] * \
              argument['FCAM'] * 1.e-4
    if cosbeta > 1.:
        cosbeta = 1.
    beta = math.acos(cosbeta)
    coeff[2] = -(argument['PIX'] * argument['ybin'] / argument['FCAM']) ** 2 * \
               math.sin(beta) / 2. / argument['rho'] * 1.e4
    coeff[1] = -(argument['PIX'] * argument['ybin'] / argument['FCAM']) ** 3 * \
               math.cos(beta) / 6. / argument['rho'] * 1.e4
    coeff[0] = (argument['PIX'] * argument['ybin'] / argument['FCAM']) ** 4 * \
               math.sin(beta) / 24. / argument['rho'] * 1.e4
    scoeff = pascal_shift(coeff, argument['x0'])
    print("Central Fit: Bar#, Cdisp, Coefs: "
                     "%3d  %.4f  %.2f  %.4f  %13.5e %13.5e" %
                     (b, bardisp, scoeff[4], scoeff[3], scoeff[2],
                      scoeff[1]))
    # store central values
    #centwave.append(coeff[4])
    #centdisp.append(coeff[3])
    # Store results
    return b, scoeff, coeff[4], coeff[3]

class fit_center(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.action.args.centcoeff=[]

    def _perform(self):
        """ Fit central region

        At this point we have the offsets between bars and the approximate
        offset from the reference bar to the atlas spectrum and the approximate
        dispersion.
        """
        self.logger.info("Finding wavelength solution for central region")
        # Are we interactive?
        #if KcwiConf.INTER >= 2:
        #    do_inter = True
        #    pl.ion()
        #else:
        #    do_inter = False
        do_inter=False
        # image label
        imlab = "Img # %d (%s) Sl: %s Fl: %s Gr: %s" % \
                (self.action.args.ccddata.header['FRAMENO'], self.action.args.illum,
                 self.action.args.ifuname, self.action.args.filter,
                 self.action.args.grating)
        # y binning
        ybin = self.action.args.ybinsize
        # let's populate the 0 points vector
        p0 = self.action.args.cwave + np.array(self.context.baroffs) * self.context.prelim_disp \
            - self.action.args.offset_wave
        # next we are going to brute-force scan around the preliminary
        # dispersion for a better solution. We will wander 5% away from it.
        max_ddisp = 0.05    # fraction
        # we will try nn values
        nn = (int(max_ddisp*abs(self.context.prelim_disp)/self.action.args.refdisp*(
                self.action.args.maxrow-self.action.args.minrow)/3.0))
        if nn < 10:
            nn = 10
        if nn > 25:
            nn = 25
        self.logger.info("N disp. samples: %d" % nn)
        # dispersions to try
        disps = self.context.prelim_disp * (1.0 + max_ddisp *
                                    (np.arange(0, nn+1) - nn/2.) * 2.0 / nn)
        # containers for bar-specific values
        bardisp = []
        barshift = []
        centwave = []
        centdisp = []

        # values for central fit
        subxvals = self.action.args.xvals[self.action.args.minrow:self.action.args.maxrow]
        # loop over bars


        my_arguments = []
        for b, bs in enumerate(self.context.arcs):
            arguments={}
            arguments['b']=b
            arguments['bs'] = bs
            arguments['minrow'] = self.action.args.minrow
            arguments['maxrow'] = self.action.args.maxrow
            arguments['disps'] = disps
            arguments['p0'] = p0
            arguments['PIX']= self.context.config.instrument.PIX
            arguments['ybin'] = ybin
            arguments['rho'] = self.action.args.rho
            arguments['FCAM'] = self.context.config.instrument.FCAM
            arguments['xvals'] = self.action.args.xvals
            arguments['refwave'] = self.action.args.refwave
            arguments['reflux'] = self.action.args.reflux
            arguments['taperfrac'] = self.context.config.instrument.TAPERFRAC
            arguments['refdisp'] = self.action.args.refdisp
            arguments['subxvals'] = subxvals
            arguments['nn']= nn
            arguments['x0']= self.action.args.x0
            my_arguments.append(arguments)

        results = []
        centcoeff = []
        centwave = []
        centdisp = []

        p=Pool()
        results = p.map(myhelper, list(my_arguments))
        p.close()

        for result in results:
            b = result[0]
            scoeff = result[1]
            _centwave = result[2]
            _centdisp = result[3]
            centcoeff.append({b: scoeff})
            centwave.append(_centwave)
            centdisp.append(_centdisp)

        self.action.args.centcoeff = centcoeff

        #for b, bs in enumerate(self.context.arcs):
        if False:

            # wavelength coefficients
            coeff = [0., 0., 0., 0., 0.]
            # container for maxima, shifts
            maxima = []
            shifts = []
            # get sub spectrum for this bar
            subspec = bs[self.action.args.minrow:self.action.args.maxrow]
            # now loop over dispersions
            for di, disp in enumerate(disps):
                # populate the coefficients
                coeff[4] = p0[b]
                coeff[3] = disp
                cosbeta = disp / (self.context.config.instrument.PIX*ybin) * self.action.args.rho * \
                    self.context.config.instrument.FCAM * 1.e-4
                if cosbeta > 1.:
                    cosbeta = 1.
                beta = math.acos(cosbeta)
                coeff[2] = -(self.context.config.instrument.PIX * ybin / self.context.config.instrument.FCAM) ** 2 * \
                    math.sin(beta) / 2. / self.action.args.rho * 1.e4
                coeff[1] = -(self.context.config.instrument.PIX * ybin / self.context.config.instrument.FCAM) ** 3 * \
                    math.cos(beta) / 6. / self.action.args.rho * 1.e4
                coeff[0] = (self.context.config.instrument.PIX * ybin / self.context.config.instrument.FCAM) ** 4 * \
                    math.sin(beta) / 24. / self.action.args.rho * 1.e4
                # what are the min and max wavelengths to consider?
                wl0 = np.polyval(coeff, self.action.args.xvals[self.action.args.minrow])
                wl1 = np.polyval(coeff, self.action.args.xvals[self.action.args.maxrow])
                minwvl = np.nanmin([wl0, wl1])
                maxwvl = np.nanmax([wl0, wl1])
                # where will we need to interpolate to cross-correlate?
                minrw = [i for i, v in enumerate(self.action.args.refwave)
                         if v >= minwvl][0]
                maxrw = [i for i, v in enumerate(self.action.args.refwave)
                         if v <= maxwvl][-1]
                subrefwvl = self.action.args.refwave[minrw:maxrw]
                subrefspec = self.action.args.reflux[minrw:maxrw]
                # get bell cosine taper to avoid nasty edge effects
                tkwgt = signal.windows.tukey(len(subrefspec),
                                             alpha=self.context.config.instrument.TAPERFRAC)
                # apply taper to atlas spectrum
                subrefspec *= tkwgt
                # adjust wavelengths
                waves = np.polyval(coeff, subxvals)
                # interpolate the bar spectrum
                obsint = interpolate.interp1d(waves, subspec, kind='cubic',
                                              bounds_error=False,
                                              fill_value='extrapolate')
                intspec = obsint(subrefwvl)
                # apply taper to bar spectrum
                intspec *= tkwgt
                # get a label
                # cross correlate the interpolated spectrum with the atlas spec
                nsamp = len(subrefwvl)
                offar = np.arange(1 - nsamp, nsamp)
                # Cross-correlate
                xcorr = np.correlate(intspec, subrefspec, mode='full')
                # Get central region
                x0c = int(len(xcorr) / 3)
                x1c = int(2 * (len(xcorr) / 3))
                xcorr_central = xcorr[x0c:x1c]
                offar_central = offar[x0c:x1c]
                # Calculate offset
                maxima.append(xcorr_central[xcorr_central.argmax()])
                shifts.append(offar_central[xcorr_central.argmax()])
            # Get interpolations
            int_max = interpolate.interp1d(disps, maxima, kind='cubic',
                                           bounds_error=False,
                                           fill_value='extrapolate')
            int_shift = interpolate.interp1d(disps, shifts, kind='cubic',
                                             bounds_error=False,
                                             fill_value='extrapolate')
            xdisps = np.linspace(min(disps), max(disps), num=nn*100)
            # get peak values
            maxima_res = int_max(xdisps)
            shifts_res = int_shift(xdisps) * self.action.args.refdisp
            bardisp.append(xdisps[maxima_res.argmax()])
            barshift.append(shifts_res[maxima_res.argmax()])
            # update coeffs
            coeff[4] = p0[b] - barshift[-1]
            coeff[3] = bardisp[-1]
            cosbeta = coeff[3] / (self.context.config.instrument.PIX * ybin) * self.action.args.rho * \
                self.context.config.instrument.FCAM * 1.e-4
            if cosbeta > 1.:
                cosbeta = 1.
            beta = math.acos(cosbeta)
            coeff[2] = -(self.context.config.instrument.PIX * ybin / self.context.config.instrument.FCAM) ** 2 * \
                math.sin(beta) / 2. / self.action.args.rho * 1.e4
            coeff[1] = -(self.context.config.instrument.PIX * ybin / self.context.config.instrument.FCAM) ** 3 * \
                math.cos(beta) / 6. / self.action.args.rho * 1.e4
            coeff[0] = (self.context.config.instrument.PIX * ybin / self.context.config.instrument.FCAM) ** 4 * \
                math.sin(beta) / 24. / self.action.args.rho * 1.e4
            scoeff = pascal_shift(coeff, self.action.args.x0)
            self.logger.info("Central Fit: Bar#, Cdisp, Coefs: "
                          "%3d  %.4f  %.2f  %.4f  %13.5e %13.5e" %
                          (b, bardisp[-1], scoeff[4], scoeff[3], scoeff[2],
                           scoeff[1]))
            # store central values
            centwave.append(coeff[4])
            centdisp.append(coeff[3])
            # Store results
            self.action.args.centcoeff.append(coeff)

            if self.context.config.plot_level >= 1:
                # plot maxima
                p = figure(title="Bar %d, Slice %d" % (b, int(b/5)),
                           x_axis_label="Central dispersion (Ang/px)", y_axis_label="X-Corr Peak Value")

                p.scatter(disps, maxima, color='red')
                p.line(xdisps, int_max(xdisps))
                ylim_min = min(maxima)
                ylim_max = max(maxima)
                p.line([bardisp[-1], bardisp[-1]], [ylim_min, ylim_max], color='green')
                bokeh_plot(p)
                if do_inter:
                    q = input("<cr> - Next, q to quit: ")
                    if 'Q' in q.upper():
                        do_inter = False
                else:
                    time.sleep(0.01)

        if self.context.config.plot_level >= 1:
            # Plot results
            p = figure(title=imlab, x_axis_label="Bar #", y_axis_label="Central Wavelength (A)")
            x = range(len(centwave))
            p.scatter(x,centwave, marker='x')
            ylim = [min(centwave), max(centwave)]
            for ix in range(1, 24):
                 sx = ix*5 - 0.5
                 p.line([sx, sx], ylim, color='black', line_dash = 'dotted')
            p.x_range = Range1d(-1, 120)
            bokeh_plot(p)
            if self.context.config.plot_level >= 2:
                 input("Next? <cr>: ")
            else:
                 time.sleep(self.context.config.instrument.plot_pause)
            p = figure(title=imlab, x_axis_label="Bar #", y_axis_label="Central Dispersion (A)")
            x = range(len(centdisp))
            p.scatter(x, centdisp, marker='x')
            ylim = [min(centdisp), max(centdisp)]
            for ix in range(1, 24):
                sx = ix * 5 - 0.5
                p.line([sx, sx], ylim,  color='black', line_dash = 'dotted')
            p.x_range = Range1d(-1, 120)
            bokeh_plot(p)
            if self.context.config.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.context.config.instrument.plot_pause)

        #print(self.action.args.centcoeff)
        return self.action.args

class process_arc(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)

    def _perform(self):
        return self.action.args