from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from kcwidrp.core.bokeh_plotting import bokeh_plot, bokeh_save

from bokeh.plotting import figure, show
from bokeh.models import Range1d, LinearAxis
import numpy as np
import math
from scipy.interpolate import interpolate
from multiprocessing import Pool
from scipy import signal
from scipy.stats import sigmaclip
import time


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
                usecoeff[len(coef) - (ic + 1)] = coef[ic]
    # get reference values
    x01 = x0
    x02 = x0 ** 2
    x03 = x0 ** 3
    x04 = x0 ** 4
    x05 = x0 ** 5
    x06 = x0 ** 6
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
    # END: def pascal_shift()


def bar_fit_helper(argument):
    logger = argument['logger']

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
        cosbeta = disp / (argument['PIX'] * argument['ybin']) * \
                  argument['rho'] * argument['FCAM'] * 1.e-4
        if cosbeta > 1.:
            cosbeta = 1.
        beta = math.acos(cosbeta)
        coeff[2] = -(argument['PIX'] * argument['ybin'] /
                     argument['FCAM']) ** 2 * math.sin(beta) / 2. / \
                   argument['rho'] * 1.e4
        coeff[1] = -(argument['PIX'] * argument['ybin'] /
                     argument['FCAM']) ** 3 * math.cos(beta) / 6. / \
                   argument['rho'] * 1.e4
        coeff[0] = (argument['PIX'] * argument['ybin'] /
                    argument['FCAM']) ** 4 * math.sin(beta) / 24. / \
                   argument['rho'] * 1.e4
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
    xdisps = np.linspace(min(argument['disps']), max(argument['disps']),
                         num=argument['nn'] * 100)
    # get peak values
    maxima_res = int_max(xdisps)
    shifts_res = int_shift(xdisps) * argument['refdisp']
    bardisp = xdisps[maxima_res.argmax()]
    barshift = shifts_res[maxima_res.argmax()]
    # update coeffs
    coeff[4] = argument['p0'][b] - barshift
    coeff[3] = bardisp
    cosbeta = coeff[3] / (argument['PIX'] * argument['ybin']) * \
              argument['rho'] * argument['FCAM'] * 1.e-4
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
    logger.info("Central Fit: Bar# %3d, Cdisp %.4f, Coefs: %.2f  %.4f  %13.5e %13.5e" %
                (b, bardisp, scoeff[4], scoeff[3], scoeff[2], scoeff[1]))
    # Return results
    return b, scoeff, coeff[4], coeff[3], maxima
    # END: def bar_fit_helper()


class FitCenter(BasePrimitive):
    """ Fit central region"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        self.action.args.twkcoeff = []

    def _perform(self):
        """At this point we have the offsets between bars and the approximate
        offset from the reference bar to the atlas spectrum and the approximate
        dispersion.
        """
        self.logger.info("Finding wavelength solution for central region")
        # Are we interactive?
        do_inter = (self.config.instrument.plot_level >= 2)

        # y binning
        ybin = self.action.args.ybinsize
        # let's populate the 0 points vector
        p0 = self.action.args.cwave + np.array(self.context.baroffs) * \
             self.context.prelim_disp - self.action.args.offset_wave
        # next we are going to brute-force scan around the preliminary
        # dispersion for a better solution. We will wander 5% away from it.
        max_ddisp = 0.05  # fraction
        # we will try nn values
        nn = (int(max_ddisp * abs(self.context.prelim_disp) /
                  self.action.args.refdisp * (self.action.args.maxrow -
                                              self.action.args.minrow) / 3.0))
        if nn < 10:
            nn = 10
        if nn > 25:
            nn = 25
        self.logger.info("N disp. samples: %d" % nn)
        # dispersions to try
        disps = self.context.prelim_disp * (1.0 + max_ddisp *
                                            (np.arange(0, nn + 1) - nn / 2.) *
                                            2.0 / nn)
        # values for central fit
        subxvals = self.action.args.xvals[
                   self.action.args.minrow:self.action.args.maxrow]

        # loop over bars and assemble input arguments
        my_arguments = []
        for b, bs in enumerate(self.context.arcs):
            arguments = {
                'b': b,
                'bs': bs,
                'minrow': self.action.args.minrow,
                'maxrow': self.action.args.maxrow,
                'disps': disps,
                'p0': p0,
                'PIX': self.config.instrument.PIX,
                'ybin': ybin,
                'rho': self.action.args.rho,
                'FCAM': self.config.instrument.FCAM,
                'xvals': self.action.args.xvals,
                'refwave': self.action.args.refwave,
                'reflux': self.action.args.reflux,
                'taperfrac': self.config.instrument.TAPERFRAC,
                'refdisp': self.action.args.refdisp,
                'subxvals': subxvals,
                'nn': nn, 'x0': self.action.args.x0,
                'logger': self.logger,
            }
            my_arguments.append(arguments)

        twkcoeff = {}
        centwave = []
        centdisp = []

        p = Pool()
        results = p.map(bar_fit_helper, list(my_arguments))
        p.close()

        for result in results:
            b = result[0]
            scoeff = result[1]
            _centwave = result[2]
            _centdisp = result[3]
            twkcoeff[b] = scoeff
            centwave.append(_centwave)
            centdisp.append(_centdisp)
            maxima = result[4]
            if do_inter:
                # plot maxima
                p = figure(title=self.action.args.plotlabel +
                                 "CENTRAL DISPERSION FIT for Bar: %d Slice: %d" %
                                 (b, int(b / 5)),
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height,
                           x_axis_label="Central dispersion (Ang/px)",
                           y_axis_label="X-Corr Peak Value")
                p.scatter(disps, maxima, color='red', legend_label="Data")
                p.line(disps, maxima, color='blue', legend_label="Data")
                ylim = [min(maxima), max(maxima)]
                p.line([_centdisp, _centdisp], ylim, color='green',
                       legend_label="Fit Disp")
                p.line([self.context.prelim_disp, self.context.prelim_disp],
                       ylim, color='red', legend_label="Calc Disp")
                bokeh_plot(p, self.context.bokeh_session)
                q = input("Next? <cr>, q to quit: ")
                if 'Q' in q.upper():
                    do_inter = False

        self.action.args.twkcoeff = twkcoeff
        # Plot results
        if self.config.instrument.plot_level >= 1:
            # Plot central wavelength
            p = figure(title=self.action.args.plotlabel + "CENTRAL VALUES",
                       x_axis_label="Bar #",
                       y_axis_label="Central Wavelength (A)",
                       plot_width=self.config.instrument.plot_width,
                       plot_height=self.config.instrument.plot_height)
            x = range(len(centwave))
            p.scatter(x, centwave, marker='x', legend_label='bar wave')
            p.line([0, 120], [self.action.args.cwave, self.action.args.cwave],
                   color='red', legend_label='CWAVE')
            ylim = [min(centwave), max(centwave)]
            for ix in range(1, 24):
                sx = ix * 5 - 0.5
                p.line([sx, sx], ylim, color='black')
            p.x_range = Range1d(-1, 120)
            p.legend.location = "top_center"
            bokeh_plot(p, self.context.bokeh_session)
            if self.config.instrument.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.config.instrument.plot_pause)
            # Plot central dispersion
            p = figure(title=self.action.args.plotlabel + "CENTRAL VALUES",
                       x_axis_label="Bar #",
                       y_axis_label="Central Dispersion (A)",
                       plot_width=self.config.instrument.plot_width,
                       plot_height=self.config.instrument.plot_height)
            x = range(len(centdisp))
            p.scatter(x, centdisp, marker='x', legend_label='bar disp')
            p.line([0, 120], [self.context.prelim_disp,
                              self.context.prelim_disp], color='red',
                   legend_label='Calc Disp')
            ylim = [min(centdisp), max(centdisp)]
            for ix in range(1, 24):
                sx = ix * 5 - 0.5
                p.line([sx, sx], ylim, color='black')
            p.x_range = Range1d(-1, 120)
            p.legend.location = "bottom_center"
            bokeh_plot(p, self.context.bokeh_session)
            if self.config.instrument.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.config.instrument.plot_pause)

        logstr = FitCenter.__module__ + "." + FitCenter.__qualname__
        self.action.args.ccddata.header['HISTORY'] = logstr
        self.logger.info(logstr)

        return self.action.args
    # END: class FitCenter()
