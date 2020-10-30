from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import get_plot_lims, oplot_slices, \
    set_plot_lims

from bokeh.plotting import figure
import numpy as np
import math
from scipy.interpolate import interpolate
from multiprocessing import get_context
from scipy import signal
import time


def pascal_shift(coefficients=None, x0=None):
    """Shift coefficients to a new reference value (X0)

    This should probably go somewhere else, but will be needed here.
    """
    if not coefficients:
        print("Error, no coefficients for pascal_shift.")
        return None
    if not x0:
        print("Warning, no reference value (x0) supplied")
        return coefficients
    if len(coefficients) == 7:
        usecoeff = list(reversed(coefficients))
        fincoeff = [0.] * 7
    else:
        if len(coefficients) > 7:
            print("Warning - this routine only handles up to 7 coefficients.")
            usecoeff = list(reversed(coefficients[0:7]))
            fincoeff = [0.] * len(coefficients)
        else:
            usecoeff = [0.] * 7
            fincoeff = usecoeff
            for ic, c in enumerate(coefficients):
                usecoeff[len(coefficients) - (ic + 1)] = coefficients[ic]
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
    if len(coefficients) < 7:
        fincoeff = fincoeff[0:len(coefficients)]
    # Reverse for python
    return list(reversed(fincoeff))
    # END: def pascal_shift()


def bar_fit_helper(argument):

    b = argument['b']
    bs = argument['bs']
    # wavelength coefficients
    coefficients = [0., 0., 0., 0., 0.]
    # container for maxima, shifts
    maxima = []
    shifts = []
    # get sub spectrum for this bar
    sub_spectrum = bs[argument['minrow']:argument['maxrow']]
    # now loop over dispersions
    for di, dispersion in enumerate(argument['disps']):
        # populate the coefficients
        coefficients[4] = argument['p0'][b]
        coefficients[3] = dispersion
        cosbeta = dispersion / (argument['PIX'] * argument['ybin']) * \
            argument['rho'] * argument['FCAM'] * 1.e-4
        if cosbeta > 1.:
            cosbeta = 1.
        beta = math.acos(cosbeta)
        coefficients[2] = -(argument['PIX'] * argument['ybin'] /
                            argument['FCAM']) ** 2 * math.sin(beta) / 2. / \
            argument['rho'] * 1.e4
        coefficients[1] = -(argument['PIX'] * argument['ybin'] /
                            argument['FCAM']) ** 3 * math.cos(beta) / 6. / \
            argument['rho'] * 1.e4
        coefficients[0] = (argument['PIX'] * argument['ybin'] /
                           argument['FCAM']) ** 4 * math.sin(beta) / 24. / \
            argument['rho'] * 1.e4
        # what are the min and max wavelengths to consider?
        wl0 = np.polyval(coefficients, argument['xvals'][argument['minrow']])
        wl1 = np.polyval(coefficients, argument['xvals'][argument['maxrow']])
        minimum_wavelength = np.nanmin([wl0, wl1])
        maximum_wavelength = np.nanmax([wl0, wl1])
        # where will we need to interpolate to cross-correlate?
        minrw = [i for i, v in enumerate(argument['refwave'])
                 if v >= minimum_wavelength][0]
        maxrw = [i for i, v in enumerate(argument['refwave'])
                 if v <= maximum_wavelength][-1]
        ref_wave_of_sub_spectrum = argument['refwave'][minrw:maxrw]
        ref_flux_of_sub_spectrum = argument['reflux'][minrw:maxrw]
        # get bell cosine taper to avoid nasty edge effects
        tkwgt = signal.windows.tukey(len(ref_flux_of_sub_spectrum),
                                     alpha=argument['taperfrac'])
        # apply taper to atlas spectrum
        ref_flux_of_sub_spectrum *= tkwgt
        # adjust wavelengths
        waves = np.polyval(coefficients, argument['subxvals'])
        # interpolate the bar spectrum
        obsint = interpolate.interp1d(waves, sub_spectrum, kind='cubic',
                                      bounds_error=False,
                                      fill_value='extrapolate')
        intspec = obsint(ref_wave_of_sub_spectrum)
        # apply taper to bar spectrum
        intspec *= tkwgt
        # get a label
        # cross correlate the interpolated spectrum with the atlas spec
        samples_number = len(ref_wave_of_sub_spectrum)
        offsets_array = np.arange(1 - samples_number, samples_number)
        # Cross-correlate
        crosscorrelation = np.correlate(intspec, ref_flux_of_sub_spectrum,
                                        mode='full')
        # Get central region
        x0c = int(len(crosscorrelation) / 3)
        x1c = int(2 * (len(crosscorrelation) / 3))
        central_crosscorrelation = crosscorrelation[x0c:x1c]
        central_offsets_array = offsets_array[x0c:x1c]
        # Calculate offset
        maxima.append(central_crosscorrelation[
                          central_crosscorrelation.argmax()])
        shifts.append(central_offsets_array[central_crosscorrelation.argmax()])
    # Get interpolations
    int_max = interpolate.interp1d(argument['disps'], maxima, kind='cubic',
                                   bounds_error=False,
                                   fill_value='extrapolate')
    int_shift = interpolate.interp1d(argument['disps'], shifts, kind='cubic',
                                     bounds_error=False,
                                     fill_value='extrapolate')
    xdisps = np.linspace(min(argument['disps']), max(argument['disps']),
                         num=argument['nn'] * 100)
    # get central region
    x0c = int(len(xdisps) / 3)
    x1c = int(2 * (len(xdisps) / 3))
    # get peak values
    central_xdisps = xdisps[x0c:x1c]
    maxima_res = int_max(central_xdisps)
    shifts_res = int_shift(central_xdisps) * argument['refdisp']
    bardisp = central_xdisps[maxima_res.argmax()]
    barshift = shifts_res[maxima_res.argmax()]
    # update coeffs
    coefficients[4] = argument['p0'][b] - barshift
    coefficients[3] = bardisp
    cosbeta = coefficients[3] / (argument['PIX'] * argument['ybin']) * \
        argument['rho'] * argument['FCAM'] * 1.e-4
    if cosbeta > 1.:
        cosbeta = 1.
    beta = math.acos(cosbeta)
    coefficients[2] = -(argument['PIX'] * argument['ybin'] /
                        argument['FCAM']) ** 2 * \
        math.sin(beta) / 2. / argument['rho'] * 1.e4
    coefficients[1] = -(argument['PIX'] * argument['ybin'] /
                        argument['FCAM']) ** 3 * \
        math.cos(beta) / 6. / argument['rho'] * 1.e4
    coefficients[0] = (argument['PIX'] * argument['ybin'] /
                       argument['FCAM']) ** 4 * \
        math.sin(beta) / 24. / argument['rho'] * 1.e4
    shifted_coefficients = pascal_shift(coefficients, argument['x0'])
    print("Bar#: %3d, Cdisp: %.4f" % (b, bardisp))

    # Return results
    return b, shifted_coefficients, coefficients[4], coefficients[3], \
           maxima, bardisp
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
        y_binning = self.action.args.ybinsize
        # let's populate the 0 points vector
        p0 = self.action.args.cwave + np.array(self.context.bar_offsets) * \
            self.context.prelim_disp - self.action.args.offset_wave
        # next we are going to brute-force scan around the preliminary
        # dispersion for a better solution. We will wander 5% away from it.
        maximum_dispersion_deviation = 0.05  # fraction
        # we will try nn values
        self.logger.info("prelim disp = %.3f, refdisp = %.3f,"
                         " min,max rows = %d, %d" % (self.context.prelim_disp,
                                                     self.action.args.refdisp,
                                                     self.action.args.minrow,
                                                     self.action.args.maxrow))
        number_of_values_to_try = (int(maximum_dispersion_deviation *
                                       abs(self.context.prelim_disp) /
                                   self.action.args.refdisp *
                                       (self.action.args.maxrow -
                                        self.action.args.minrow) / 2.0))
        if number_of_values_to_try < 10:
            number_of_values_to_try = 10
        if number_of_values_to_try > 50:
            number_of_values_to_try = 50
        self.logger.info("N disp. samples: %d" % number_of_values_to_try)
        # dispersions to try
        disps = self.context.prelim_disp * (
                1.0 + maximum_dispersion_deviation *
                (np.arange(0, number_of_values_to_try + 1) -
                 number_of_values_to_try / 2.) *
                2.0 / number_of_values_to_try)
        # values for central fit
        subxvals = self.action.args.xvals[
                   self.action.args.minrow:self.action.args.maxrow]

        # log taperfrac: important!
        self.logger.info("Using TAPERFRAC = %.3f" %
                         self.config.instrument.TAPERFRAC)
        self.action.args.ccddata.header['TAPFRAC'] = (
            self.config.instrument.TAPERFRAC, "taper fraction for central fit")
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
                'ybin': y_binning,
                'rho': self.action.args.rho,
                'FCAM': self.config.instrument.FCAM,
                'xvals': self.action.args.xvals,
                'refwave': self.action.args.refwave,
                'reflux': self.action.args.reflux,
                'taperfrac': self.config.instrument.TAPERFRAC,
                'refdisp': self.action.args.refdisp,
                'subxvals': subxvals,
                'nn': number_of_values_to_try, 'x0': self.action.args.x0
            }
            my_arguments.append(arguments)

        twkcoeff = {}
        centwave = []
        centdisp = []

        p = get_context("spawn").Pool()
        results = p.map(bar_fit_helper, list(my_arguments))
        p.close()

        next_bar_to_plot = 0
        for ir, result in enumerate(results):
            b = result[0]
            shifted_coefficients = result[1]
            _centwave = result[2]
            _centdisp = result[3]
            twkcoeff[b] = shifted_coefficients
            centwave.append(_centwave)
            centdisp.append(_centdisp)
            maxima = result[4]
            bardisp = result[5]
            self.logger.info("Central Fit: Bar# %3d, Cdisp %.4f, "
                             "Coefs: %.2f %.4f %13.5e %13.5e" %
                             (b, bardisp, shifted_coefficients[4],
                              shifted_coefficients[3], shifted_coefficients[2],
                              shifted_coefficients[1]))
            if do_inter and ir == next_bar_to_plot:
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
                q = input("Next? <int> or <cr>, q to quit: ")
                if 'Q' in q.upper():
                    do_inter = False
                else:
                    try:
                        next_bar_to_plot = int(q)
                    except ValueError:
                        next_bar_to_plot = ir + 1

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
            xlim = [-1, 120]
            ylim = get_plot_lims(centwave)
            p.xgrid.grid_line_color = None
            oplot_slices(p, ylim)
            p.legend.location = "top_center"
            set_plot_lims(p, xlim=xlim, ylim=ylim)
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
            xlim = [-2, 121]
            ylim = get_plot_lims(centdisp)
            p.xgrid.grid_line_color = None
            oplot_slices(p, ylim)
            p.legend.location = "bottom_center"
            set_plot_lims(p, xlim=xlim, ylim=ylim)
            bokeh_plot(p, self.context.bokeh_session)
            if self.config.instrument.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.config.instrument.plot_pause)

        log_string = FitCenter.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class FitCenter()
