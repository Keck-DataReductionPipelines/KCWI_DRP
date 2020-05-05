from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_correct_extin import kcwi_correct_extin
from bokeh.plotting import figure
from scipy.signal import find_peaks

import os
import pkg_resources
import numpy as np


class MakeInvsens(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """Checks if object is a standard star"""
        self.logger.info("Checking precondition for MakeInvsens")

        stdfile = None
        if 'object' in self.action.args.imtype.lower():
            obname = self.action.args.ccddata.header['OBJECT'].lower()
            path = 'data/stds/%s.fits' % obname
            package = __name__.split('.')[0]
            full_path = pkg_resources.resource_filename(package, path)
            if os.path.exists(full_path):
                self.logger.info("Found std file: %s" % full_path)
                stdfile = full_path
            else:
                self.logger.info("Not found in data/stds: %s" % full_path)
        else:
            self.logger.warning("Not object type: %s" %
                                self.action.args.imtype)
        self.action.args.stdfile = stdfile
        if self.action.args.stdfile is not None:
            return True
        else:
            return False

    def _perform(self):
        self.logger.info("Making inverse sensitivity curve")

        suffix = 'invsens'

        # get size
        sz = self.action.args.ccddata.data.shape
        print(sz)
        # default pixel ranges
        z = np.arange(sz[0])
        z0 = 175
        z1 = sz[0] - 175
        # get exposure time
        expt = self.action.args.ccddata.header['XPOSURE']
        if expt == 0.:
            self.logger.warning("No exposure time found, setting to 1s")
            expt = 1.
        else:
            self.logger.info("Using exposure time of %.1f" % expt)
        # get wavelength scale
        w0 = self.action.args.ccddata.header['CRVAL3']
        dw = self.action.args.ccddata.header['CD3_3']
        crpixw = self.action.args.ccddata.header['CRPIX3']
        # get all good wavelength range
        wgoo0 = self.action.args.ccddata.header['WAVGOOD0']
        if wgoo0 < 3650:
            wgoo0 = 3650.
        wgoo1 = self.action.args.ccddata.header['WAVGOOD1']
        # get all inclusive wavelength range
        wall0 = self.action.args.ccddata.header['WAVALL0']
        wall1 = self.action.args.ccddata.header['WAVALL1']
        # get DAR padding in y
        pad_y = self.action.args.ccddata.header['DARPADY']
        # get sky subtraction status
        skycor = self.action.args.ccddata.header['SKYCOR']
        # get telescope and atm. correction
        tel = self.action.args.ccddata.header['TELESCOP']
        if 'Keck' in tel:
            area = 760000.0
        else:
            area = -1.0
        tlab = tel
        # compute good y pixel ranges
        if w0 > 0. and dw > 0. and wgoo0 > 0. and wgoo1 > 0.:
            z0 = int((wgoo0 - w0) / dw) + 10
            z1 = int((wgoo1 - w0) / dw) - 10
        gz = [i for i, v in enumerate(z) if z0 <= v <= z1]
        # wavelength scale
        w = w0 + z * dw
        # good spatial range
        gy0 = pad_y if pad_y > 1 else 1
        gy1 = sz[1] - (pad_y if pad_y > 2 else 2)
        # log results
        self.logger.info("Invsen. Pars: Y0, Y1, Z0, Z1, Wav0, Wav1: "
                         "%d, %d, %d, %d, %.3f, %.3f" %
                         (gy0, gy1, z0, z1, w[z0], w[z1]))
        # central wavelength
        cwv = self.action.args.cwave
        # find standard in slices
        # sum over wavelength
        tot = np.sum(self.action.args.ccddata.data[z0:z1, gy0:gy1, :], 0)
        yy = np.arange(gy1-gy0) + gy0
        mxsl = -1.
        mxsg = 0.
        # for each slice
        for i in range(sz[2]):
            tstd = float(np.nanstd(tot[:, i]))
            if tstd > mxsg:
                mxsg = tstd
                mxsl = i

        # relevant slices
        sl0 = (mxsl - 3) if mxsl >= 3 else 0
        sl1 = (mxsl + 3) if (mxsl + 3) <= sz[2]-1 else sz[2]-1
        # get y position of std
        cy, _ = find_peaks(tot[:, mxsl], height=np.nanmean(tot[:, mxsl]))
        cy = int(cy[0]) + gy0
        # log results
        self.logger.info("Std lices: max, sl0, sl1, spatial cntrd: "
                         "%d, %d, %d, %.2f" % (mxsl, sl0, sl1, cy))
        # get dwave spectrum
        dwspec = np.zeros(sz[0]) + dw
        # copy of input cube
        scub = self.action.args.ccddata.data.copy()
        # sky window width in pixels
        skywin = int(self.config.instrument.psfwid / self.action.args.xbinsize)
        # do sky subtraction, if needed
        if not skycor:
            self.logger.warning("Sky should have been subtraced already")
        # apply extinction correction
        ucub = scub.copy()  # uncorrected cube
        kcwi_correct_extin(scub, self.action.args.ccddata.header,
                           logger=self.logger)

        # get slice spectra y limits
        sy0 = (cy - skywin) if (cy - skywin) > 0 else 0
        sy1 = (cy + skywin) if (cy + skywin) < (sz[1]-1) else (sz[1]-1)
        # sum over y range
        slspec = np.sum(scub[:, sy0:sy1, :], 1)
        ulspec = np.sum(ucub[:, sy0:sy1, :], 1)
        # sum over slices
        obsspec = np.sum(slspec[:, sl0:sl1], 1)
        ubsspec = np.sum(ulspec[:, sl0:sl1], 1)
        # convert to e-/second
        obsspec /= expt
        ubsspec /= expt

        p = figure(
            title=self.action.args.plotlabel + ' Sum',
            x_axis_label='Wave (A)',
            y_axis_label='Flux (e-/s)',
            plot_width=self.config.instrument.plot_width,
            plot_height=self.config.instrument.plot_height)
        p.line(w, ubsspec, line_color='red', legend_label='Raw')
        p.line(w, obsspec, line_color='green', legend_label='ExtCor')
        bokeh_plot(p, self.context.bokeh_session)
        input("Next? <cr>: ")

        log_string = MakeInvsens.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class MakeInvsens()
