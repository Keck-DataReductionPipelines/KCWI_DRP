from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer

import numpy as np
from astroscrappy import detect_cosmics


class RemoveCosmicRays(BasePrimitive):
    """
    Remove cosmic rays and generate a flag image recording their location.

    Uses astroscrappy to detect and flag cosmic rays.  Updates the following
    FITS header keywords:

        * CRCLEAN: set to ``True`` if operation performed.
        * NCRCLEAN: set to the number of cosmic ray pixels cleaned.

    Uses the following configuration parameters:

        * saveintims: if set to ``True`` write out a CR cleaned version of image in \*_crr.fits.  Defaults to ``False``.
        * CRR_MINEXPTIME - exposure time below which no CR cleaning is done.
        * CRR\_\* - see kcwi.cfg file for parameters controlling CR cleaning.

    Updates image in returned arguments with CR cleaned image and adds flags
    indicating where changes were made.

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        # TODO: implement parameter options from kcwi_stage1.pro
        self.logger.info("Cleaning and flagging cosmic rays")

        # Header keyword to update
        key = 'CRCLEAN'
        keycom = 'cosmic rays cleaned?'

        header = self.action.args.ccddata.header

        exptime = header['TELAPSE']
        nshuf = header['NSHFUP']
        ttime = header['TTIME']
        if nshuf * ttime > exptime:
            exptime = nshuf * ttime

        if exptime >= self.config.instrument.CRR_MINEXPTIME:

            namps = header['NVIDINP']
            bsec, dsec, tsec, direc, amps, aoff = self.action.args.map_ccd
            read_noise = 0.
            if len(amps) != namps:
                self.logger.warning("Amp count disagreement!")
            for ia in amps:
                if 'BIASRN%d' % ia in header:
                    read_noise += header['BIASRN%d' % ia]
                elif 'OSCNRN%d' % ia in header:
                    read_noise += header['OSCNRN%d' % ia]
                else:
                    read_noise += 3.
            read_noise /= float(namps)

            # Set sigclip according to image parameters
            sigclip = self.config.instrument.CRR_SIGCLIP
            if 'FLATLAMP' in self.action.args.ccddata.header['IMTYPE']:
                if self.action.args.nasmask:
                    sigclip = 10.
                else:
                    sigclip = 7.
            if 'OBJECT' in self.action.args.ccddata.header['IMTYPE']:
                if self.action.args.ccddata.header['TTIME'] < 300.:
                    sigclip = 10.

            mask, clean = detect_cosmics(
                self.action.args.ccddata.data, gain=1.0, readnoise=read_noise,
                psffwhm=self.config.instrument.CRR_PSFFWHM,
                sigclip=sigclip,
                sigfrac=self.config.instrument.CRR_SIGFRAC,
                objlim=self.config.instrument.CRR_OBJLIM,
                fsmode=self.config.instrument.CRR_FSMODE,
                psfmodel=self.config.instrument.CRR_PSFMODEL,
                verbose=self.config.instrument.CRR_VERBOSE,
                sepmed=self.config.instrument.CRR_SEPMED,
                cleantype=self.config.instrument.CRR_CLEANTYPE)

            self.logger.info("Astroscrappy: cleaned cosmic rays")
            header['history'] = "Astroscrappy: cleaned cosmic rays"
            header['history'] = \
                "Astroscrappy params: sigclip=%5.2f sigfrac=%5.2f " \
                "objlim=%5.2f" % (
                self.config.instrument.CRR_SIGCLIP,
                self.config.instrument.CRR_SIGFRAC,
                self.config.instrument.CRR_OBJLIM)
            header['history'] = \
                "Astroscrappy params: fsmode=%s psfmodel=%s psffwhm=%5.2f" % (
                self.config.instrument.CRR_FSMODE,
                self.config.instrument.CRR_PSFMODEL,
                self.config.instrument.CRR_PSFFWHM)
            header['history'] = "Astroscrappy params: sepmed=%s " \
                                "minexptime=%f" % (
                self.config.instrument.CRR_SEPMED,
                self.config.instrument.CRR_MINEXPTIME)
            # header['history'] = "LA CosmicX run on %s" % time.strftime("%c")

            # update arrays
            mask = np.cast["bool"](mask)
            fmask = np.where(mask)
            try:
                self.action.args.ccddata.flags[fmask] += 4
            except AttributeError:
                self.logger.warning("Flags array not found!")
            n_crs = mask.sum()
            # DN 2023-may-28: commenting out because it causes bad things
            # self.action.args.ccddata.mask += mask
            self.action.args.ccddata.data = clean
            # update header
            header[key] = (True, keycom)
            header['NCRCLEAN'] = (n_crs, "number of cosmic ray pixels")

        else:
            self.logger.info("Astroscrappy: exptime < minexptime=%.1f" %
                             self.config.instrument.CRR_MINEXPTIME)
            header['history'] = \
                "Astroscrappy: exptime < minexptime=%.1f" % \
                self.config.instrument.CRR_MINEXPTIME
            header[key] = (False, keycom)
            header['NCRCLEAN'] = (0, "number of cosmic ray pixels")

        log_string = RemoveCosmicRays.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        if self.config.instrument.saveintims:
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name,
                             output_dir=self.config.instrument.output_directory,
                             suffix="crr")

        return self.action.args
    # END: class RemoveCosmicRays()
