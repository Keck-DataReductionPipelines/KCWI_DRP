from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
        get_master_name
import os


class SubtractDark(BasePrimitive):
    """
    Subtract the master dark frame.

    Checks for existence of master dark frame and, if it exists, performs
    subtraction and records processing in the header.

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):

        # Header keyword to update
        key = 'DARKSUB'
        keycom = 'master dark subtracted?'
        target_type = 'MDARK'

        self.logger.info("Subtracting master dark")
        tab = self.context.proctab.search_proctab(
            frame=self.action.args.ccddata, target_type=target_type,
            nearest=True)
        self.logger.info("%d master dark frames found" % len(tab))

        if len(tab) > 0:
            mdname = get_master_name(tab, target_type)
            self.logger.info("*************** READING IMAGE: %s" % mdname)
            mdark = kcwi_fits_reader(
                os.path.join(self.config.instrument.cwd, 'redux',
                             mdname))[0]
            # scale by exposure time
            fac = 1.0
            if 'TTIME' in mdark.header and \
               'TTIME' in self.action.args.ccddata.header:
                fac = float(self.action.args.ccddata.header['TTIME']) / \
                      float(mdark.header['TTIME'])
                self.logger.info("dark scaled by %.3f" % fac)
            else:
                self.logger.warn("unable to scale dark by exposure time")

            # do the subtraction
            self.action.args.ccddata.data -= mdark.data * fac

            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['MDFILE'] = (mdname,
                                                         "Master dark filename")
            self.action.args.ccddata.header['DARKSCL'] = (fac,
                                                          "dark scale factor")
        else:
            self.logger.info("No master dark frame available, skipping")
            self.action.args.ccddata.header[key] = (False, keycom)

        log_string = SubtractDark.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class SubtractDark()
