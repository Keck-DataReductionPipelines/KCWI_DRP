from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer
import os


class SubtractSky(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if a master sky exists to subtract
        :return:
        """
        self.logger.info("Checking precondition for CorrectIllumination")
        target_type = 'SKY'
        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type=target_type,
                                             nearest=True)
        self.logger.info("pre condition got %d master sky, expected 1"
                         % len(tab))
        if len(tab) <= 0:
            return False
        else:
            return True

    def _perform(self):
        self.logger.info("Subtracting sky background")

        # Header keyword to update
        key = 'SKYCOR'
        keycom = 'sky corrected?'
        target_type = 'SKY'

        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type=target_type,
                                             nearest=True)
        self.logger.info("%d master sky frames found" % len(tab))

        if len(tab) > 0:
            msname = tab['OFNAME'][0].split('.')[0] + '_' + \
                     target_type.lower() + ".fits"
            self.logger.info("Reading image: %s" % msname)
            msky = kcwi_fits_reader(
                os.path.join(os.path.dirname(self.action.args.name), 'redux',
                             msname))[0]

            # do the subtraction
            self.action.args.ccddata.data -= msky.data

            # update header keywords
            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['SKYMAST'] = (msname,
                                                          "Master sky filename")

        else:
            # update header keywords
            self.action.args.ccddata.header[key] = (False, keycom)

        log_string = SubtractSky.__module__ + "." + SubtractSky.__qualname__
        self.action.args.ccddata.header['HISTORY'] = log_string

        # write out int image
        kcwi_fits_writer(self.action.args.ccddata,
                         table=self.action.args.table,
                         output_file=self.action.args.name, suffix="intk")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="intk")
        self.context.proctab.write_proctab()

        self.logger.info(log_string)

        return self.action.args
    # END: class SubtractSky()
