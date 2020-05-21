from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import parse_imsec, \
    kcwi_fits_writer


class CorrectGain(BasePrimitive):
    """Convert raw data numbers to electrons"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        # Header keyword to update
        key = 'GAINCOR'
        keycom = 'Gain corrected?'
        # print(self.action.args.ccddata.header)
        number_of_amplifiers = self.action.args.namps
        for amplifier in range(number_of_amplifiers):
            # get amp section
            section = self.action.args.ccddata.header['ATSEC%d' %
                                                      (amplifier + 1)]
            parsed_section, read_forward = parse_imsec(section)
            # get gain for this amp
            gain = self.action.args.ccddata.header['GAIN%d' % (amplifier + 1)]
            self.logger.info(
                "Applying gain correction of %.3f in section %s" %
                (gain, self.action.args.ccddata.header['ATSEC%d' %
                                                       (amplifier + 1)]))
            self.action.args.ccddata.data[
                parsed_section[0]:(parsed_section[1]+1),
                parsed_section[2]:(parsed_section[3]+1)] *= gain

        self.action.args.ccddata.header[key] = (True, keycom)
        self.action.args.ccddata.header['BUNIT'] = ('electron', 'Pixel units')
        self.action.args.ccddata.unit = 'electron'

        log_string = CorrectGain.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        if self.config.instrument.saveintims:
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name,
                             output_dir=self.config.instrument.output_directory,
                             suffix="gain")
        return self.action.args
    # END: class CorrectGain()
