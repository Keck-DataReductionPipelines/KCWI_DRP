from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from kcwidrp.primitives.kcwi_file_primitives import parse_imsec


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
        namps = self.action.args.namps
        for ia in range(namps):
            # get amp section
            section = self.action.args.ccddata.header['ATSEC%d' % (ia + 1)]
            sec, rfor = parse_imsec(section)
            # get gain for this amp
            gain = self.context.data_set.get_info_column(
                self.action.args.name, 'GAIN%d' % (ia + 1))
            self.logger.info(
                "Applying gain correction of %.3f in section %s" %
                (gain, self.action.args.ccddata.header['ATSEC%d' % (ia + 1)]))
            self.action.args.ccddata.data[sec[0]:(sec[1]+1),
                                          sec[2]:(sec[3]+1)] *= gain

        self.action.args.ccddata.header[key] = (True, keycom)
        self.action.args.ccddata.header['BUNIT'] = ('electron', 'Pixel units')
        self.action.args.ccddata.unit = 'electron'

        logstr = CorrectGain.__module__ + "." + CorrectGain.__qualname__
        self.action.args.ccddata.header['HISTORY'] = logstr
        self.logger.info(logstr)

        if self.config.instrument.saveintims:
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name, suffix="gain")
        return self.action.args
    # END: class CorrectGain()
