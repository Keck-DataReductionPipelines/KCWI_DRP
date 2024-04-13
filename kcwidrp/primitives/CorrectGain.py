from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import parse_imsec, \
    kcwi_fits_writer


class CorrectGain(BasePrimitive):
    """
    Convert raw data numbers to electrons.

    Uses the ATSECn FITS header keywords to divide image into amp regions and
    then corrects each region with the corresponding GAINn keyword.  Updates the
    following FITS header keywords:

        * GAINCOR: sets to ``True`` if operation performed.
        * BUNIT: sets to `electron`.
        * HISTORY: records the operation.

    Uses the following configuration parameter:

        * saveintims: if set to ``True`` write out a \*_gain.fits file with gain corrected.  Default is ``False``.

    Updates image in returned arguments.

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        # Header keyword to update
        key = 'GAINCOR'
        keycom = 'Gain corrected?'
        number_of_amplifiers = self.action.args.namps
        bsec, dsec, tsec, direc, amps, aoff = self.action.args.map_ccd
        namps = len(amps)
        if namps != number_of_amplifiers:
            self.logger.warning("Amp count disagreement!")
        for amplifier in amps:
            # get amp section
            section = self.action.args.ccddata.header['ATSEC%d' %
                                                      amplifier]
            parsed_section, read_forward = parse_imsec(section)
            # get gain for this amp
            gain = self.action.args.ccddata.header['GAIN%d' % amplifier]
            self.logger.info(
                "Applying gain correction of %.3f in section %s" %
                (gain, self.action.args.ccddata.header['ATSEC%d' %
                                                       amplifier]))
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
