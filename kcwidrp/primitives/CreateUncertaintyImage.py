from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from kcwidrp.primitives.kcwi_file_primitives import parse_imsec

from astropy.nddata import VarianceUncertainty


class CreateUncertaintyImage(BasePrimitive):
    """Generate a variance image based on Poisson noise plus readnoise"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        """Assumes units of image are electron"""

        # Header keyword to update
        key = 'UNCVAR'
        keycom = 'variance created?'

        self.logger.info("Create uncertainty image")
        # start with Poisson noise
        self.action.args.ccddata.uncertainty = VarianceUncertainty(
            self.action.args.ccddata.data, unit='electron^2', copy=True)
        # add readnoise, if known
        if 'BIASRN1' in self.action.args.ccddata.header:
            namps = self.action.args.ccddata.header['NVIDINP']
            for ia in range(namps):
                # get amp parameters
                biasrn = self.action.args.ccddata.header['BIASRN%d' % (ia + 1)]
                section = self.action.args.ccddata.header['ATSEC%d' % (ia + 1)]
                sec, rfor = parse_imsec(section)
                self.action.args.ccddata.uncertainty.array[
                    sec[0]:(sec[1]+1), sec[2]:(sec[3]+1)] += biasrn
        else:
            self.logger.warn("Readnoise undefined, uncertainty Poisson only")
        # document variance image creation
        self.action.args.ccddata.header[key] = (True, keycom)

        logstr = CreateUncertaintyImage.__module__ + \
            "." + CreateUncertaintyImage.__qualname__
        self.action.args.ccddata.header['HISTORY'] = logstr
        self.logger.info(logstr)

        return self.action.args
    # END: class CreateUncertaintyImage()


