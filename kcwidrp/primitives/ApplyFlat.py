from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments


class ApplyFlat(BasePrimitive):
    """
    Generic routine to apply flat fielding. Not yet implemented.
    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Applying flat field (not yet implemented)")

        log_string = ApplyFlat.__module__ + "." + ApplyFlat.__qualname__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class ApplyFlat()

