from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments


class MakeInvsens(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger0

    def _perform(self):
        self.logger.info("Making inverse sensitivity curve "
                         "(not yet implemented)")

        log_string = MakeInvsens.__module__ + "." + MakeInvsens.__qualname__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class MakeInvsens()

