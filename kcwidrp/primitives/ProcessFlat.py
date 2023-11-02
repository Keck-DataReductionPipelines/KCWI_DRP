from keckdrpframework.primitives.base_primitive import BasePrimitive


class ProcessFlat(BasePrimitive):
    """
    Preliminary processing of flat images.
    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        return self.action.args
