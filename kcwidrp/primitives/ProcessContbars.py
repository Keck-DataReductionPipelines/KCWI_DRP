from keckdrpframework.primitives.base_primitive import BasePrimitive


class ProcessContbars(BasePrimitive):
    """
    Preliminary processing of continuum bars images.
    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        return self.action.args
