from keckdrpframework.primitives.base_primitive import BasePrimitive


class ProcessObject(BasePrimitive):
    """
    Preliminary processing of object images.
    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        return self.action.args
