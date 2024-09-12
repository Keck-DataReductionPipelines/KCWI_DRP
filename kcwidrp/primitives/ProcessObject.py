from keckdrpframework.primitives.base_primitive import BasePrimitive


class ProcessObject(BasePrimitive):
    """
    Preliminary processing of object images.
    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        import pdb
        print("ProcessObject new_event:")
        print(self.action.new_event)
        pdb.set_trace()
        return self.action.args
