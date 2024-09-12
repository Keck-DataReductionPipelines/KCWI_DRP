from keckdrpframework.primitives.base_primitive import BasePrimitive



class StopPipeline(BasePrimitive):
    """
    Stop the pipeline, if requested by the user
    """

    def __init__(self, action, context):
        self.logger = context.pipeline_logger

    def _pre_condition(self):
       return True

    def _perform(self):

        if self.action.args.stop_pipeline:
            self.logger.info("User requested pipeline stop")    
            exit()
        
        return self.action.args