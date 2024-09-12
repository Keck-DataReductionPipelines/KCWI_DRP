from keckdrpframework.primitives.base_primitive import BasePrimitive



class StopPipeline(BasePrimitive):
    """
    Stop the pipeline, if requested by the user
    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
       return True

    def _perform(self):

        self.logger.info("Checking for pipeline stop request")
        self.logger.info(f"Stop pipeline: {self.action.args.stop_pipeline}")


        if self.action.args.stop_pipeline:
            self.logger.info("User requested pipeline stop")    
            # This overrides the event table, so the pipeline will stop processing
            # this file, but still continue with the next file
            self.action.new_event = None        
        return self.action.args