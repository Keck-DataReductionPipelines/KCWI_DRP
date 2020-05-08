from keckdrpframework.primitives.base_primitive import BasePrimitive


class ProcessArc(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        if self.action.args.calibration_lamp == \
                self.config.instrument.default_arc_lamp:
            self.logger.info("Using default arc lamp: %s" %
                             self.config.instrument.default_arc_lamp)
            return True
        else:
            self.logger.info("Skipping arc lamp file %s (not the default lamp)"
                             % self.action.args.name)
            self.logger.info("** Resetting next event to None to stop "
                             "processing the file ** ")
            self.action.new_event = None
            return False

    def _perform(self):
        return self.action.args
