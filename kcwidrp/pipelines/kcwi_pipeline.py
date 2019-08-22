"""
KCWI

@author: lrizzi
"""

from keckdrpframework.pipelines.base_pipeline import Base_pipeline

from ..primitives.kcwi_primitives import *
from ..primitives.kcwi_file_primitives import *



class Kcwi_pipeline(Base_pipeline):
    """
    Pipeline to process KCWI data

    """

    event_table = {
        #"next_file": ("ingest_file", "file_ingested", None),
        "next_file": ("ingest_file", "file_ingested", "file_ingested"),
        "file_ingested": ("action_planner", None, None),
        # BIAS
        "process_bias": ("process_bias", None, None),
        # CONTBARS PROCESSING
        "process_contbars": ("process_contbars", "contbars_processing_started", "contbar_subtract_overscan"),
        "contbar_subtract_overscan": ("subtract_overscan", "subtract_overscan_started", "contbar_trim_overscan"),
        "contbar_trim_overscan": ("trim_overscan", "trim_overscan_started", "contbar_correct_gain"),
        "contbar_correct_gain": ("correct_gain", "gain_correction_started", "contbar_find_bars"),
        "contbar_find_bars": ("find_bars", "find_bars_started", "contbar_trace_bars"),
        "contbar_trace_bars": ("trace_bars", "trace_bars_started", None),
        # ARCS PROCESSING
        "process_arc": ("process_arc", "arcs_processing_started", "arcs_subtract_overscan"),
        "arcs_subtract_overscan": ("subtract_overscan", "subtract_overscan_started", "arcs_trim_overscan"),
        "arcs_trim_overscan": ("trim_overscan", "trim_overscan_started", "arcs_correct_gain"),
        "arcs_correct_gain": ("correct_gain", "gain_correction_started", "arcs_extract_arcs"),
        "arcs_extract_arcs": ("extract_arcs", "extract_arcs_started", "arcs_arc_offsets"),
        "arcs_arc_offsets":  ("arc_offsets", "arc_offset_started", "arcs_calc_prelim_disp"),
        "arcs_calc_prelim_disp": ("calc_prelim_disp", "prelim_disp_started", "arcs_read_atlas"),
        "arcs_read_atlas": ("read_atlas", "read_atlas_started", "arcs_fit_center"),
        "arcs_fit_center": ("fit_center", "fit_center_started", None),
        # FLAT
        "process_flat": ("process_flat", None, None),
        # General processing
        "next_file_stop": ("ingest_file", "file_ingested", None)
        #"process_object": ("process_object", None, "save_png"),
        #"save_png": ("save_png", None, None)
    }

    #event_table = kcwi_event_table



    def __init__(self):
        """
        Constructor
        """
        Base_pipeline.__init__(self)
        self.cnt = 0

    def action_planner (self, action, context):
        self.logger.info("******* FILE TYPE DETERMINED AS %s" % action.args.imtype)
        groupid = action.args.groupid
        self.logger.info("******* GROUPID is %s " % action.args.groupid)
        if action.args.imtype == "BIAS":
            bias_args = Arguments(name="bias_args",
                                  groupid = groupid,
                                  want_type="BIAS",
                                  new_type="MASTER_BIAS",
                                  min_files=context.config.instrument.bias_min_nframes,
                                  new_file_name="master_bias_%s.fits" % groupid)
            context.push_event("process_bias", bias_args)
        elif "CONTBARS" in action.args.imtype:
            context.push_event("process_contbars", action.args)
        elif "FLAT" in action.args.imtype:
            context.push_event("process_flat", action.args)
        elif "ARCLAMP" in action.args.imtype:
            context.push_event("process_arc", action.args)
        elif "OBJECT" in action.args.imtype:
            context.push_event("process_object", action.args)



if __name__ == "__main__":
    """
    Standalone test
    """
    pass