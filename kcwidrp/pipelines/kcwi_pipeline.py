"""
KCWI

@author: lrizzi
"""

from keckdrpframework.pipelines.base_pipeline import BasePipeline

from ..primitives.kcwi_primitives import *
from ..primitives.kcwi_file_primitives import *
from ..primitives.kcwi_bokeh import *
from ..core.kcwi_proctab import Proctab


class Kcwi_pipeline(BasePipeline):
    """
    Pipeline to process KCWI data

    """

    event_table = {
        "start_bokeh": ("start_bokeh", None, None),
        # "next_file": ("ingest_file", "file_ingested", None),
        "next_file": ("ingest_file", "file_ingested", "file_ingested"),
        "file_ingested": ("action_planner", None, None),
        # BIAS PROCESSING
        "process_bias": ("process_bias", None, None),
        # CONTBARS PROCESSING
        "process_contbars": ("process_contbars",
                             "contbars_processing_started",
                             "contbar_subtract_overscan"),
        "contbar_subtract_overscan": ("subtract_overscan",
                                      "subtract_overscan_started",
                                      "contbar_trim_overscan"),
        "contbar_trim_overscan": ("trim_overscan",
                                  "trim_overscan_started",
                                  "contbar_correct_gain"),
        "contbar_correct_gain": ("correct_gain",
                                 "gain_correction_started",
                                 "contbar_rectify_image"),
        "contbar_rectify_image": ("rectify_image",
                                  "rectification_started",
                                  "contbar_find_bars"),
        "contbar_find_bars": ("find_bars",
                              "find_bars_started",
                              "contbar_trace_bars"),
        "contbar_trace_bars": ("trace_bars", "trace_bars_started", None),
        # ARCS PROCESSING
        "process_arc": ("process_arc",
                        "arcs_processing_started",
                        "arcs_subtract_overscan"),
        "arcs_subtract_overscan": ("subtract_overscan",
                                   "subtract_overscan_started",
                                   "arcs_trim_overscan"),
        "arcs_trim_overscan": ("trim_overscan",
                               "trim_overscan_started",
                               "arcs_correct_gain"),
        "arcs_correct_gain": ("correct_gain",
                              "gain_correction_started",
                              "arcs_rectify_image"),
        "arcs_rectify_image": ("rectify_image",
                               "rectification_started",
                               "arcs_extract_arcs"),
        "arcs_extract_arcs": ("extract_arcs",
                              "extract_arcs_started",
                              "arcs_arc_offsets"),
        "arcs_arc_offsets":  ("arc_offsets",
                              "arc_offset_started",
                              "arcs_calc_prelim_disp"),
        "arcs_calc_prelim_disp": ("calc_prelim_disp",
                                  "prelim_disp_started",
                                  "arcs_read_atlas"),
        "arcs_read_atlas": ("read_atlas",
                            "read_atlas_started",
                            "arcs_fit_center"),
        "arcs_fit_center": ("fit_center",
                            "fit_center_started", None),
        # FLAT PROCESSING
        "process_flat": ("process_flat", None, None),
        # OBJECT PROCESSING
        "process_object": ("process_object",
                           "object_processing_started",
                           "object_subtract_bias"),
        "object_subtract_bias": ("subtract_bias",
                                 "subtract_bias started",
                                 "object_subtract_overscan"),
        "object_subtract_overscan": ("subtract_overscan",
                                     "subtract_overscan_started",
                                     "object_trim_overscan"),
        "object_trim_overscan": ("trim_overscan",
                                 "trim_overscan_started",
                                 "object_correct_gain"),
        "object_correct_gain": ("correct_gain",
                                "gain_correction_started",
                                "object_rectify_image"),
        "object_rectify_image": ("rectify_image",
                                 "rectification_started",
                                 None),
        "next_file_stop": ("ingest_file", "file_ingested", None)
        # "save_png": ("save_png", None, None)
    }

    # event_table = kcwi_event_table

    def __init__(self):
        """
        Constructor
        """
        BasePipeline.__init__(self)
        self.cnt = 0

    def action_planner(self, action, context):
        self.logger.info("******* FILE TYPE DETERMINED AS %s" %
                         action.args.imtype)
        groupid = action.args.groupid
        self.logger.info("******* GROUPID is %s " % action.args.groupid)
        if action.args.imtype == "BIAS":
            bias_args = action.args
            bias_args.groupid = groupid
            bias_args.want_type = "BIAS"
            bias_args.new_type = "MASTER_BIAS"
            bias_args.min_files = context.config.instrument.bias_min_nframes
            bias_args.new_file_name = "master_bias_%s.fits" % groupid
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
