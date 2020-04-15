"""
KCWI

@author: lrizzi
"""

from keckdrpframework.pipelines.base_pipeline import BasePipeline
from keckdrpframework.models.processing_context import ProcessingContext
#from kcwidrp.primitives.kcwi_primitives import *
from kcwidrp.primitives.kcwi_file_primitives import *
#from kcwidrp.primitives.StartBokeh import *
from kcwidrp.core.kcwi_proctab import Proctab


class Kcwi_pipeline(BasePipeline):
    """
    Pipeline to process KCWI data

    """
    name = 'KCWI-DRP'

    event_table = {
        "start_bokeh":               ("StartBokeh", None, None),
        # For every file do this
        "next_file":                 ("ingest_file",
                                      "ingest_file_started",
                                      "file_ingested"),
        "file_ingested":             ("action_planner", None, None),
        # BIAS PROCESSING
        "process_bias":              ("ProcessBias",
                                      "bias_processing_started",
                                      "bias_make_master"),
        "bias_make_master":          ("MakeMasterBias",
                                      "master_bias_started",
                                      None),
        # DARK PROCESSING
        "process_dark":              ("ProcessDark",
                                      "dark_processing_started",
                                      "dark_subtract_bias"),
        "dark_subtract_bias":        ("SubtractBias",
                                      "subtract_bias_started",
                                      "dark_subtract_overscan"),
        "dark_subtract_overscan":    ("SubtractOverscan",
                                      "subtract_overscan_started",
                                      "dark_trim_overscan"),
        "dark_trim_overscan":        ("TrimOverscan",
                                      "trim_overscan_started",
                                      "dark_correct_gain"),
        "dark_correct_gain":         ("CorrectGain",
                                      "gain_correction_started",
                                      "dark_rectify_image"),
        "dark_rectify_image":        ("RectifyImage",
                                      "rectification_started",
                                      "dark_make_master"),
        "dark_make_master":          ("MakeMasterDark",
                                      "master_dark_started",
                                      None),
        # CONTBARS PROCESSING
        "process_contbars":          ("ProcessContbars",
                                      "contbars_processing_started",
                                      "contbar_subtract_overscan"),
        "contbar_subtract_overscan": ("SubtractOverscan",
                                      "subtract_overscan_started",
                                      "contbar_trim_overscan"),
        "contbar_trim_overscan":     ("TrimOverscan",
                                      "trim_overscan_started",
                                      "contbar_correct_gain"),
        "contbar_correct_gain":      ("CorrectGain",
                                      "gain_correction_started",
                                      "contbar_rectify_image"),
        "contbar_rectify_image":     ("RectifyImage",
                                      "rectification_started",
                                      "contbar_find_bars"),
        "contbar_find_bars":         ("FindBars",
                                      "find_bars_started",
                                      "contbar_trace_bars"),
        "contbar_trace_bars":        ("TraceBars",
                                      "trace_bars_started",
                                      None),
        # ARCS PROCESSING
        "process_arc":               ("ProcessArc",
                                      "arcs_processing_started",
                                      "arcs_subtract_overscan"),
        "arcs_subtract_overscan":    ("SubtractOverscan",
                                      "subtract_overscan_started",
                                      "arcs_trim_overscan"),
        "arcs_trim_overscan":        ("TrimOverscan",
                                      "trim_overscan_started",
                                      "arcs_correct_gain"),
        "arcs_correct_gain":         ("CorrectGain",
                                      "gain_correction_started",
                                      "arcs_rectify_image"),
        "arcs_rectify_image":        ("RectifyImage",
                                      "rectification_started",
                                      "arcs_extract_arcs"),
        "arcs_extract_arcs":         ("ExtractArcs",
                                      "extract_arcs_started",
                                      "arcs_arc_offsets"),
        "arcs_arc_offsets":          ("ArcOffsets",
                                      "arc_offset_started",
                                      "arcs_calc_prelim_disp"),
        "arcs_calc_prelim_disp":     ("CalcPrelimDisp",
                                      "prelim_disp_started",
                                      "arcs_read_atlas"),
        "arcs_read_atlas":           ("ReadAtlas",
                                      "read_atlas_started",
                                      "arcs_fit_center"),
        "arcs_fit_center":           ("FitCenter",
                                      "fit_center_started",
                                      "arcs_get_atlas_lines"),
        "arcs_get_atlas_lines":      ("GetAtlasLines",
                                      "getting_atlas_lines_started",
                                      "atlas_solve_arcs"),
        "atlas_solve_arcs":          ("SolveArcs",
                                      "solving_arcs_started",
                                      "arcs_solve_geom"),
        "arcs_solve_geom":           ("SolveGeom",
                                      "solving_geom_started",
                                      "arcs_generate_maps"),
        "arcs_generate_maps":        ("GenerateMaps",
                                      "generating_maps_started",
                                      None),
        # FLAT PROCESSING
        "process_flat":              ("ProcessFlat",
                                      "flat_processing_started",
                                      "flat_subtract_bias"),
        "flat_subtract_bias":        ("SubtractBias",
                                      "subtract_bias started",
                                      "flat_subtract_overscan"),
        "flat_subtract_overscan":    ("SubtractOverscan",
                                      "subtract_overscan_started",
                                      "flat_trim_overscan"),
        "flat_trim_overscan":        ("TrimOverscan",
                                      "trim_overscan_started",
                                      "flat_correct_defects"),
        "flat_correct_defects":      ("CorrectDefects",
                                      "defect_correction_started",
                                      "flat_remove_crs"),
        "flat_remove_crs":           ("RemoveCosmicRays",
                                      "remove_crs_started",
                                      "flat_rectify_image"),
        "flat_rectify_image":        ("RectifyImage",
                                      "rectification_started",
                                      "flat_subtract_dark"),
        "flat_subtract_dark":        ("SubtractDark",
                                      "subtract_dark_started",
                                      "flat_subtract_scat"),
        "flat_subtract_scat":        ("SubtractScatteredLight",
                                      "scat_subtract_started",
                                      "flat_make_stack"),
        "flat_make_stack":           ("StackFlats",
                                      "stack_flats_started",
                                      "flat_make_master"),
        "flat_make_master":          ("MakeMasterFlat",
                                      "master_flat_started",
                                      None),
        # OBJECT PROCESSING
        "process_object":            ("ProcessObject",
                                      "object_processing_started",
                                      "object_subtract_bias"),
        "object_subtract_bias":      ("SubtractBias",
                                      "subtract_bias started",
                                      "object_subtract_overscan"),
        "object_subtract_overscan":  ("SubtractOverscan",
                                      "subtract_overscan_started",
                                      "object_trim_overscan"),
        "object_trim_overscan":      ("TrimOverscan",
                                      "trim_overscan_started",
                                      "object_correct_gain"),
        "object_correct_gain":       ("CorrectGain",
                                      "gain_correction_started",
                                      "object_correct_defects"),
        "object_correct_defects":    ("CorrectDefects",
                                      "defect_correction_started",
                                      "object_remove_crs"),
        "object_remove_crs":         ("RemoveCosmicRays",
                                      "remove_crs_started",
                                      "object_create_unc"),
        "object_create_unc":         ("CreateUncertaintyImage",
                                      "create_unc_started",
                                      "object_rectify_image"),
        "object_rectify_image":      ("RectifyImage",
                                      "rectification_started",
                                      "object_subtract_dark"),
        "object_subtract_dark":      ("SubtractDark",
                                      "subtract_dark started",
                                      "object_subtract_scat"),
        "object_subtract_scat":      ("SubtractScatteredLight",
                                      "scat_subtract_started",
                                      "object_make_cube"),
        "object_make_cube":          ("MakeCube",
                                      "making_cube_started",
                                      "object_correct_dar"),
        "object_correct_dar":        ("CorrectDar",
                                      "correcting_dar_started",
                                      "object_flux_calibrate"),
        "object_flux_calibrate":     ("FluxCalibrate",
                                      "flux_calibration_started",
                                      None),
        "next_file_stop":            ("ingest_file", "file_ingested", None)
    }

    # event_table = kcwi_event_table

    def __init__(self, context: ProcessingContext):
        """
        Constructor
        """
        BasePipeline.__init__(self, context)
        self.cnt = 0

    def action_planner(self, action, context):
        try:
            self.logger.info("******* FILE TYPE DETERMINED AS %s" %
                             action.args.imtype)
        except:
            return

        groupid = action.args.groupid
        self.logger.info("******* GROUPID is %s " % action.args.groupid)
        if action.args.in_proctab:
            self.logger.warn("Already processed (already in proctab")
        if action.args.in_proctab and not context.config.instrument.clobber:
            self.logger.warn("Pushing noop to queue")
            context.push_event("noop", action.args)
        elif action.args.imtype == "BIAS":
            bias_args = action.args
            bias_args.groupid = groupid
            bias_args.want_type = "BIAS"
            bias_args.new_type = "MBIAS"
            bias_args.min_files = context.config.instrument.bias_min_nframes
            bias_args.new_file_name = "master_bias_%s.fits" % groupid
            context.push_event("process_bias", bias_args)
        elif action.args.imtype == "DARK":
            dark_args = action.args
            dark_args.groupid = groupid
            dark_args.want_type = "DARK"
            dark_args.new_type = "MDARK"
            dark_args.min_files = context.config.instrument.dark_min_nframes
            dark_args.new_file_name = "master_dark_%s.fits" % groupid
            dark_args.in_directory = "redux"
            context.push_event("process_dark", dark_args)
        elif "CONTBARS" in action.args.imtype:
            context.push_event("process_contbars", action.args)
        elif "FLATLAMP" in action.args.imtype:
            flat_args = action.args
            flat_args.groupid = groupid
            flat_args.want_type = "FLATLAMP"
            flat_args.new_type = "MFLAT"
            flat_args.min_files = context.config.instrument.flat_min_nframes
            flat_args.new_file_name = "master_flat_%s.fits" % groupid
            flat_args.in_directory = "redux"
            context.push_event("process_flat", flat_args)
        elif "ARCLAMP" in action.args.imtype:
            context.push_event("process_arc", action.args)
        elif "OBJECT" in action.args.imtype:
            context.push_event("process_object", action.args)


if __name__ == "__main__":
    """
    Standalone test
    """
    pass
