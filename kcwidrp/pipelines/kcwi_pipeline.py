"""
Class definition for the KCWI Pipeline

The event_table defines the recipes with each type of image having a unique
entry point and trajectory in the table.

The action_planner sets up the parameters for each image type before the event
queue in the KeckDRPFramework gets started.

@author: lrizzi, jneill
"""

from keckdrpframework.pipelines.base_pipeline import BasePipeline
from keckdrpframework.models.processing_context import ProcessingContext
from kcwidrp.primitives.kcwi_file_primitives import *


class Kcwi_pipeline(BasePipeline):
    """
    Pipeline to process KCWI data

    """
    name = 'KCWI-DRP'

    event_table = {
        # this method is used with the "group" option,
        # to ingest the data without triggering any processing.
        # it is defined lower in this file
        "add_only":                  ("add_to_dataframe_only", None, None),
        #
        "start_bokeh":               ("StartBokeh", None, None),
        # For every file do this
        "next_file":                 ("ingest_file",
                                      "ingest_file_started",
                                      "file_ingested"),
        "file_ingested":             ("action_planner", None, None),
        # BIAS PROCESSING
        "process_bias":              ("ProcessBias",
                                      "bias_processing_started",
                                      "bias_subtract_overscan"),
        "bias_subtract_overscan":    ("SubtractOverscan",
                                      "subtract_overscan_started",
                                      "bias_trim_overscan"),
        "bias_trim_overscan":        ("TrimOverscan",
                                      "trim_overscan_started",      # intb
                                      "bias_make_master"),
        "bias_make_master":          ("MakeMasterBias",
                                      "master_bias_started",        # mbias
                                      None),
        # DARK PROCESSING
        "process_dark":              ("ProcessDark",
                                      "dark_processing_started",
                                      "dark_subtract_overscan"),
        "dark_subtract_overscan":    ("SubtractOverscan",
                                      "subtract_overscan_started",
                                      "dark_trim_overscan"),
        "dark_trim_overscan":        ("TrimOverscan",
                                      "trim_overscan_started",
                                      "dark_subtract_bias"),
        "dark_subtract_bias":        ("SubtractBias",
                                      "subtract_bias_started",
                                      "dark_correct_gain"),
        "dark_correct_gain":         ("CorrectGain",
                                      "gain_correction_started",
                                      "dark_correct_defects"),
        "dark_correct_defects":      ("CorrectDefects",
                                      "defect_correction_started",
                                      "dark_remove_crs"),
        "dark_remove_crs":           ("RemoveCosmicRays",
                                      "remove_crs_started",
                                      "dark_create_unc"),
        "dark_create_unc":           ("CreateUncertaintyImage",
                                      "create_unc_started",
                                      "dark_rectify_image"),
        "dark_rectify_image":        ("RectifyImage",
                                      "rectification_started",      # int
                                      "dark_make_master"),
        "dark_make_master":          ("MakeMasterDark",
                                      "master_dark_started",        # mdark
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
                                      "rectification_started",      # int
                                      "contbar_make_master"),
        "contbar_make_master":       ("MakeMasterContbars",
                                      "master_contbar_started",     # mcont
                                      "contbar_find_bars"),
        "contbar_find_bars":         ("FindBars",
                                      "find_bars_started",
                                      "contbar_trace_bars"),
        "contbar_trace_bars":        ("TraceBars",
                                      "trace_bars_started",         # trace
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
                                      "arcs_correct_defects"),
        "arcs_correct_defects":      ("CorrectDefects",
                                      "defect_correction_started",
                                      "arcs_create_unc"),
        "arcs_create_unc":           ("CreateUncertaintyImage",
                                      "create_unc_started",
                                      "arcs_rectify_image"),
        "arcs_rectify_image":        ("RectifyImage",
                                      "rectification_started",      # int
                                      "arcs_make_master"),
        "arcs_make_master":          ("MakeMasterArc",
                                      "master_arcs_started",        # marc
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
                                      "solving_geom_started",       # geom
                                      "arcs_generate_maps"),
        "arcs_generate_maps":        ("GenerateMaps",
                                      "generating_maps_started",    # maps
                                      "arc_make_cube"),
        "arc_make_cube":             ("MakeCube",
                                      "making_cube_started",        # icube
                                      "arc_make_cubeimage"),
        "arc_make_cubeimage":        ("CubeImage",
                                      "making_cubeimage_started",   # icube_2d
                                      None),
        # FLAT PROCESSING
        "process_flat":              ("ProcessFlat",
                                      "flat_processing_started",
                                      "flat_subtract_overscan"),
        "flat_subtract_overscan":    ("SubtractOverscan",
                                      "subtract_overscan_started",
                                      "flat_trim_overscan"),
        "flat_trim_overscan":        ("TrimOverscan",
                                      "trim_overscan_started",
                                      "flat_subtract_bias"),
        "flat_subtract_bias":        ("SubtractBias",
                                      "subtract_bias started",
                                      "flat_correct_gain"),
        "flat_correct_gain":         ("CorrectGain",
                                      "gain_correction_started",
                                      "flat_correct_defects"),
        "flat_correct_defects":      ("CorrectDefects",
                                      "defect_correction_started",
                                      "flat_remove_crs"),
        "flat_remove_crs":           ("RemoveCosmicRays",
                                      "remove_crs_started",
                                      "flat_create_unc"),
        "flat_create_unc":           ("CreateUncertaintyImage",
                                      "create_unc_started",
                                      "flat_rectify_image"),
        "flat_rectify_image":        ("RectifyImage",
                                      "rectification_started",      # int
                                      "flat_subtract_dark"),
        "flat_subtract_dark":        ("SubtractDark",
                                      "subtract_dark_started",
                                      "flat_subtract_scat"),
        "flat_subtract_scat":        ("SubtractScatteredLight",
                                      "scat_subtract_started",      # intd
                                      "flat_make_stack"),
        "flat_make_stack":           ("StackFlats",
                                      "stack_flats_started",        # sflat
                                      "flat_make_master"),
        "flat_make_master":          ("MakeMasterFlat",
                                      "master_flat_started",        # mflat
                                      "flat_correct_illumination"),
        "flat_correct_illumination": ("CorrectIllumination",
                                      "illumination_correction_started",  # intf
                                      "flat_make_cube"),
        "flat_make_cube":            ("MakeCube",
                                      "making_flat_cube_started",   # icube
                                      None),
        # OBJECT PROCESSING
        "process_object":            ("ProcessObject",
                                      "object_processing_started",
                                      "object_flag_saturation"),
        "object_flag_saturation":    ("FlagSaturation",
                                      "object_flag_saturation_started",
                                      "object_subtract_overscan"),
        "object_subtract_overscan":  ("SubtractOverscan",
                                      "subtract_overscan_started",
                                      "object_trim_overscan"),
        "object_trim_overscan":      ("TrimOverscan",
                                      "trim_overscan_started",
                                      "object_subtract_bias"),
        "object_subtract_bias":      ("SubtractBias",
                                      "subtract_bias started",
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
                                      "rectification_started",      # int
                                      "object_subtract_dark"),
        "object_subtract_dark":      ("SubtractDark",
                                      "subtract_dark started",
                                      "object_subtract_scat"),
        "object_subtract_scat":      ("SubtractScatteredLight",
                                      "scat_subtract_started",      # intd
                                      "object_correct_illum"),
        "object_correct_illum":      ("CorrectIllumination",
                                      "illumination_correction_started",  # intf
                                      "object_make_master"),
        "object_make_master":        ("MakeMasterObject",
                                      "master_object_started",      # mobj
                                      "object_make_sky"),
        "object_make_sky":           ("MakeMasterSky",
                                      "making_master_sky_started",
                                      "object_subtract_sky"),
        "object_subtract_sky":       ("SubtractSky",
                                      "subtracting_sky_started",    # intk
                                      "object_make_cube"),
        "object_make_cube":          ("MakeCube",
                                      "making_cube_started",        # icube
                                      "object_wavelengthcorr"),
        "object_wavelengthcorr":     ("WavelengthCorrections",
                                      "wavelength_correction_started",  # icubew
                                      "object_correct_dar"),
        "object_correct_dar":        ("CorrectDar",
                                      "correcting_dar_started",     # icubed
                                      "object_make_invsens"),
        "object_make_invsens":       ("MakeInvsens",
                                      "make_invsens_started",
                                      "object_flux_calibrate"),
        "object_flux_calibrate":     ("FluxCalibrate",
                                      "flux_calibration_started",   # icubes
                                      None),
        # NOD AND SHUFFLE OBJECT PROCESSING
        "process_nandshuff":         ("ProcessObject",
                                      "nandshuff_processing_started",
                                      "nandshuff_subtract_overscan"),
        "nandshuff_subtract_overscan": ("SubtractOverscan",
                                        "subtract_overscan_started",
                                        "nandshuff_trim_overscan"),
        "nandshuff_trim_overscan":   ("TrimOverscan",
                                      "trim_overscan_started",
                                      "nandshuff_subtract_bias"),
        "nandshuff_subtract_bias":   ("SubtractBias",
                                      "subtract_bias started",
                                      "nandshuff_correct_gain"),
        "nandshuff_correct_gain":    ("CorrectGain",
                                      "gain_correction_started",
                                      "nandshuff_correct_defects"),
        "nandshuff_correct_defects": ("CorrectDefects",
                                      "defect_correction_started",
                                      "nandshuff_remove_crs"),
        "nandshuff_remove_crs":      ("RemoveCosmicRays",
                                      "remove_crs_started",
                                      "nandshuff_create_unc"),
        "nandshuff_create_unc":      ("CreateUncertaintyImage",
                                      "create_unc_started",
                                      "nandshuff_rectify_image"),
        "nandshuff_rectify_image":   ("RectifyImage",
                                      "rectification_started",      # int
                                      "nandshuff_subtract_sky"),
        "nandshuff_subtract_sky":    ("NandshuffSubtractSky",
                                      "nandshuff_skysub_started",   # intk
                                      "nandshuff_correct_illum"),
        "nandshuff_correct_illum":   ("CorrectIllumination",
                                      "illumination_correction_started",  # intf
                                      "nandshuff_make_cube"),
        "nandshuff_make_cube":       ("MakeCube",
                                      "making_cube_started",        # icube
                                      "nandshuff_wavelengthcorr"),
        "nandshuff_wavelengthcorr":  ("WavelengthCorrections",
                                      "wavelength_correction_started",  # icubew
                                      "nandshuff_correct_dar"),
        "nandshuff_correct_dar":     ("CorrectDar",
                                      "correcting_dar_started",     # icubed
                                      "nandshuff_flux_calibrate"),
        "nandshuff_flux_calibrate":  ("FluxCalibrate",
                                      "flux_calibration_started",   # icubes
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

    def add_to_dataframe_only(self, action, context):
        self.context.pipeline_logger.info("******* ADD to DATAFRAME ONLY: %s" %
                                          action.args.name)
        return action.args

    def action_planner(self, action, context):
        try:
            self.context.pipeline_logger.info(
                "******* FILE TYPE DETERMINED AS %s" % action.args.imtype)
        except (AttributeError, TypeError, ValueError):
            self.context.pipeline_logger.warn(
                "******* FILE TYPE is NOT determined. "
                "No processing is possible.")
            return False

        groupid = action.args.groupid
        camera = action.args.ccddata.header['CAMERA'].upper()
        self.context.pipeline_logger.info("******* GROUPID is %s " %
                                          action.args.groupid)
        self.context.pipeline_logger.info(
            "******* STATEID is %s (%s) " %
            (action.args.ccddata.header["STATENAM"],
             action.args.ccddata.header["STATEID"]))
        self.context.pipeline_logger.info("******* CAMERA is %s " % camera)
        if action.args.in_proctab:
            if len(action.args.last_suffix) > 0:
                self.context.pipeline_logger.warn(
                    "Already processed (already in proctab up to %s)" %
                    action.args.last_suffix)
            else:
                self.context.pipeline_logger.warn(
                    "Already processed (already in proctab)")
        if action.args.in_proctab and not context.config.instrument.clobber:
            self.context.pipeline_logger.warn("Pushing noop to queue")
            context.push_event("noop", action.args)
        elif "BIAS" in action.args.imtype:
            if action.args.ttime > 0:
                self.context.pipeline_logger.warn(
                    f"BIAS frame with exposure time = {action.args.ttime} "
                    f"> 0. Discarding.")
                return False
            bias_args = action.args
            bias_args.groupid = groupid
            bias_args.want_type = "BIAS"
            bias_args.new_type = "MBIAS"
            bias_args.min_files = context.config.instrument.bias_min_nframes
            bias_args.new_file_name = "master_bias_%s.fits" % groupid
            context.push_event("process_bias", bias_args)
        elif "DARK" in action.args.imtype:
            dark_args = action.args
            dark_args.groupid = groupid
            dark_args.want_type = "DARK"
            dark_args.new_type = "MDARK"
            dark_args.min_files = context.config.instrument.dark_min_nframes
            dark_args.new_file_name = "master_dark_%s.fits" % groupid
            dark_args.in_directory = "redux"
            context.push_event("process_dark", dark_args)
        elif "CONTBARS" in action.args.imtype:
            contbars_args = action.args
            contbars_args.groupid = groupid
            contbars_args.want_type = "CONTBARS"
            contbars_args.new_type = "MCBARS"
            contbars_args.min_files = \
                context.config.instrument.contbars_min_nframes
            contbars_args.new_file_name = "master_contbars_%s.fits" % groupid
            contbars_args.in_directory = "redux"
            context.push_event("process_contbars", contbars_args)
        elif "FLATLAMP" in action.args.imtype:
            flat_args = action.args
            flat_args.groupid = groupid
            flat_args.want_type = "FLATLAMP"
            flat_args.stack_type = "SFLAT"
            flat_args.new_type = "MFLAT"
            flat_args.min_files = context.config.instrument.flat_min_nframes
            flat_args.new_file_name = "master_flat_%s.fits" % groupid
            flat_args.in_directory = "redux"
            context.push_event("process_flat", flat_args)
        elif "DOMEFLAT" in action.args.imtype:
            flat_args = action.args
            flat_args.groupid = groupid
            flat_args.want_type = "DOMEFLAT"
            flat_args.stack_type = "SDOME"
            flat_args.new_type = "MDOME"
            flat_args.min_files = context.config.instrument.dome_min_nframes
            flat_args.new_file_name = "master_flat_%s.fits" % groupid
            flat_args.in_directory = "redux"
            context.push_event("process_flat", flat_args)
        elif "TWIFLAT" in action.args.imtype:
            flat_args = action.args
            flat_args.groupid = groupid
            flat_args.want_type = "TWIFLAT"
            flat_args.stack_type = "STWIF"
            flat_args.new_type = "MTWIF"
            flat_args.min_files = context.config.instrument.twiflat_min_nframes
            flat_args.new_file_name = "master_flat_%s.fits" % groupid
            flat_args.in_directory = "redux"
            context.push_event("process_flat", flat_args)
        elif "ARCLAMP" in action.args.imtype:
            arc_args = action.args
            arc_args.groupid = groupid
            arc_args.want_type = "ARCLAMP"
            arc_args.new_type = "MARC"
            arc_args.min_files = context.config.instrument.arc_min_nframes
            arc_args.new_file_name = "master_arc_%s.fits" % groupid
            arc_args.in_directory = "redux"
            context.push_event("process_arc", arc_args)
        elif "OBJECT" in action.args.imtype:
            if action.args.nasmask and action.args.numopen > 1:
                context.push_event("process_nandshuff", action.args)
            else:
                object_args = action.args
                # object_args.new_type = "SKY"
                object_args.new_type = "MOBJ"
                object_args.min_files = \
                    context.config.instrument.object_min_nframes
                object_args.in_directory = "redux"
                context.push_event("process_object", object_args)
        return True


if __name__ == "__main__":
    """
    Standalone test
    """
    pass
