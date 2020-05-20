from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer

import os


class CorrectIllumination(BasePrimitive):
    """Subtract master bias frame"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can correct illumination based on the processing table
        :return:
        """
        self.action.args.master_flat = None
        self.logger.info("Checking precondition for CorrectIllumination")
        # first check for internal flat
        target_type = 'MFLAT'
        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type=target_type,
                                             nearest=True)
        if len(tab) <= 0:
            # next look for twilight flat
            target_type = 'MTWIF'
            tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                                 target_type=target_type,
                                                 nearest=True)
            if len(tab) <= 0:
                # finally look for dome flat
                target_type = 'MDOME'
                tab = self.context.proctab.n_proctab(
                    frame=self.action.args.ccddata,
                    target_type=target_type,
                    nearest=True)
                if len(tab) <= 0:
                    precondition = False
                else:
                    precondition = True
            else:
                precondition = True
        else:
            precondition = True

        self.logger.info("pre condition got %d %s flats, expected 1"
                         % (len(tab), target_type))
        if precondition:
            self.action.args.master_flat = tab['OFNAME'][0].split('.')[0] + \
                                           '_' + target_type.lower() + ".fits"
        return precondition

    def _perform(self):

        # Header keyword to update
        key = 'FLATCOR'
        keycom = 'flat corrected?'
        # obj, sky
        obj = None
        sky = None

        self.logger.info("Correcting Illumination")
        if self.action.args.master_flat:
            mflat = kcwi_fits_reader(
                os.path.join(os.path.dirname(self.action.args.name),
                             self.config.instrument.output_directory,
                             self.action.args.master_flat))[0]

            # do the correction
            self.action.args.ccddata.data *= mflat.data

            # update header keywords
            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['MFFILE'] = (
                self.action.args.master_flat, "Master flat filename")

            # check for obj, sky observations
            if self.action.args.nasmask and self.action.args.numopen > 1:
                ofn = self.action.args.ccddata.header['OFNAME']

                objfn = ofn.split('.')[0] + '_obj.fits'
                full_path = os.path.join(
                    os.path.dirname(self.action.args.name),
                    self.config.instrument.output_directory, objfn)
                if os.path.exists(full_path):
                    obj = kcwi_fits_reader(full_path)[0]
                    # correction
                    obj.data *= mflat.data
                    # update header
                    obj.header[key] = (True, keycom)
                    obj.header['MFFILE'] = (
                        self.action.args.master_flat, 'Master flat filename')
                else:
                    obj = None

                skyfn = ofn.split('.')[0] + '_sky.fits'
                full_path = os.path.join(
                    os.path.dirname(self.action.args.name),
                    self.config.instrument.output_directory, skyfn)
                if os.path.exists(full_path):
                    sky = kcwi_fits_reader(full_path)[0]
                    # correction
                    sky.data *= mflat.data
                    # update header
                    sky.header[key] = (True, keycom)
                    sky.header['MFFILE'] = (
                        self.action.args.master_flat, 'Master flat filename')
                else:
                    sky = None
        else:
            self.logger.error("No master flat found, "
                              "cannot correct illumination.")
            self.action.args.ccddata.header[key] = (False, keycom)

        log_string = CorrectIllumination.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string

        # write out intf image
        kcwi_fits_writer(self.action.args.ccddata,
                         table=self.action.args.table,
                         output_file=self.action.args.name,
                         output_dir=self.config.instrument.output_directory,
                         suffix="intf")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="intf")
        self.context.proctab.write_proctab()

        # check for obj, sky images
        if obj is not None:
            kcwi_fits_writer(obj, output_file=self.action.args.name,
                             output_dir=self.config.instrument.output_directory,
                             suffix="objf")
        if sky is not None:
            kcwi_fits_writer(sky, output_file=self.action.args.name,
                             output_dir=self.config.instrument.output_directory,
                             suffix="skyf")

        self.logger.info(log_string)

        return self.action.args
    # END: class CorrectIllumination()
