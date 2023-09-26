from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer, \
    strip_fname

import numpy as np

from astropy.nddata import CCDData


class NandshuffSubtractSky(BasePrimitive):
    """
    Locate object and sky panels and perform a nod-and-shuffle sky subtraction.

    With a nod-and-shuffle observation, the resulting image will have a sky and
    an object panel with interleaved exposures of the same total duration.  This
    routine identifies those panels and does a subtraction that should eliminate
    the sky accurately, as long as the interleaving occurs at a rate that tracks
    the changes in sky throughout the exposure.

    Writes out a \*_obj.fits and \*_sky.fits image with the panels aligned, but
    not subtracted.  These are subsequently processed through the pipeline to
    allow the confirmation that resulting features are not the result of sky
    subtraction.

    Also writes out the sky-subtracted image in a \*_intk.fits image and adds
    an entry in the proc table.

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if it is a nod-and-shuffle observation
        """
        self.logger.info("Checking precondition for NandshuffSubtractSky")
        if self.action.args.nasmask and self.action.args.numopen > 1:
            self.logger.info("Preconditions for NandshuffSubtractSky met.")
            return True
        else:
            self.logger.warning("Precondition not met: "
                                "not a nod-and-shuffle observation.")
            return False

    def _perform(self):
        self.logger.info("Subtracting nod-and shuffle sky background")

        # Header keyword to update
        key = 'NASSUB'
        keycom = 'Nod-and-shuffle subtraction done?'
        target_type = 'SKY'

        # get header values
        ofn = self.action.args.name
        shrows = self.action.args.ccddata.header['SHUFROWS']
        nshfup = self.action.args.ccddata.header['NSHFUP']
        nshfdn = self.action.args.ccddata.header['NSHFDN']
        camera = self.action.args.ccddata.header['CAMERA'].upper()
        ybin = int(self.action.args.ccddata.header['BINNING'].split(',')[1])

        # units
        u_out = self.action.args.ccddata.unit

        # BLUE camera handling
        if 'BLUE' in camera:
            # BLUE nominal conditions (sky on bottom, object in middle)
            skyrow0 = 0
            skyrow1 = shrows - 1
            objrow0 = shrows
            objrow1 = shrows + shrows - 1
            # record subtraction rows
            subrow0 = objrow0
            subrow1 = objrow1

            # aborted script with inverted panels (sky in middle, object above)
            if nshfdn != nshfup + 1:
                self.logger.warning("Inverted N&S obs detected!")
                skyrow0 = shrows
                skyrow1 = shrows + shrows - 1
                objrow0 = skyrow1 + 1
                objrow1 = objrow0 + shrows - 1
                # record subtraction rows
                subrow0 = skyrow0
                subrow1 = skyrow1

        # RED camera handling
        elif 'RED' in camera:
            # RED uses unbinned pixels for shrows
            if ybin == 2:
                shrows = int(shrows / 2) + 1
            # RED nominal conditions (object on bottom, sky in middle)
            objrow0 = 0
            objrow1 = shrows - 1
            skyrow0 = objrow1 + 1
            skyrow1 = skyrow0 + shrows - 1
            # record subtraction rows
            subrow0 = skyrow0
            subrow1 = skyrow1

            # aborted script with inverted panels (object in middle sky above)
            if nshfdn != nshfup:
                self.logger.warning("Inverted N&S obs detected!")
                objrow0 = shrows
                objrow1 = objrow0 + shrows - 1
                skyrow0 = objrow1 + 1
                skyrow1 = skyrow0 + shrows - 1
                # record subtraction rows
                subrow0 = objrow0
                subrow1 = objrow1

        else:
            self.logger.error("CAMERA cannot be determined for N&S Sub")
            return self.action.args

        # check limits
        self.logger.info(
            "Nod-and-shuffle rows (sky0, 1, obj0,1): %d, %d, %d, %d" %
            (skyrow0, skyrow1, objrow0, objrow1))
        if (skyrow1-skyrow0) != (objrow1-objrow0):
            self.logger.error("Nod-and-shuffle panel mis-match error")
            return self.action.args
        else:
            self.logger.info("Nod-and-shuffle panels match")

        # record in header (and convert to 1-bias, for IRAF)
        self.action.args.ccddata.header['OBJROW0'] = (objrow0+1, "Object row 0")
        self.action.args.ccddata.header['OBJROW1'] = (objrow1+1, "Object row 1")
        self.action.args.ccddata.header['SKYROW0'] = (skyrow0+1, "Sky row 0")
        self.action.args.ccddata.header['SKYROW1'] = (skyrow1+1, "Sky row 1")
        self.action.args.ccddata.header['SUBROW0'] = (subrow0 + 1,
                                                      "Subtraction row 1")
        self.action.args.ccddata.header['SUBROW1'] = (subrow1 + 1,
                                                      "Subtraction row 1")

        # create intermediate images and headers
        sky = self.action.args.ccddata.data.copy()
        obj = self.action.args.ccddata.data.copy()
        std = self.action.args.ccddata.uncertainty.array.copy()
        msk = self.action.args.ccddata.mask.copy()
        flg = self.action.args.ccddata.flags.copy()
        skyhdr = self.action.args.ccddata.header.copy()
        objhdr = self.action.args.ccddata.header.copy()

        # BLUE camera handling
        if 'BLUE' in camera:
            # nominal condition
            if skyrow0 < 10:
                self.logger.info("standard nod-and-shuffle configuration")
                skystd = self.action.args.ccddata.uncertainty.array.copy()
                skymsk = self.action.args.ccddata.mask.copy()
                skyflg = self.action.args.ccddata.flags.copy()
                # move sky to object position
                sky[objrow0:objrow1+1, :] = obj[skyrow0:skyrow1+1, :]
                skystd[objrow0:objrow1+1, :] = std[skyrow0:skyrow1+1, :]
                skymsk[objrow0:objrow1+1, :] = msk[skyrow0:skyrow1+1, :]
                skyflg[objrow0:objrow1+1, :] = flg[skyrow0:skyrow1+1, :]
                # do subtraction
                self.action.args.ccddata.data -= sky
                self.action.args.ccddata.uncertainty.array = np.sqrt(
                    std ** 2 + skystd ** 2)
                self.action.args.ccddata.mask += skymsk
                self.action.args.ccddata.flags |= skyflg
                # clean images
                self.action.args.ccddata.data[skyrow0:skyrow1+1, :] = 0.
                self.action.args.ccddata.data[objrow1+1:-1, :] = 0.
                self.action.args.ccddata.uncertainty.array[
                    skyrow0:skyrow1+1, :] = 0.
                self.action.args.ccddata.uncertainty.array[objrow1+1:-1, :] = 0.
                self.action.args.ccddata.mask[skyrow0:skyrow1+1, :] = 1
                self.action.args.ccddata.mask[objrow1+1:-1, :] = 1
                self.action.args.ccddata.flags[skyrow0:skyrow1+1, :] = 64
                self.action.args.ccddata.flags[objrow1+1:-1, :] = 64
                sky[skyrow0:skyrow1+1, :] = 0.
                sky[objrow1+1:-1, :] = 0.
                obj[skyrow0:skyrow1+1, :] = 0.
                obj[objrow1+1:-1, :] = 0.
            else:
                self.logger.warning(
                    "non-standard nod-and-shuffle configuration")
                skyscl = -1.
                while skyscl < 0.:
                    if self.config.instrument.plot_level >= 2:
                        q = input("Enter scale factor for sky to match obj "
                                  "(float): ")
                        try:
                            skyscl = float(q)
                        except ValueError:
                            self.logger.warning(
                                "Invalid input: %s, try again" % q)
                            skyscl = -1.0
                    else:
                        skyscl = 1.0
                self.logger.info("Sky scaling used = %.2f" % skyscl)
                objstd = self.action.args.ccddata.uncertainty.array.copy()
                objmsk = self.action.args.ccddata.mask.copy()
                objflg = self.action.args.ccddata.flags.copy()
                # move object to sky position
                obj[skyrow0:skyrow1+1, :] = obj[objrow0:objrow1+1, :]
                objstd[skyrow0:skyrow1+1, :] = std[objrow0:objrow1+1, :]
                objmsk[skyrow0:skyrow1+1, :] = msk[objrow0:objrow1+1, :]
                objflg[skyrow0:skyrow1+1, :] = flg[objrow0:objrow1+1, :]
                # do subtraction
                sky *= skyscl
                self.action.args.ccddata.data = obj - sky
                self.action.args.ccddata.uncertainty.array = np.sqrt(
                    std ** 2 + objstd ** 2)
                self.action.args.ccddata.mask += objmsk
                self.action.args.ccddata.flags |= objflg
                # clean images
                self.action.args.ccddata.data[objrow0:-1, :] = 0.
                self.action.args.ccddata.data[0:skyrow0+1, :] = 0.
                self.action.args.ccddata.uncertainty.array[objrow0:-1, :] = 0.
                self.action.args.ccddata.uncertainty.array[0:skyrow0+1, :] = 0.
                self.action.args.ccddata.mask[objrow0:-1, :] = 1
                self.action.args.ccddata.mask[0:skyrow0+1, :] = 1
                self.action.args.ccddata.flags[objrow0:-1, :] = 64
                self.action.args.ccddata.flags[0:skyrow0+1, :] = 64
                sky[objrow0:-1, :] = 0.
                sky[0:skyrow0+1, :] = 0.
                obj[objrow0:-1, :] = 0.
                obj[0:skyrow0+1, :] = 0.
                cmnt = 'Aborted nod-and-shuffle observations'
                objhdr['COMMENT'] = cmnt
                skyhdr['COMMENT'] = cmnt
                self.action.args.ccddata.header['COMMENT'] = cmnt
                skyhdr['NASSCL'] = (skyscl, 'Scale factor applied to sky panel')
                self.action.args.ccddata.header['NASSCL'] = (
                    skyscl, 'Scale factor applied to sky panel')
        elif 'RED' in camera:
            # nominal condition
            if objrow0 < 10:
                self.logger.info("standard nod-and-shuffle configuration")
                objstd = self.action.args.ccddata.uncertainty.array.copy()
                objmsk = self.action.args.ccddata.mask.copy()
                objflg = self.action.args.ccddata.flags.copy()
                # move object to sky position
                obj[skyrow0:skyrow1+1, :] = obj[objrow0:objrow1+1, :]
                objstd[objrow0:objrow1+1, :] = std[skyrow0:skyrow1+1, :]
                objmsk[objrow0:objrow1+1, :] = msk[skyrow0:skyrow1+1, :]
                objflg[objrow0:objrow1+1, :] = flg[skyrow0:skyrow1+1, :]
                # do subtraction
                self.action.args.ccddata.data = obj - sky
                self.action.args.ccddata.uncertainty.array = np.sqrt(
                    std ** 2 + objstd ** 2)
                self.action.args.ccddata.mask += objmsk
                self.action.args.ccddata.flags |= objflg
                # clean images
                self.action.args.ccddata.data[objrow0:objrow1+1, :] = 0.
                self.action.args.ccddata.data[skyrow1+1:-1, :] = 0.
                self.action.args.ccddata.uncertainty.array[
                    objrow0:objrow1+1, :] = 0.
                self.action.args.ccddata.uncertainty.array[skyrow1+1:-1, :] = 0.
                self.action.args.ccddata.mask[objrow0:objrow1+1, :] = 1
                self.action.args.ccddata.mask[skyrow1+1:-1, :] = 1
                self.action.args.ccddata.flags[objrow0:objrow1+1, :] = 64
                self.action.args.ccddata.flags[skyrow1+1:-1, :] = 64
                obj[objrow0:objrow1+1, :] = 0.
                obj[skyrow1+1:-1, :] = 0.
                sky[objrow0:objrow1+1, :] = 0.
                sky[skyrow1+1:-1, :] = 0.
            else:
                self.logger.warning(
                    "non-standard nod-and-shuffle configuration")
                skyscl = -1.
                while skyscl < 0.:
                    if self.config.instrument.plot_level >= 2:
                        q = input("Enter scale factor for sky to match obj "
                                  "(float): ")
                        try:
                            skyscl = float(q)
                        except ValueError:
                            self.logger.warning(
                                "Invalid input: %s, try again" % q)
                            skyscl = -1.0
                    else:
                        skyscl = 1.0
                self.logger.info("Sky scaling used = %.2f" % skyscl)
                skystd = self.action.args.ccddata.uncertainty.array.copy()
                skymsk = self.action.args.ccddata.mask.copy()
                skyflg = self.action.args.ccddata.flags.copy()
                # move sky to object position
                sky[objrow0:objrow1+1, :] = obj[skyrow0:skyrow1+1, :]
                skystd[objrow0:objrow1+1, :] = std[skyrow0:skyrow1+1, :]
                skymsk[objrow0:objrow1+1, :] = msk[skyrow0:skyrow1+1, :]
                skyflg[objrow0:objrow1+1, :] = flg[skyrow0:skyrow1+1, :]
                # do subtraction
                self.action.args.ccddata.data -= sky
                self.action.args.ccddata.uncertainty.array = np.sqrt(
                    std ** 2 + skystd ** 2)
                self.action.args.ccddata.mask += skymsk
                self.action.args.ccddata.flags |= skyflg
                # clean images
                self.action.args.ccddata.data[0:objrow0, :] = 0.
                self.action.args.ccddata.data[objrow1+1:-1, :] = 0.
                self.action.args.ccddata.uncertainty.array[0:objrow0, :] = 0.
                self.action.args.ccddata.uncertainty.array[objrow1+1:-1:] = 0.
                self.action.args.ccddata.mask[0:objrow0, :] = 1
                self.action.args.ccddata.mask[objrow1+1:-1, :] = 1
                self.action.args.ccddata.flags[0:objrow0, :] = 64
                self.action.args.ccddata.flags[objrow1+1:-1, :] = 64
                sky[0, objrow0, :] = 0.
                sky[objrow1+1:-1, :] = 0.
                obj[0, objrow0, :] = 0.
                obj[objrow1+1:-1, :] = 0.
                cmnt = 'Aborted nod-and-shuffle observations'
                objhdr['COMMENT'] = cmnt
                skyhdr['COMMENT'] = cmnt
                self.action.args.ccddata.header['COMMENT'] = cmnt
                skyhdr['NASSCL'] = (skyscl, 'Scale factor applied to sky panel')
                self.action.args.ccddata.header['NASSCL'] = (
                    skyscl, 'Scale factor applied to sky panel')
        else:
            self.logger.error("CAMERA cannot be determined for N&S Sub")
            return self.action.arg

        # log
        self.logger.info("Nod-and-shuffle subtracted")
        # update headers
        log_string = NandshuffSubtractSky.__module__

        objhdr[key] = (False, keycom)
        objhdr['HISTORY'] = log_string

        skyhdr[key] = (False, keycom)
        skyhdr['SKYOBS'] = (True, 'Sky observation?')
        skyhdr['HISTORY'] = log_string

        self.action.args.ccddata.header[key] = (True, keycom)
        self.action.args.ccddata.header['HISTORY'] = log_string

        # write out sky image
        msname = strip_fname(ofn) + '_' + target_type.lower() + '.fits'
        out_sky = CCDData(sky, meta=skyhdr, unit=u_out)
        kcwi_fits_writer(out_sky, output_file=msname,
                         output_dir=self.config.instrument.output_directory)
        # write out object image
        obname = strip_fname(ofn) + '_obj.fits'
        out_obj = CCDData(obj, meta=objhdr, unit=u_out)
        kcwi_fits_writer(out_obj, output_file=obname,
                         output_dir=self.config.instrument.output_directory)

        # update header keywords
        self.action.args.ccddata.header['SKYMAST'] = (msname,
                                                      "Master sky filename")

        # write out int image
        kcwi_fits_writer(self.action.args.ccddata,
                         table=self.action.args.table,
                         output_file=self.action.args.name,
                         output_dir=self.config.instrument.output_directory,
                         suffix="intk")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="intk",
                                            filename=self.action.args.name)
        self.context.proctab.write_proctab(tfil=self.config.instrument.procfile)

        self.logger.info(log_string)

        return self.action.args
    # END: class NandshuffSubtractSky()
