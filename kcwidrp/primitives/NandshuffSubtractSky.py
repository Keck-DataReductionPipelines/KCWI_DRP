from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer

import numpy as np
import os

from astropy.nddata import CCDData
from astropy import units as u


class NandshuffSubtractSky(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if it is a nod-and-shuffle observation
        :return:
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
        ofn = self.action.args.ccddata.header['OFNAME']
        shrows = self.action.args.ccddata.header['SHUFROWS']
        nshfup = self.action.args.ccddata.header['NSHFUP']
        nshfdn = self.action.args.ccddata.header['NSHFDN']

        # units
        u_out = self.action.args.ccddata.unit

        # nominal conditions (sky on bottom, object in middle)
        skyrow0 = 0
        skyrow1 = shrows - 1
        objrow0 = shrows
        objrow1 = shrows + shrows - 1

        # aborted script with inverted panels (sky in middle, object above)
        if nshfdn != nshfup + 1:
            skyrow0 = shrows
            skyrow1 = shrows + shrows - 1
            objrow0 = skyrow1
            objrow1 = objrow0 + shrows - 1

        # check limits
        if (skyrow1-skyrow0) != (objrow1-objrow0):
            self.logger.error("Nod-and-shuffle row limits error")
            return self.action.args

        # create intermediate images and headers
        sky = self.action.args.ccddata.data.copy()
        obj = self.action.args.ccddata.data.copy()
        std = self.action.args.ccddata.uncertainty.array.copy()
        msk = self.action.args.ccddata.mask.copy()
        flg = self.action.args.ccddata.flags.copy()
        skyhdr = self.action.args.ccddata.header.copy()
        objhdr = self.action.args.ccddata.header.copy()

        # nominal condition
        if skyrow0 < 10:
            self.logger.info("standard nod-and-shuffle configuration")
            skystd = self.action.args.ccddata.uncertainty.array.copy()
            skymsk = self.action.args.ccddata.mask.copy()
            skyflg = self.action.args.ccddata.flags.copy()
            # move sky to object position
            sky[objrow0:objrow1, :] = obj[skyrow0:skyrow1, :]
            skystd[objrow0:objrow1, :] = std[skyrow0:skyrow1, :]
            skymsk[objrow0:objrow1, :] = msk[skyrow0:skyrow1, :]
            skyflg[objrow0:objrow1, :] = flg[skyrow0:skyrow1, :]
            # do subtraction
            self.action.args.ccddata.data -= sky
            self.action.args.ccddata.uncertainty.array = np.sqrt(std ** 2 +
                                                                 skystd ** 2)
            self.action.args.ccddata.mask += skymsk
            self.action.args.ccddata.flags |= skyflg
            # clean images
            self.action.args.ccddata.data[skyrow0:skyrow1, :] = 0.
            self.action.args.ccddata.data[(objrow1+1):-1, :] = 0.
            self.action.args.ccddata.uncertainty.array[skyrow0:skyrow1, :] = 0.
            self.action.args.ccddata.uncertainty.array[(objrow1+1):-1, :] = 0.
            self.action.args.ccddata.mask[skyrow0:skyrow1, :] = 1
            self.action.args.ccddata.mask[(objrow1 + 1):-1, :] = 1
            self.action.args.ccddata.flags[skyrow0:skyrow1, :] = 64
            self.action.args.ccddata.flags[(objrow1 + 1):-1, :] = 64
            sky[skyrow0:skyrow1, :] = 0.
            sky[(objrow1+1), :] = 0.
            obj[skyrow0:skyrow1, :] = 0.
            obj[(objrow1+1), :] = 0.
        else:
            self.logger.warning("non-standard nod-and-shuffle configuration")
            skyscl = -1.
            while skyscl < 0.:
                if self.config.instrument.plot_level >= 2:
                    q = input("Enter scale factor for sky to match obj "
                              "(float): ")
                    try:
                        skyscl = float(q)
                    except ValueError:
                        self.logger.warning("Invalid input: %s, try again" % q)
                        skyscl = -1.0
                else:
                    skyscl = 1.0
            self.logger.info("Sky scaling used = %.2f" % skyscl)
            objstd = self.action.args.ccddata.uncertainty.array.copy()
            objmsk = self.action.args.ccddata.mask.copy()
            objflg = self.action.args.ccddata.flags.copy()
            # move object to sky position
            obj[skyrow0:skyrow1, :] = obj[objrow0:objrow1, :]
            objstd[skyrow0:skyrow1, :] = std[objrow0:objrow1, :]
            objmsk[skyrow0:skyrow1, :] = msk[objrow0:objrow1, :]
            objflg[skyrow0:skyrow1, :] = flg[objrow0:objrow1, :]
            # do subtraction
            sky *= skyscl
            self.action.args.ccddata.data = obj - sky
            self.action.args.ccddata.uncertainty.array = np.sqrt(std ** 2 +
                                                                 objstd ** 2)
            self.action.args.ccddata.mask += objmsk
            self.action.args.ccddata.flags |= objflg
            # clean images
            self.action.args.ccddata.data[objrow0:objrow1, :] = 0.
            self.action.args.ccddata.data[0:skyrow0, :] = 0.
            self.action.args.ccddata.uncertainty.array[objrow0:objrow1, :] = 0.
            self.action.args.ccddata.uncertainty.array[0:skyrow0, :] = 0.
            self.action.args.ccddata.mask[objrow0:objrow1, :] = 1
            self.action.args.ccddata.mask[0:skyrow0, :] = 1
            self.action.args.ccddata.flags[objrow0:objrow1, :] = 64
            self.action.args.ccddata.flags[0:skyrow0, :] = 64
            sky[objrow0:objrow1, :] = 0.
            sky[0:skyrow0, :] = 0.
            obj[objrow0:objrow1, :] = 0.
            obj[0:skyrow0, :] = 0.
            cmnt = 'Aborted nod-and-shuffle observations'
            objhdr['COMMENT'] = cmnt
            skyhdr['COMMENT'] = cmnt
            self.action.args.ccddata.header['COMMENT'] = cmnt
            skyhdr['NASSCL'] = (skyscl, 'Scale factor applied to sky panel')
            self.action.args.ccddata.header['NASSCL'] = (
                skyscl, 'Scale factor applied to sky panel')

        # log
        self.logger.info("nod-and-shuffle subtracted, rows (sky0, 1, obj0,1): "
                         "%d, %d, %d, %d" % (skyrow0, skyrow1,
                                             objrow0, objrow1))
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
        msname = ofn.split('.')[0] + '_' + target_type.lower() + '.fits'
        out_sky = CCDData(sky, meta=skyhdr, unit=u_out)
        kcwi_fits_writer(out_sky, output_file=msname,
                         output_dir=self.config.instrument.output_directory)
        # write out object image
        obname = ofn.split('.')[0] + '_obj.fits'
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
                                            suffix="intk")
        self.context.proctab.write_proctab()

        self.logger.info(log_string)

        return self.action.args
    # END: class SubtractSky()
