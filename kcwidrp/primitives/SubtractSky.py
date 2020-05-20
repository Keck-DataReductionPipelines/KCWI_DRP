from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer
import os


class SubtractSky(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if a master sky exists to subtract
        :return:
        """
        self.logger.info("Checking precondition for SubtractSky")
        skyfile = None
        skymask = None
        # check if kcwi.sky exists
        if os.path.exists('kcwi.sky'):
            f = open('kcwi.sky')
            skyproc = f.readlines()
            f.close()
            # is our file in the list?
            ofn = self.action.args.ccddata.header['OFNAME']
            for row in skyproc:
                if ofn in row.split()[0]:
                    skyfile = row.split()[1]
                    if len(row.split()) > 2:
                        skymask = row.split()[2]
            if skyfile:
                if not os.path.exists(skyfile):
                    skyfile = None
            if skymask:
                if not os.path.exists(skymask):
                    skymask = None
        self.action.args.skyfile = skyfile
        self.action.args.skymask = skymask
        if skyfile:
            self.logger.info("pre condition got 1 master sky, expected 1")
            return True
        else:

            target_type = 'SKY'
            tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                                 target_type=target_type,
                                                 nearest=True)
            self.logger.info("pre condition got %d master sky, expected 1"
                             % len(tab))
            if len(tab) <= 0:
                return False
            else:
                return True

    def _perform(self):
        self.logger.info("Subtracting sky background")

        # Header keyword to update
        key = 'SKYCOR'
        keycom = 'sky corrected?'
        target_type = 'SKY'

        skyfile = self.action.args.skyfile
        skymask = self.action.args.skymask

        if not self.action.args.skyfile:
            tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                                 target_type=target_type,
                                                 nearest=True)
            self.logger.info("%d master sky frames found" % len(tab))

            if len(tab) > 0:
                skyfile = tab['OFNAME'][0]

        msname = skyfile.split('.')[0] + '_' + target_type.lower() + ".fits"
        if os.path.exists(os.path.join(os.path.dirname(self.action.args.name),
                                       'redux', msname)):
            self.logger.info("Reading image: %s" % msname)
            msky = kcwi_fits_reader(
                os.path.join(os.path.dirname(self.action.args.name), 'redux',
                             msname))[0]

            # scale the sky?
            obtime = self.action.args.ccddata.header['XPOSURE']
            sktime = msky.header['XPOSURE']

            if obtime <= 0. or sktime <= 0.:
                self.logger.warning("Bad exposure times (obj, sky): %.1f, %1f"
                                    % (obtime, sktime))
                skscl = 1.
            else:
                skscl = obtime / sktime
            self.logger.info("Sky scale factor = %.3f" % skscl)

            # do the subtraction
            self.action.args.ccddata.data -= msky.data * skscl

            # update header keywords
            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['SKYMAST'] = (msname,
                                                          "Master sky filename")
            self.action.args.ccddata.header['SKYSCL'] = (skscl,
                                                         'sky scale factor')
            if skymask:
                self.action.args.ccddata.header['SKYMSKF'] = (skymask,
                                                              'sky mask file')

        else:
            # update header keywords
            self.action.args.ccddata.header[key] = (False, keycom)

        log_string = SubtractSky.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string

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
