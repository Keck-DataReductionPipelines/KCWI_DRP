from astropy.table import Table, unique
import os
import logging


class Proctab:

    def __init__(self, logger=None):
        self.log = logger if logger is not None else logging.getLogger('KCWI')
        self.proctab = None

    def new_proctab(self):
        cnames = ('FRAMENO', 'CID', 'DID', 'TYPE', 'GRPID', 'TTIME', 'CAM',
                  'IFU', 'GRAT', 'GANG', 'CWAVE', 'BIN', 'FILT', 'MJD',
                  'STAGE', 'SUFF', 'OFNAME', 'TARGNAME', 'filename')
        dtypes = ('int32', 'S24', 'int64', 'S9', 'S12', 'float64', 'S4',
                  'S6', 'S5', 'float64', 'float64', 'S4', 'S5', 'float64',
                  'int32', 'S5', 'S25', 'S25', 'S25')
        meta = {'KCWI DRP PROC TABLE': 'new table'}
        self.proctab = Table(names=cnames, dtype=dtypes, meta=meta)
        # prevent string column truncation
        for col in self.proctab.itercols():
            if col.dtype.kind in 'SU':
                self.proctab.replace_column(col.name, col.astype('object'))

    def read_proctab(self, tfil='kcwi.proc'):
        if os.path.isfile(tfil):
            self.log.info("reading proc table file: %s" % tfil)
            self.proctab = Table.read(tfil, format='ascii.fixed_width')
            self.proctab.dtypes = ('int32', 'S24', 'int64', 'S9', 'S12',
                                   'float64', 'S4', 'S6', 'S5', 'float64',
                                   'float64', 'S4', 'S5', 'float64', 'int32',
                                   'S5', 'S25', 'S25', 'S25')
        else:
            self.log.info("proc table file not found: %s" % tfil)
            self.new_proctab()
        if len(self.proctab) == 0:
            self.new_proctab()
        # format columns
        self.proctab['GANG'].format = '7.2f'
        self.proctab['CWAVE'].format = '8.2f'
        self.proctab['MJD'].format = '15.6f'
        # prevent string column truncation
        for col in self.proctab.itercols():
            if col.dtype.kind in 'SU':
                self.proctab.replace_column(col.name, col.astype('object'))

    def write_proctab(self, tfil='kcwi.proc'):
        if self.proctab is not None:
            if os.path.isfile(tfil):
                over_write = True
            else:
                over_write = False

            self.proctab.write(filename=tfil, format='ascii.fixed_width',
                               overwrite=over_write)
            self.log.info("writing proc table file: %s" % tfil)
        else:
            self.log.info("no proc table to write")

    def update_proctab(self, frame, suffix='raw', newtype=None, filename=""):
        if filename == "":
            self.log.error(f"No filename given for {frame.header['OFNAME']}")
        if frame is not None and self.proctab is not None:
            stages = {'RAW': 0,
                      'mbias': 1,
                      'int': 1,
                      'intd': 2,
                      'mcbars': 3,
                      'marc': 3,
                      'mflat': 4,
                      'sflat': 4,
                      'intf': 4,
                      'sky': 5,
                      'intk': 5,
                      'icube': 6,
                      'icubed': 7,
                      'invsens': 8,
                      'icubes': 8}
            if suffix in stages:
                stage = stages[suffix]
            else:
                stage = 9
            if newtype is not None:
                frame.header['IMTYPE'] = newtype
            # new row for proc table
            if frame.header['STATEID'].strip() == '0':
                frame.header['STATEID'] = 'NONE'
            if 'GROUPID' not in frame.header:
                frame.header['GROUPID'] = "NONE"
            else:
                grpid = frame.header['GROUPID'].strip()
                if len(grpid) <= 0:
                    frame.header['GROUPID'] = "NONE"
            #    dto = frame.header['DATE-OBS']
            #    fno = frame.header['FRAMENO']
            #    frame.header['GROUPID'] = "%s-%s" % (dto, fno)
            cam = frame.header['CAMERA'].upper()
            if 'BLUE' in cam:
                grnam = frame.header['BGRATNAM']
                grang = frame.header['BGRANGLE']
                cwave = frame.header['BCWAVE']
                fltnm = frame.header['BFILTNAM']
            elif 'RED' in cam:
                grnam = frame.header['RGRATNAM']
                grang = frame.header['RGRANGLE']
                cwave = frame.header['RCWAVE']
                fltnm = None
            else:
                grnam = None
                grang = None
                cwave = None
                fltnm = None
            trgnm = frame.header['TARGNAME'].replace(" ", "")
            if len(trgnm) <= 0:
                trgnm = frame.header['OBJECT'].replace(" ", "")
            new_row = [frame.header['FRAMENO'],
                       frame.header['STATEID'],
                       frame.header['CCDCFG'],
                       frame.header['IMTYPE'],
                       frame.header['GROUPID'],
                       frame.header['TTIME'],
                       cam,
                       frame.header['IFUNAM'],
                       grnam,
                       grang,
                       cwave,
                       frame.header['BINNING'],
                       fltnm,
                       frame.header['MJD'],
                       stage,
                       suffix,
                       frame.header['OFNAME'],
                       trgnm,
                       filename]
        else:
            new_row = None
        self.proctab.add_row(new_row)
        self.proctab = unique(self.proctab, keys=['CID', 'FRAMENO', 'TYPE'],
                              keep='last')
        self.proctab.sort('FRAMENO')
        self.log.info(
            f"proctable updated with {frame.header['OFNAME']} and {filename}")

    def search_proctab(self, frame, target_type=None, target_group=None,
                       nearest=False):
        if target_type is not None and self.proctab is not None:
            self.log.info('Looking for %s frames' % target_type)
            # get relevant camera (blue or red)
            self.log.info('Camera is %s' % frame.header['CAMERA'])
            tab = self.proctab[(self.proctab['CAM'] ==
                                frame.header['CAMERA'].upper().strip())]
            # get target type images
            tab = tab[(self.proctab['TYPE'] == target_type)]
            self.log.info('Target type is %s' % target_type)
            # BIASES must have the same CCDCFG
            if 'BIAS' in target_type:
                self.log.info('Looking for frames with CCDCFG = %s' %
                              frame.header['CCDCFG'])
                tab = tab[(tab['DID'] == int(frame.header['CCDCFG']))]
                if target_group is not None:
                    self.log.info('Looking for frames with GRPID = %s' %
                                  target_group)
                    tab = tab[(tab['GRPID'] == target_group)]
            # raw DARKS must have the same CCDCFG and TTIME
            elif target_type == 'DARK':
                self.log.info('Looking for frames with CCDCFG = %s and '
                              'TTIME = %f' % (frame.header['CCDCFG'],
                                              frame.header['TTIME']))
                tab = tab[tab['DID'] == int(frame.header['CCDCFG'])]
                tab = tab[tab['TTIME'] == float(frame.header['TTIME'])]
                if target_group is not None:
                    self.log.info('Looking for frames with GRPID = %s' %
                                  target_group)
                    tab = tab[tab['GRPID'] == target_group]
            # MDARKS must have the same CCDCFG, will be scaled to match TTIME
            elif target_type == 'MDARK':
                self.log.info('Looking for frames with CCDCFG = %s' %
                              frame.header['CCDCFG'])
                tab = tab[(tab['DID'] == int(frame.header['CCDCFG']))]
            elif target_type == 'OBJECT':
                self.log.info('Looking for frames with GRPID = %s' %
                              target_group)
                tab = tab[tab['GRPID'] == target_group]
            else:
                self.log.info('Looking for frames with STATEID = %s (%s)' %
                              (frame.header['STATEID'],
                               frame.header['STATENAM']))
                tab = tab[(tab['CID'] == frame.header['STATEID'])]
            # Check if nearest entry is requested
            if nearest and len(tab) > 1:
                tfno = frame.header['MJD']
                minoff = 99999
                trow = None
                for row in tab:
                    off = abs(row['MJD'] - tfno)
                    if off < minoff:
                        minoff = off
                        trow = row
                if trow is not None:
                    tab = tab[(tab['MJD'] == trow['MJD'])]
        else:
            if target_type is None:
                self.log.warning("No target for proctab")
            if self.proctab is None:
                self.log.warning("Proctab is empty")
            tab = None
        return tab

    def in_proctab(self, frame):
        # get relevant camera (blue or red)
        tab = self.proctab[(self.proctab['CAM'] ==
                            frame.header['CAMERA'].upper().strip())]
        imno_list = tab['MJD']
        if frame.header['MJD'] in imno_list:
            return True
        else:
            return False

    def last_suffix(self, frame):
        # get last suffix if currently in proctab
        # check for master object first
        tab_mobj = self.proctab[(self.proctab['TYPE'] == 'MOBJ')]
        tab_mjd = tab_mobj[(tab_mobj['MJD'] == frame.header['MJD'])]
        tab = tab_mjd[(tab_mjd['FRAMENO'] == frame.header['FRAMENO'])]

        if len(tab) == 1:
            return tab['SUFF'].value[0]
        # No master object, so check simple object
        elif len(tab) <= 0:
            tab_obj = self.proctab[(self.proctab['TYPE'] == 'OBJECT')]
            tab_mjd = tab_obj[(tab_obj['MJD'] == frame.header['MJD'])]
            tab = tab_mjd[(tab_mjd['FRAMENO'] == frame.header['FRAMENO'])]
            if len(tab) == 1:
                return tab['SUFF'].value[0]
            else:
                return ""
        else:
            self.log.warning("Ambiguous entries for frame number %d" %
                             frame.header['FRAMENO'])
            return ""
