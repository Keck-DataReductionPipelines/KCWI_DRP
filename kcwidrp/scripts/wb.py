#!/usr/bin/env python
from astropy.io import fits as pf


def get_log_string(ifile, batch=False):
    try:
        ff = pf.open(ifile)
    except IOError:
        print("***ERROR*** empty or corrupt fits file: %s" % ifile)
        quit()

    header = ff[0].header
    header['FNAME'] = ifile
    if 'CAMERA' in header:
        if 'BLUE' in header['CAMERA']:
            if 'OFNAME' not in header:
                header['OFNAME'] = ifile
            if 'AMPMODE' not in header:
                header['AMPMODE'] = '-'
            if 'BINNING' not in header:
                header['BINNING'] = '-'
            if 'CCDMODE' not in header:
                header['CCDMODE'] = -1
            if 'GAINMUL' not in header:
                header['GAINMUL'] = -1
            if 'NUMOPEN' not in header:
                header['NUMOPEN'] = -1
            if 'XPOSURE' not in header:
                header['XPOSURE'] = -1.
            if 'AIRMASS' not in header:
                header['AIRMASS'] = -1.
            if 'TELAPSE' not in header:
                header['TELAPSE'] = -1.
            if 'BFILTNAM' not in header:
                header['BFILTNAM'] = '-'
            if 'BGRATNAM' not in header:
                header['BGRATNAM'] = '-'
            if 'BGROTNAM' not in header:
                header['BGROTNAM'] = '-'
            if 'BCWAVE' not in header:
                header['BCWAVE'] = -1.
            if 'CALPNAM' not in header:
                header['CALPNAM'] = '-'
            if 'CALLANG' not in header:
                header['CALLANG'] = -1.
            if 'BARTANG' not in header:
                header['BARTANG'] = -1.
            if 'BNASNAM' not in header:
                header['BNASNAM'] = '-'
            if 'BFOCMM' not in header:
                header['BFOCMM'] = -1.
            if 'CALTYPE' not in header:
                header['CALTYPE'] = '-'
            if 'IFUNAM' not in header:
                header['IFUNAM'] = '-'
            else:
                header['IFUNAM'] = header['IFUNAM'][:3]
            if 'OBJECT' not in header:
                header['OBJECT'] = '-'
            if 'TARGNAME' not in header:
                header['TARGNAME'] = '-'
            if 'CALXNAM' not in header:
                header['CALXNAM'] = '-'
            if 'object' not in header['CALTYPE']:
                header['OBJECT'] = header['CALXNAM']
            else:
                header['OBJECT'] = header['TARGNAME']
            if header['TELAPSE'] > header['XPOSURE']:
                header['EXPTIME'] = header['TELAPSE']
            else:
                header['EXPTIME'] = header['XPOSURE']
            header['ILLUME'] = '-'
            if header['LMP0STAT'] == 1:
                if header['LMP0SHST'] == 1:
                    header['ILLUME'] = header['LMP0NAM']
            if header['LMP1STAT'] == 1:
                if header['LMP1SHST'] == 1:
                    header['ILLUME'] = header['LMP1NAM']
            if header['LMP2STAT'] == 1:
                if header['LMP2SHST'] == 1:
                    header['ILLUME'] = header['LMP2NAM']
            if header['LMP3STAT'] == 1:
                header['ILLUME'] = header['LMP3NAM'][0:6]
            if 'object' in header['CALTYPE']:
                if 'on' in header['FLSPECTR'] or 'on' in header['FLIMAGIN']:
                    header['ILLUME'] = 'DOME'
                    if not batch:
                        header['OBJECT'] = 'DOME'
            if not batch:
                if 'object' not in header['CALTYPE']:
                    header['OBJECT'] = header['OBJECT'] + header['ILLUME']
            try:
                lstring = "%(OFNAME)19s (%(AMPMODE)3s/%(BINNING)3s/%(CCDMODE)1d/" \
                          "%(GAINMUL)2d/%(NUMOPEN)2d/%(EXPTIME)6.1f s), (%(IFUNAM)3s/" \
                          "%(BFILTNAM)5s/%(BGRATNAM)4s/%(BGROTNAM)9s dg/" \
                          "%(BCWAVE)6.1f/%(CALPNAM)5s/%(CALLANG)5.1f dg), " \
                          "(%(BARTANG)5.1f/%(BNASNAM)4s/%(BFOCMM)6.3f) %(AIRMASS)5.3f: %(IMTYPE)7s/" \
                          "%(ILLUME)6s/%(OBJECT)s" % header
            except:
                lstring = "%28s : ?" % ifile

            if header['EXPTIME'] <= 0.0:
                cstr = "%(BINNING)3s:%(AMPMODE)3s:%(CCDMODE)1d:%(GAINMUL)2d:BIAS" % \
                       header
            else:
                if batch:
                    cstr = "%(BINNING)3s:%(BGRATNAM)s:%(IFUNAM)s:%(BCWAVE).1f" \
                           % header
                else:
                    cstr = "%(BINNING)3s:%(BGRATNAM)s:%(IFUNAM)s:%(BCWAVE).1f:" \
                           "%(EXPTIME)6.1f:%(OBJECT)s" % header
        else:
            lstring = ifile + ' FPC image'
            cstr = None

        return lstring, cstr
    else:
        print("ERROR - Camera can not be determined.")


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        print("Usage - wb <fspec>")
    else:
        configs = []
        for ifl in sys.argv[1:]:
            logstr, cfgstr = get_log_string(ifl, batch=True)
            print(logstr)
            if cfgstr:
                configs.append(cfgstr)

        # Unique configs
        uconfigs = sorted(set(configs))
        print("Number of unique configurations = %d" % len(uconfigs))
        for c in uconfigs:
            print(c)
