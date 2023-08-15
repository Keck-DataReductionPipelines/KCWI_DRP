#!/usr/bin/env python
import sys
if sys.version_info[0] > 2:
    from astropy.io import fits as pf
else:
    import pyfits as pf


def get_cal_list_file(hdr):
    """Return list file name given configuration in header"""

    fpre = {"ARCLAMP": "arcs", "CONTBARS": "cbars", "FLATLAMP": "cflat",
            "DOMEFLAT": "dflat", "TWIFLAT": "tflat"}

    imtype = hdr['IMTYPE']
    lfname = None
    if imtype in fpre:
        lfname = fpre[imtype]
        lfname += hdr['BINNING'].replace(',', 'x')
        lfname += hdr['IFUNAM'][:3]
        lfname += hdr['RGRATNAM']
        if 'arcs' in lfname:
            if 'FeAr' in hdr['ILLUME']:
                lfname += 'FeAr'
            else:
                lfname += 'ThAr'
        else:
            lfname += '_'
        lfname += "%.0f" % hdr['RCWAVE']
        lfname += "_%.1f" % hdr['EXPTIME']
    elif 'BIAS' in imtype:
        lfname = 'bias'
        lfname += hdr['BINNING'].replace(',', 'x')
        lfname += hdr['AMPMODE']
        lfname += str(hdr['CDSSPEED'])
        lfname += str(hdr['ADCGAINS'])
    elif 'OBJECT' in imtype:
        lfname = hdr['TARGNAME']
        lfname += hdr['BINNING'].replace(',', 'x')
        lfname += hdr['IFUNAM'][:3]
        lfname += hdr['RGRATNAM']
        lfname += "%.0f" % hdr['RCWAVE']
    return lfname


def get_log_string(ifile, batch=False):
    try:
        ff = pf.open(ifile)
    except IOError:
        print("***ERROR*** empty or corrupt fits file: %s" % ifile)
        return None, None

    header = ff[0].header
    header['FNAME'] = ifile
    if 'CAMERA' in header:
        if 'RED' in header['CAMERA'].upper():
            is_bias = False
            if 'OFNAME' not in header:
                header['OFNAME'] = ifile
            if 'AMPMODE' not in header:
                header['AMPMODE'] = '-'
            if 'BINNING' not in header:
                header['BINNING'] = '-'
            if 'CDSSPEED' not in header:
                header['CDSSPEED'] = -1
            if 'ADCGAINS' not in header:
                header['ADCGAINS'] = -1
            if 'NUMOPEN' not in header:
                header['NUMOPEN'] = -1
            if 'XPOSURE' not in header:
                header['XPOSURE'] = -1.
            if 'AIRMASS' not in header:
                header['AIRMASS'] = -1.
            if 'TELAPSE' not in header:
                header['TELAPSE'] = -1.
            if 'RFILTNAM' not in header:
                header['RFILTNAM'] = '-'
            if 'RGRATNAM' not in header:
                header['RGRATNAM'] = '-'
            if 'RGROTNAM' not in header:
                header['RGROTNAM'] = '-'
            if 'RCWAVE' not in header:
                header['RCWAVE'] = -1.
            if 'CALMNAM' not in header:
                header['CALMNAM'] = '-'
            if 'CALPNAM' not in header:
                header['CALPNAM'] = '-'
            if 'CALLANG' not in header:
                header['CALLANG'] = -1.
            if 'RARTANG' not in header:
                header['RARTANG'] = -1.
            if 'RNASNAM' not in header:
                header['RNASNAM'] = '-'
            if 'RFOCMM' not in header:
                header['RFOCMM'] = -1.
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
            # if header['TELAPSE'] > header['XPOSURE']:
            #    header['EXPTIME'] = header['TELAPSE']
            # else:
            nshuf = header['NSHFUP']
            ttime = header['TTIME']
            if nshuf > 0:
                header['EXPTIME'] = (nshuf * ttime, "N&S frame exptime")
            else:
                header['EXPTIME'] = (header['XPOSURE'], "Shutter open time")
            if header['EXPTIME'] <= 0.:
                is_bias = True
            header['ILLUME'] = '-'
            try:
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
                if 'BIAS' in header['IMTYPE']:
                    if not is_bias:
                        header['IMTYPE'] = 'DARK'
                if not batch:
                    if 'object' not in header['CALTYPE']:
                        header['OBJECT'] = header['OBJECT'] + header['ILLUME']
            except KeyError:
                pass
            try:
                lstring = "%(OFNAME)19s (%(AMPMODE)8s/%(BINNING)3s/%(CDSSPEED)1d/" \
                          "%(ADCGAINS)1d/%(NUMOPEN)2d/%(EXPTIME)6.1f s), (%(IFUNAM)3s/" \
                          "%(RFILTNAM)5s/%(RGRATNAM)4s/%(RGROTNAM)9s dg/" \
                          "%(RCWAVE)6.1f/%(CALMNAM)5s/%(CALPNAM)5s/%(CALLANG)5.1f dg), " \
                          "(%(RARTANG)5.1f/%(RNASNAM)4s/%(RFOCMM)6.3f) %(AIRMASS)5.3f: %(IMTYPE)7s/" \
                          "%(ILLUME)6s/%(TARGNAME)s:%(OBJECT)s" % header
            except:
                lstring = "%19s : ?" % ifile

            if header['EXPTIME'] <= 0.0:
                cstr = "%(BINNING)3s:%(AMPMODE)8s:%(CDSSPEED)1d:%(ADCGAINS)1d:BIAS" % \
                       header
            else:
                if batch:
                    cstr = "%(BINNING)3s:%(RGRATNAM)s:%(IFUNAM)s:%(RCWAVE).1f" \
                           % header
                else:
                    cstr = "%(BINNING)3s:%(RGRATNAM)s:%(IFUNAM)s:%(RCWAVE).1f:" \
                           "%(EXPTIME)6.1f:%(OBJECT)s" % header
            lfn = get_cal_list_file(header)
        else:
            lstring = "%19s : NOT a RED image!" % ifile
            cstr = None
            lfn = None

    else:
        if not batch:
            print("ERROR - Camera can not be determined.")
        lstring = "%19s : No CAMERA keyword!" % ifile
        cstr = None
        lfn = None

    return lstring, cstr, lfn


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        print("Usage - wr <fspec>")
    else:
        configs = []
        fnames = {"allr": []}
        for ifl in sys.argv[1:]:
            logstr, cfgstr, lsfn = get_log_string(ifl, batch=True)
            print(logstr)
            fnames['allr'].append(ifl)
            if lsfn:
                if lsfn in fnames:
                    fnames[lsfn].append(ifl)
                else:
                    fnames[lsfn] = [ifl]
            if cfgstr:
                configs.append(cfgstr)
        # Unique configs
        uconfigs = sorted(set(configs))
        print("Number of unique configurations = %d" % len(uconfigs))
        for c in uconfigs:
            print(c)

        for cal in fnames:
            with open(cal+".txt", 'w') as cal_list:
                for f in fnames[cal]:
                    cal_list.write(f + "\n")
