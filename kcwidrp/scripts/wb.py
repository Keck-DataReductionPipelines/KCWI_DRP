#!/usr/bin/env python
import sys
if sys.version_info[0] > 2:
    from astropy.io import fits as pf
else:
    import pyfits as pf


def wb_main():
    """
    Generate summary log and group list files for BLUE channel images.

    Call get_log_string to create summary entry for each image and group them
    according to processing group.  Write out unique processing group lists
    in \*.txt files.  These files can be input to the pipeline with the -l
    command line parameter to allow processing of groups one at a time.  For
    example, 2x2 Blue biases taken with the TUP amp configuration in slow
    readout with gainmul 10 will end up in the file bias2x2TUP010_0.txt.  A
    master bias can be created by issuing the following command:

        >>> reduce_kcwi -b -l bias2x2TUP010_0.txt

    These group files are generated for biases, darks, continuum bars, arcs,
    flats, and all objects.  The filenames are all appended with the last four
    characters of the STATEID header keyword, so identical configurations from
    different states can be distinguished.

    Always good to type out the list file before processing it.

    Examples:

        >>> wb kb*.fits > whatb.list

        This will generate a summary log file along with associated group list
        files that can be used as inputs to the `reduce_kcwi` command with the
        `-l` parameter.  An example of the resulting \*.txt files is below::

            SN2023ixf2x2MedKBlueBL4500_75fe.txt      bias2x2TUP010_0.txt
            allb.txt                                 cbars2x2MedKBlueBL_4500_0.7_75fe.txt
            arcs2x2MedKBlueBLFeAr4500_10.0_75fe.txt  cflat2x2MedKBlueBL_4500_0.7_75fe.txt
            arcs2x2MedKBlueBLThAr4500_20.0_75fe.txt  dflat2x2MedKBlueBL_4500_14.0_75fe.txt
            bd26d26062x2MedKBlueBL4500_75fe.txt

        One can proceed through processing steps like this:

        >>> reduce_kcwi -b -l bias2x2TUP010_0.txt
        >>> reduce_kcwi -b -l cbars2x2MedKBlueBL_4500_0.7_75fe.txt
        >>> reduce_kcwi -b -l arcs2x2MedKBlueBLThAr4500_20.0_75fe.txt
        >>> reduce_kcwi -b -l cflat2x2MedKBlueBL_4500_0.7_75fe.txt
        >>> reduce_kcwi -b -l bd26d26062x2MedKBlueBL4500_75fe.txt
        >>> reduce_kcwi -b -l SN2023ixf2x2MedKBlueBL4500_75fe.txt

    """
    import sys

    if len(sys.argv) < 2:
        print("Usage - wb <fspec>")
    else:
        configs = []
        fnames = {"allb": []}
        for ifl in sys.argv[1:]:
            logstr, cfgstr, lsfn = get_log_string(ifl, batch=True)
            print(logstr)
            fnames['allb'].append(ifl)
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


def get_cal_list_file(hdr):
    """
    Return list file name given configuration in header.

    Generates the group file name apropriate for the input image header.

    Args:
        hdr (FITS header): the current image header

    :returns:
        (str): group filename

    :meta private:
    """

    fpre = {"ARCLAMP": "arcs", "CONTBARS": "cbars", "FLATLAMP": "cflat",
            "DOMEFLAT": "dflat", "TWIFLAT": "tflat"}

    imtype = hdr['IMTYPE']
    lfname = None
    if imtype in fpre:
        lfname = fpre[imtype]
        lfname += hdr['BINNING'].replace(',', 'x')
        lfname += hdr['IFUNAM'][:3]
        lfname += hdr['BFILTNAM']
        lfname += hdr['BGRATNAM']
        if 'arcs' in lfname:
            if 'FeAr' in hdr['ILLUME']:
                lfname += 'FeAr'
            else:
                lfname += 'ThAr'
        else:
            lfname += '_'
        lfname += "%.0f" % hdr['BCWAVE']
        lfname += "_%.1f" % hdr['EXPTIME']
        lfname += "_" + hdr['CONFIGID']
    elif 'BIAS' in imtype:
        lfname = 'bias'
        lfname += hdr['BINNING'].replace(',', 'x')
        lfname += hdr['AMPMODE']
        lfname += str(hdr['CCDMODE'])
        lfname += str(hdr['GAINMUL'])
        lfname += "_" + hdr['CONFIGID']
    elif 'DARK' in imtype:
        lfname = 'dark'
        lfname += hdr['BINNING'].replace(',', 'x')
        lfname += hdr['AMPMODE']
        lfname += str(hdr['CCDMODE'])
        lfname += str(hdr['GAINMUL'])
        lfname += "_%.1f" % hdr['EXPTIME']
        lfname += "_" + hdr['CONFIGID']
    elif 'OBJECT' in imtype:
        lfname = hdr['TARGNAME']
        lfname += hdr['BINNING'].replace(',', 'x')
        lfname += hdr['IFUNAM'][:3]
        lfname += hdr['BFILTNAM']
        lfname += hdr['BGRATNAM']
        lfname += "%.0f" % hdr['BCWAVE']
        lfname += "_" + hdr['CONFIGID']
    return lfname


def get_log_string(ifile, batch=False):
    """
    Generate log entry from BLUE FITS header keywords.

    Attempt to encapsulate the instrument configuration for each image by
    summarizing FITS header keyword values tersely.

    Args:
        ifile (str): filename of FITS image to summarize
        batch (bool): set to ``True`` for an abreviated record.  Defaults to ``False``.

    :returns:
        (str): Configuration summary string for the input FITS image file.

    :meta private:
    """

    try:
        ff = pf.open(ifile)
    except IOError:
        print("***ERROR*** empty or corrupt fits file: %s" % ifile)
        return None, None

    header = ff[0].header
    header['FNAME'] = ifile
    if 'CAMERA' in header:
        if 'BLUE' in header['CAMERA'].upper():
            is_bias = False
            if 'OFNAME' not in header:
                header['OFNAME'] = ifile
            if 'STATEID' not in header:
                header['CONFIGID'] = '----'
            else:
                header['CONFIGID'] = header['STATEID'][-4:]
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
            if 'CALMNAM' not in header:
                header['CALMNAM'] = '-'
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
            if header['TELAPSE'] > header['XPOSURE']:
                header['EXPTIME'] = header['TELAPSE']
            else:
                header['EXPTIME'] = header['XPOSURE']
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
                if 'DARK' in header['IMTYPE']:
                    if is_bias:
                        header['IMTYPE'] = 'BIAS'
                if not batch:
                    if 'object' not in header['CALTYPE']:
                        header['OBJECT'] = header['OBJECT'] + header['ILLUME']
            except KeyError:
                pass
            try:
                lstring = "%(OFNAME)19s %(CONFIGID)4s (%(AMPMODE)3s/%(BINNING)3s/%(CCDMODE)1d/" \
                          "%(GAINMUL)2d/%(NUMOPEN)2d/%(EXPTIME)6.1f s), (%(IFUNAM)3s/" \
                          "%(BFILTNAM)5s/%(BGRATNAM)4s/%(BGROTNAM)9s dg/" \
                          "%(BCWAVE)6.1f/%(CALMNAM)s/%(CALPNAM)5s/%(CALLANG)5.1f dg), " \
                          "(%(BARTANG)5.1f/%(BNASNAM)4s/%(BFOCMM)6.3f) %(AIRMASS)5.3f: %(IMTYPE)7s/" \
                          "%(ILLUME)6s/%(TARGNAME)s:%(OBJECT)s" % header
            except:
                lstring = "%19s : ?" % ifile

            if header['EXPTIME'] <= 0.0:
                cstr = "%(BINNING)3s:%(AMPMODE)3s:%(CCDMODE)1d:%(GAINMUL)2d:BIAS" % \
                       header
            else:
                if batch:
                    cstr = "%(BINNING)3s:%(BFILTNAM)s:%(BGRATNAM)s:%(IFUNAM)s:%(BCWAVE).1f" \
                           % header
                else:
                    cstr = "%(BINNING)3s:%(BFILTNAM)s:%(BGRATNAM)s:%(IFUNAM)s:%(BCWAVE).1f:" \
                           "%(EXPTIME)6.1f:%(OBJECT)s" % header
            lfn = get_cal_list_file(header)
        else:
            lstring = "%19s : NOT a BLUE image!" % ifile
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
    wb_main()
