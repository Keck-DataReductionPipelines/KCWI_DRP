# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the misc directory in idlutils.
"""
import numpy as np
from . import PydlutilsException


def decode_mixed(x):
    """Convert bytes in Numpy arrays into strings.  Leave other stuff alone.

    Parameters
    ----------
    x : object
        Input object.

    Returns
    -------
    object
        If `x` has a ``decode()`` method, ``x.decode()`` will be returned.
        Otherwise `x` will be returned unchanged.
    """
    try:
        return x.decode()
    except:
        return x


def djs_laxisgen(dims, iaxis=0):
    """Returns an integer array where each element of the array is set
    equal to its index number along the specified axis.

    Parameters
    ----------
    dims : :class:`list`
        Dimensions of the array to return.
    iaxis : :class:`int`, optional
        Index along this dimension.

    Returns
    -------
    :class:`numpy.ndarray`
        An array of indexes with ``dtype=int32``.

    Raises
    ------
    :exc:`ValueError`
        If `iaxis` is greater than or equal to the number of dimensions.

    Notes
    -----
    For two or more dimensions, there is no difference between this routine
    and :func:`~pydl.pydlutils.misc.djs_laxisnum`.

    Examples
    --------
    >>> from pydl.pydlutils.misc import djs_laxisgen
    >>> print(djs_laxisgen([4,4]))
    [[0 0 0 0]
     [1 1 1 1]
     [2 2 2 2]
     [3 3 3 3]]
    """
    ndimen = len(dims)
    if ndimen == 1:
        return np.arange(dims[0], dtype='i4')
    return djs_laxisnum(dims, iaxis)


def djs_laxisnum(dims, iaxis=0):
    """Returns an integer array where each element of the array is set equal
    to its index number in the specified axis.

    Parameters
    ----------
    dims : :class:`list`
        Dimensions of the array to return.
    iaxis : :class:`int`, optional
        Index along this dimension.

    Returns
    -------
    :class:`numpy.ndarray`
        An array of indexes with ``dtype=int32``.

    Raises
    ------
    :exc:`ValueError`
        If `iaxis` is greater than or equal to the number of dimensions, or
        if number of dimensions is greater than three.

    Notes
    -----
    For two or more dimensions, there is no difference between this routine
    and :func:`~pydl.pydlutils.misc.djs_laxisgen`.

    Examples
    --------
    >>> from pydl.pydlutils.misc import djs_laxisnum
    >>> print(djs_laxisnum([4,4]))
    [[0 0 0 0]
     [1 1 1 1]
     [2 2 2 2]
     [3 3 3 3]]
    """
    ndimen = len(dims)
    result = np.zeros(dims, dtype='i4')
    if ndimen == 1:
        pass
    elif ndimen == 2:
        if iaxis == 0:
            for k in range(dims[0]):
                result[k, :] = k
        elif iaxis == 1:
            for k in range(dims[1]):
                result[:, k] = k
        else:
            raise ValueError("Bad value for iaxis: {0:d}".format(iaxis))
    elif ndimen == 3:
        if iaxis == 0:
            for k in range(dims[0]):
                result[k, :, :] = k
        elif iaxis == 1:
            for k in range(dims[1]):
                result[:, k, :] = k
        elif iaxis == 2:
            for k in range(dims[2]):
                result[:, :, k] = k
        else:
            raise ValueError("Bad value for iaxis: {0:d}".format(iaxis))
    else:
        raise ValueError("{0:d} dimensions not supported.".format(ndimen))
    return result


def struct_print(array, filename=None, formatcodes=None, alias=None,
                 fdigit=5, ddigit=7, html=False, no_head=False,
                 silent=False):
    """Print a NumPy record array (analogous to an IDL structure) in a
    nice way.

    Parameters
    ----------
    array : :class:`numpy.ndarray`
        A record array to print.
    filename : :class:`str` or file-like, optional
        If supplied, write to this file.
    formatcodes : :class:`dict`, optional
        If supplied, use explicit format for certain columns.
    alias : :class:`dict`, optional
        If supplied, use this mapping of record array column names to printed
        column names.
    fdigit : :class:`int`, optional
        Width of 32-bit floating point columns, default 5.
    ddigit : :class:`int`, optional
        Width of 64-bit floating point columns, default 7.
    html : :class:`bool`, optional
        If ``True``, print an html table.
    no_head : :class:`bool`, optional
        If ``True``, *don't* print a header line.
    silent : :class:`bool`, optional
        If ``True``, do not print the table, just return it.

    Returns
    -------
    :func:`tuple`
        A tuple containing a list of the lines in the table.  If `html` is
        ``True``, also returns a list of lines of CSS for formatting the
        table.

    Examples
    --------
    import numpy as np
    from pydl.pydlutils.misc import struct_print
    struct_print(np.array([(1,2.34,'five'),(2,3.456,'seven'),(3,4.5678,'nine')],
        dtype=[('a','i4'),('bb','f4'),('ccc','S5')]),silent=True)
    (['a bb          ccc  ', '- ----------- -----', '1        2.34 five ',
         '2       3.456 seven', '3      4.5678 nine '], [])
    """
    if html:
        headstart = '<tr><th>'
        headsep = '</th><th>'
        headend = '</th></tr>'
        colstart = '<tr><td>'
        colsep = '</td><td>'
        colend = '</td></tr>'
        css = ['<style type="text/css">',
               'table {',
               '    border-collapse: collapse;',
               '}',
               'th {',
               '    padding: 2px;',
               '    text-align: right;',
               '    border: 1px solid black;',
               '    font-weight: bold;',
               '}',
               'td {',
               '    padding: 2px;',
               '    text-align: right;',
               '    border: 1px solid black;',
               '}',
               '</style>']
    else:
        headstart = ''
        headsep = ' '
        headend = ''
        colstart = ''
        colsep = ' '
        colend = ''
        css = list()
    #
    # Alias should be a dictionary that maps structure names to column names
    #
    if alias is None:
        #
        # Create a dummy alias dictionary
        #
        alias = dict(list(zip(array.dtype.names, array.dtype.names)))
    else:
        #
        # Fill in any missing values of the alias dictionary
        #
        for tag in array.dtype.names:
            if tag not in alias:
                alias[tag] = tag
    #
    # Formatcodes allows an override for certain columns.
    #
    if formatcodes is None:
        formatcodes = dict()
    #
    # This dictionary will hold the number of characters in each column
    #
    nchar = dict()
    #
    # Construct format codes for each column
    #
    for k, tag in enumerate(array.dtype.names):
        if tag in formatcodes:
            thiscode = formatcodes[tag]
            thisn = len(thiscode.format(array[tag][0]))
        else:
            d = array.dtype.fields[tag][0]
            if d.kind == 'i' or d.kind == 'u':
                thisn = max(max(len(str(array[tag].min())),
                            len(str(array[tag].max()))), len(tag))
                thiscode = "{{{0:d}:{1:d}d}}".format(k, thisn)
            elif d.kind == 'f':
                if d.itemsize == 8:
                    prec = ddigit
                else:
                    prec = fdigit
                thisn = prec + 6
                if array[tag].min() < 0:
                    thisn += 1
                thiscode = "{{{0:d}:{1:d}.{2:d}g}}".format(k, thisn, prec)
            elif d.kind == 'S' or d.kind == 'U':
                thisn = max(d.itemsize, len(tag))
                thiscode = "{{{0:d}:{1:d}s}}".format(k, thisn)
            else:
                raise PydlutilsException(
                        "Unsupported kind: {0}".format(d.kind))
            formatcodes[tag] = thiscode
        nchar[tag] = thisn
    #
    # Start building an array of lines
    #
    lines = list()
    #
    # Construct header lines
    #
    if html:
        lines.append('<table>')
        hdr1 = (headstart + headsep.join([alias[tag]
                for tag in array.dtype.names]) + headend)
        lines.append(hdr1)
    else:
        if not no_head:
            hdr1 = (headstart + headsep.join([("{{0:{0:d}s}}".format(
                    nchar[tag])).format(alias[tag])
                    for tag in array.dtype.names]) + headend)
            hdr2 = (headstart + headsep.join(['-' * nchar[tag]
                    for tag in array.dtype.names]) + headend)
            lines.append(hdr1)
            lines.append(hdr2)
    #
    # Create a format string for the data from the individual format codes
    #
    rowformat = (colstart + colsep.join([formatcodes[tag]
                 for tag in array.dtype.names]) + colend)
    for k in range(array.size):
        lines.append(rowformat.format(
                    *([decode_mixed(i) for i in array[k].tolist()])))
    if html:
        lines.append('</table>')
    f = None   # This variable will store a file handle
    close_file = False
    if filename is not None:
        if hasattr(filename, 'write'):
            f = filename
        else:
            f = open(filename, 'w+b')
            close_file = True
    if f is None:
        if not silent:  # pragma: no cover
            print("\n".join(lines)+"\n")
    else:
        f.write(("\n".join(lines)+"\n").encode('utf-8'))
        if close_file:
            f.close()
    return lines, css
