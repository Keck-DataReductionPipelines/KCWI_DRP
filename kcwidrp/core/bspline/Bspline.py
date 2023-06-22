# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the bspline directory in idlutils.
"""
from warnings import warn
import numpy as np
from numpy.linalg.linalg import LinAlgError
from scipy.linalg import cholesky_banded, cho_solve_banded
from . import PydlutilsUserWarning
from .math import djs_reject, flegendre
from .trace import fchebyshev
from .uniq import uniq


class Bspline(object):
    """B-spline class.

    Functions in the idlutils bspline library are implemented as methods on this
    class.

    Parameters
    ----------
    x : :class:`numpy.ndarray`
        The data.
    nord : :class:`int`, optional
        The order of the B-spline.  Default is 4, which is cubic.
    npoly : :class:`int`, optional
        Polynomial order to fit over 2nd variable, if supplied.  If not
        supplied the order is 1.
    bkpt : :class:`numpy.ndarray`, optional
        To be documented.
    bkspread : :class:`float`, optional
        To be documented.

    Attributes
    ----------
    breakpoints
        To be documented.
    nord
        To be documented.
    npoly
        To be documented.
    mask
        To be documented.
    coeff
        To be documented.
    icoeff
        To be documented.
    xmin
        To be documented.
    xmax
        To be documented.
    funcname
        To be documented.
    """

    def __init__(self, x, nord=4, npoly=1, bkpt=None, fullbkpt=None,
                 bkspread=1.0, verbose=False, **kwargs):
        """Init creates an object whose attributes are similar to the
        structure returned by the ``create_bsplineset()`` function.
        """
        #
        # Set the breakpoints.
        #
        if fullbkpt is None:
            if bkpt is None:
                startx = x.min()
                rangex = x.max() - startx
                if 'placed' in kwargs:
                    w = ((kwargs['placed'] >= startx) &
                         (kwargs['placed'] <= startx + rangex))
                    if w.sum() < 2:
                        bkpt = np.arange(2, dtype='f') * rangex + startx
                    else:
                        bkpt = kwargs['placed'][w]
                elif 'bkspace' in kwargs:
                    nbkpts = int(rangex / kwargs['bkspace']) + 1
                    if nbkpts < 2:
                        nbkpts = 2
                    tempbkspace = rangex / float(nbkpts - 1)
                    bkpt = np.arange(nbkpts, dtype='f') * tempbkspace + startx
                elif 'nbkpts' in kwargs:
                    nbkpts = kwargs['nbkpts']
                    if nbkpts < 2:
                        nbkpts = 2
                    tempbkspace = rangex / float(nbkpts - 1)
                    bkpt = np.arange(nbkpts, dtype='f') * tempbkspace + startx
                elif 'everyn' in kwargs:
                    nx = x.size
                    nbkpts = max(nx / kwargs['everyn'], 1)
                    if nbkpts == 1:
                        xspot = [0]
                    else:
                        xspot = (nx / nbkpts) * np.arange(nbkpts)
                    bkpt = np.interp(xspot, np.arange(nx), x)
                else:
                    raise ValueError('No information for bkpts.')
            imin = bkpt.argmin()
            imax = bkpt.argmax()
            if x.min() < bkpt[imin]:
                if verbose:
                    print('Lowest breakpoint does not cover lowest x value: '
                          'changing.')
                bkpt[imin] = x.min()
            if x.max() > bkpt[imax]:
                if verbose:
                    print('Highest breakpoint does not cover highest x value: '
                          'changing.')
                bkpt[imax] = x.max()
            nshortbkpt = bkpt.size
            fullbkpt = bkpt.copy()
            if nshortbkpt == 1:
                bkspace = np.float32(bkspread)
            else:
                bkspace = (bkpt[1] - bkpt[0]) * np.float32(bkspread)
            for i in np.arange(1, nord, dtype=np.float32):
                fullbkpt = np.insert(fullbkpt, 0, bkpt[0] - bkspace * i)
                fullbkpt = np.insert(fullbkpt, fullbkpt.shape[0],
                                     bkpt[nshortbkpt - 1] + bkspace * i)
        #
        # Set the attributes
        #
        nc = fullbkpt.size - nord
        self.breakpoints = fullbkpt.astype(float)
        self.nord = nord
        self.npoly = npoly
        self.mask = np.ones((fullbkpt.size,), dtype=bool)
        if npoly > 1:
            self.coeff = np.zeros((npoly, nc), dtype=float)
            self.icoeff = np.zeros((npoly, nc), dtype=float)
        else:
            self.coeff = np.zeros((nc,), dtype=float)
            self.icoeff = np.zeros((nc,), dtype=float)
        self.xmin = 0.0
        self.xmax = 1.0
        self.funcname = kwargs[
            'funcname'] if 'funcname' in kwargs else 'legendre'

        return

    def fit(self, xdata, ydata, invvar, x2=None):
        """Calculate a B-spline in the least-squares sense.

        Fit is based on two variables: `xdata` sorted and spans a large range
        where breakpoints are required `ydata` can be described with a low order
        polynomial.

        Parameters
        ----------
        xdata : :class:`numpy.ndarray`
            Independent variable.
        ydata : :class:`numpy.ndarray`
            Dependent variable.
        invvar : :class:`numpy.ndarray`
            Inverse variance of `ydata`.
        x2 : :class:`numpy.ndarray`, optional
            Orthogonal dependent variable for 2d fits.

        Returns
        -------
        :func:`tuple`
            A tuple containing an integer error code, and the evaluation of the
            b-spline at the input values.  An error code of -2 is a failure,
            -1 indicates dropped breakpoints, 0 is success, and positive
            integers indicate ill-conditioned breakpoints.
        """
        goodbk = self.mask[self.nord:]
        nn = goodbk.sum()
        if nn < self.nord:
            yfit = np.zeros(ydata.shape, dtype='f')
            return -2, yfit
        nfull = nn * self.npoly
        bw = self.npoly * self.nord
        a1, lower, upper = self.action(xdata, x2=x2)
        foo = np.tile(invvar, bw).reshape(bw, invvar.size).transpose()
        a2 = a1 * foo
        alpha = np.zeros((bw, nfull+bw), dtype='d')
        beta = np.zeros((nfull+bw,), dtype='d')
        bi = np.arange(bw, dtype='i4')
        bo = np.arange(bw, dtype='i4')
        for k in range(1, bw):
            bi = np.append(bi, np.arange(bw-k, dtype='i4')+(bw+1)*k)
            bo = np.append(bo, np.arange(bw-k, dtype='i4')+bw*k)
        for k in range(nn-self.nord+1):
            itop = k*self.npoly
            ibottom = min(itop, nfull) + bw - 1
            ict = upper[k] - lower[k] + 1
            if ict > 0:
                work = np.dot(a1[lower[k]:upper[k]+1, :].T,
                              a2[lower[k]:upper[k]+1, :])
                wb = np.dot(ydata[lower[k]:upper[k]+1],
                            a2[lower[k]:upper[k]+1, :])
                alpha.T.flat[bo+itop*bw] += work.flat[bi]
                beta[itop:ibottom+1] += wb
        min_influence = 1.0e-10 * invvar.sum() / nfull
        errb = cholesky_band(alpha, mininf=min_influence)
        if isinstance(errb[0], int) and errb[0] == -1:
            a = errb[1]
        else:
            yfit, foo = self.value(xdata, x2=x2, action=a1, upper=upper,
                                   lower=lower)
            return self.maskpoints(errb[0]), yfit
        sol = cholesky_solve(a, beta)
        if self.coeff.ndim == 2:
            # JFH made major bug fix here.
            self.icoeff[:, goodbk] = np.array(
                a[0, 0:nfull].T.reshape(self.npoly, nn, order='F'),
                dtype=a.dtype)
            self.coeff[:, goodbk] = np.array(
                sol[0:nfull].T.reshape(self.npoly, nn, order='F'),
                dtype=sol.dtype)
        else:
            self.icoeff[goodbk] = np.array(a[0, 0:nfull], dtype=a.dtype)
            self.coeff[goodbk] = np.array(sol[0:nfull], dtype=sol.dtype)
        yfit, foo = self.value(xdata, x2=x2, action=a1, upper=upper,
                               lower=lower)
        return 0, yfit

    def action(self, x, x2=None):
        """Construct banded B-spline matrix, with dimensions [ndata, bandwidth].

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Independent variable.
        x2 : :class:`numpy.ndarray`, optional
            Orthogonal dependent variable for 2d fits.

        Returns
        -------
        :func:`tuple`
            A tuple containing the B-spline action matrix; the 'lower'
            parameter, a list of pixel positions, each corresponding to the
            first occurence of position greater than breakpoint indx; and
            'upper', Same as lower, but denotes the upper pixel positions.
        """
        nx = x.size
        nbkpt = self.mask.sum()
        if nbkpt < 2*self.nord:
            return -2, 0, 0
        n = nbkpt - self.nord
        bw = self.npoly*self.nord
        lower = np.zeros((n - self.nord + 1,), dtype=int)
        upper = np.zeros((n - self.nord + 1,), dtype=int) - 1
        indx = self.intrv(x)
        bf1 = self.bsplvn(x, indx)
        action = bf1
        aa = uniq(indx, np.arange(indx.size, dtype=int))
        upper[indx[aa]-self.nord+1] = aa
        rindx = indx[::-1]
        bb = uniq(rindx, np.arange(rindx.size, dtype=int))
        lower[rindx[bb]-self.nord+1] = nx - bb - 1
        if x2 is not None:
            if x2.size != nx:
                raise ValueError('Dimensions of x and x2 do not match.')
            x2norm = 2.0 * (x2 - self.xmin) / (self.xmax - self.xmin) - 1.0
            if self.funcname == 'poly':
                temppoly = np.ones((nx, self.npoly), dtype=float)
                for i in range(1, self.npoly):
                    temppoly[:, i] = temppoly[:, i-1] * x2norm
            elif self.funcname == 'poly1':
                temppoly = np.tile(x2norm, self.npoly).reshape(nx, self.npoly)
                for i in range(1, self.npoly):
                    temppoly[:, i] = temppoly[:, i-1] * x2norm
            elif self.funcname == 'chebyshev':
                temppoly = fchebyshev(x2norm, self.npoly).T
            elif self.funcname == 'legendre':
                temppoly = flegendre(x2norm, self.npoly).T
            else:
                raise ValueError('Unknown value of funcname.')
            action = np.zeros((nx, bw), dtype=float, order='F')
            counter = -1
            for ii in range(self.nord):
                for jj in range(self.npoly):
                    counter += 1
                    action[:, counter] = bf1[:, ii]*temppoly[:, jj]
        return action, lower, upper

    def intrv(self, x):
        """Find the segment between breakpoints which contain each value
        in the array `x`.

        The minimum breakpoint is ``nbkptord - 1``, and the maximum
        is ``nbkpt - nbkptord - 1``.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Data values, assumed to be monotonically increasing.

        Returns
        -------
        :class:`numpy.ndarray`
            Position of array elements with respect to breakpoints.
        """
        gb = self.breakpoints[self.mask]
        n = gb.size - self.nord
        indx = np.zeros((x.size,), dtype='i4')
        ileft = self.nord - 1
        for i in range(x.size):
            while x[i] > gb[ileft+1] and ileft < n - 1:
                ileft += 1
            indx[i] = ileft
        return indx

    def bsplvn(self, x, ileft):
        """Calculates the value of all possibly nonzero B-splines at `x`
        of a certain order.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Independent variable.
        ileft : :class:`int`
            Breakpoint segements that contain `x`.

        Returns
        -------
        :class:`numpy.ndarray`
            B-spline values.
        """
        bkpt = self.breakpoints[self.mask]
        vnikx = np.zeros((x.size, self.nord), dtype=x.dtype, order='F')
        deltap = vnikx.copy()
        deltam = vnikx.copy()
        j = 0
        vnikx[:, 0] = 1.0
        while j < self.nord - 1:
            ipj = ileft+j+1
            deltap[:, j] = bkpt[ipj] - x
            imj = ileft-j
            deltam[:, j] = x - bkpt[imj]
            vmprev = 0.0
            for k in range(j+1):
                vm = vnikx[:, k]/(deltap[:, k] + deltam[:, j-k])
                vnikx[:, k] = vm*deltap[:, k] + vmprev
                vmprev = vm*deltam[:, j-k]
            j += 1
            vnikx[:, j] = vmprev
        return vnikx

    def value(self, x, x2=None, action=None, lower=None, upper=None):
        """Evaluate a B-spline at specified values.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Independent variable.
        x2 : :class:`numpy.ndarray`, optional
            Orthogonal dependent variable for 2d fits.
        action : :class:`numpy.ndarray`, optional
            Action matrix to use.  If not supplied it is calculated.
        lower : :class:`numpy.ndarray`, optional
            If the action parameter is supplied, this parameter must also
            be supplied.
        upper : :class:`numpy.ndarray`, optional
            If the action parameter is supplied, this parameter must also
            be supplied.

        Returns
        -------
        :func:`tuple`
            A tuple containing the results of the bspline evaluation and a
            mask indicating where the evaluation was good.
        """
        xsort = x.argsort()
        xwork = x[xsort]
        if x2 is not None:
            x2work = x2[xsort]
        else:
            x2work = None
        if action is not None:
            if lower is None or upper is None:
                raise ValueError('Must specify lower and upper '
                                 'if action is set.')
        else:
            action, lower, upper = self.action(xwork, x2=x2work)
        yfit = np.zeros(x.shape, dtype=x.dtype)
        bw = self.npoly * self.nord
        spot = np.arange(bw, dtype='i4')
        goodbk = self.mask.nonzero()[0]
        coeffbk = self.mask[self.nord:].nonzero()[0]
        n = self.mask.sum() - self.nord
        if self.npoly > 1:
            goodcoeff = self.coeff[:, coeffbk]
        else:
            goodcoeff = self.coeff[coeffbk]
        # maskthis = np.zeros(xwork.shape,dtype=xwork.dtype)
        for i in range(n-self.nord+1):
            ict = upper[i] - lower[i] + 1
            if ict > 0:
                yfit[lower[i]:upper[i] + 1] = np.dot(
                    action[lower[i]:upper[i] + 1, :],
                    (goodcoeff.flatten('F'))[i * self.npoly + spot])
        yy = yfit.copy()
        yy[xsort] = yfit
        mask = np.ones(x.shape, dtype='bool')
        gb = self.breakpoints[goodbk]
        outside = ((x < gb[self.nord-1]) | (x > gb[n]))
        if outside.any():
            mask[outside] = False
        hmm = ((np.diff(goodbk) > 2).nonzero())[0]
        for jj in range(hmm.size):
            inside = ((x >= self.breakpoints[goodbk[hmm[jj]]]) &
                      (x <= self.breakpoints[goodbk[hmm[jj]+1]-1]))
            if inside.any():
                mask[inside] = False
        return yy, mask

    def maskpoints(self, err):
        """Perform simple logic of which breakpoints to mask.

        Parameters
        ----------
        err : :class:`numpy.ndarray`
            The list of indexes returned by the cholesky routines.

        Returns
        -------
        :class:`int`
            An integer indicating the results of the masking.  -1 indicates
            that the error points were successfully masked.  -2 indicates
            failure; the calculation should be aborted.

        Notes
        -----
        The mask attribute is modified, assuming it is possible to create the
        mask.
        """
        nbkpt = self.mask.sum()
        goodbkpt = np.where(self.mask)[0]
        nbkpt = len(goodbkpt)
        if nbkpt <= 2*self.nord:
            return -2
        hmm = err[uniq(err//self.npoly)]//self.npoly
        n = nbkpt - self.nord
        if np.any(hmm >= n):
            return -2
        test = np.zeros(nbkpt, dtype='bool')
        for jj in range(-int(np.ceil(self.nord/2)), int(self.nord/2.)):
            foo = np.where((hmm+jj) > 0, hmm+jj, np.zeros(hmm.shape,
                                                          dtype=hmm.dtype))
            inside = np.where((foo+self.nord) < n-1, foo+self.nord,
                              np.zeros(hmm.shape, dtype=hmm.dtype)+n-1)
            if len(inside) > 0:
                test[inside] = True
        if test.any():
            reality = goodbkpt[test]
            if self.mask[reality].any():
                self.mask[reality] = False
                return -1
            else:
                return -2
        else:
            return -2

    def workit(self, xdata, ydata, invvar, action, lower, upper):
        """An internal routine for bspline_extract and bspline_radial which
        solve a general banded correlation matrix which is represented by the
        variable "action".  This routine only solves the linear system once,
        and stores the coefficients in sset. A non-zero return value signifies
        a failed inversion
        Parameters
        ----------
        xdata : :class:`numpy.ndarray`
            Independent variable.
        ydata : :class:`numpy.ndarray`
            Dependent variable.
        invvar : :class:`numpy.ndarray`
            Inverse variance of `ydata`.
        action : :class:`numpy.ndarray`
            Banded correlation matrix
        lower  : :class:`numpy.ndarray`
            A list of pixel positions, each corresponding to the first
            occurence of position greater than breakpoint indx
        upper  : :class:`numpy.ndarray`
            Same as lower, but denotes the upper pixel positions
        Returns
        -------
        :func:`tuple`
            A tuple containing an integer error code, and the evaluation of the
            b-spline at the input values.  The error codes are as follows:
             -2:  is a failure,
             -1:  indicates dropped breakpoints
              0:  is success
              positive integers:  indicate ill-conditioned breakpoints.
        """
        goodbk = self.mask[self.nord:]
        nn = goodbk.sum()
        if nn < self.nord:
            yfit = np.zeros(ydata.shape, dtype=float)
            return -2, yfit
        nfull = nn * self.npoly
        bw = self.npoly * self.nord
        foo = np.sqrt(np.tile(invvar, bw).reshape(bw, invvar.size).transpose())
        a2 = action * foo
        # a2 = action*np.sqrt(np.outer(invvar,np.ones(bw)))

        alpha = np.zeros((bw, nfull + bw), dtype='d')
        beta = np.zeros((nfull + bw,), dtype='d')
        bi = np.arange(bw, dtype='i4')
        bo = np.arange(bw, dtype='i4')
        for k in range(1, bw):
            bi = np.append(bi, np.arange(bw - k, dtype='i4') + (bw + 1) * k)
            bo = np.append(bo, np.arange(bw - k, dtype='i4') + bw * k)
        for k in range(nn - self.nord + 1):
            itop = k * self.npoly
            ibottom = min(itop, nfull) + bw - 1
            ict = upper[k] - lower[k] + 1
            if ict > 0:
                work = np.dot(a2[lower[k]:upper[k] + 1, :].T,
                              a2[lower[k]:upper[k] + 1, :])
                wb = np.dot(ydata[lower[k]:upper[k] + 1] * np.sqrt(
                    invvar[lower[k]:upper[k] + 1]),
                            a2[lower[k]:upper[k] + 1, :])
                alpha.T.flat[bo + itop * bw] += work.flat[bi]
                beta[itop:ibottom + 1] += wb
        min_influence = 1.0e-10 * invvar.sum() / nfull
        # Right now we are not returning the covariance,
        # although it may arise that we should
        # covariance = alpha
        errb = cholesky_band(alpha, mininf=min_influence)  # ,verbose=True)
        if isinstance(errb[0], int) and errb[0] == -1:
            a = errb[1]
        else:
            yfit, foo = self.value(xdata, x2=xdata, action=action, upper=upper,
                                   lower=lower)
            return self.maskpoints(errb[0]), yfit
        errs = cholesky_solve(a, beta)
        if isinstance(errs[0], int) and errs[0] == -1:
            sol = errs[1]
        else:
            #
            # It is not possible for this to get called, because cholesky_solve
            # has only one return statement, & that statement guarantees that
            # errs[0] == -1
            #
            yfit, foo = self.value(xdata, x2=xdata, action=action, upper=upper,
                                   lower=lower)
            return self.maskpoints(errs[0]), yfit
        if self.coeff.ndim == 2:
            self.icoeff[:, goodbk] = np.array(
                a[0, 0:nfull].T.reshape(self.npoly, nn, order='F'),
                dtype=a.dtype)
            self.coeff[:, goodbk] = np.array(
                sol[0:nfull].T.reshape(self.npoly, nn, order='F'),
                dtype=sol.dtype)
        else:
            self.icoeff[goodbk] = np.array(a[0, 0:nfull], dtype=a.dtype)
            self.coeff[goodbk] = np.array(sol[0:nfull], dtype=sol.dtype)
        yfit, foo = self.value(xdata, x2=xdata, action=action, upper=upper,
                               lower=lower)
        return 0, yfit


def cholesky_band(ndl, mininf=0.0):
    """Compute *lower* Cholesky decomposition of a banded matrix.

    This function provides informative error messages to pass back to the
    :class:`~pydl.pydlutils.bspline.Bspline` machinery; the actual
    computation is delegated to :func:`scipy.linalg.cholesky_banded`.

    Parameters
    ----------
    ndl : :class:`numpy.ndarray`
        A matrix on which to perform the Cholesky decomposition.  The
        matrix must be in a special, *lower* form described in
        :func:`scipy.linalg.cholesky_banded`.  In addition, the input
        must be padded.  If the original, square matrix has size
        :math:`N \\times N`, and the width of the band is :math:`b`,
        `l` must be :math:`b \\times (N + b)`.
    mininf : :class:`float`, optional
        Entries in the `l` matrix are considered negative if they are less
        than this value (default 0.0).

    Returns
    -------
    :func:`tuple`
        If problems were detected, the first item will be the index or
        indexes where the problem was detected, and the second item will simply
        be the input matrix.  If no problems were detected, the first item
        will be -1, and the second item will be the Cholesky decomposition.
    """
    bw, nn = ndl.shape
    n = nn - bw
    negative = ndl[0, 0:n] <= mininf
    if negative.any() or not np.all(np.isfinite(ndl)):
        warn('Bad entries: ' + str(negative.nonzero()[0]), PydlutilsUserWarning)
        return negative.nonzero()[0], ndl
    try:
        lower = cholesky_banded(ndl[:, 0:n], lower=True)
    except LinAlgError:
        #
        # Figure out where the error is.
        #
        lower = ndl.copy()
        kn = bw - 1
        spot = np.arange(kn, dtype='i4') + 1
        for j in range(n):
            lower[0, j] = np.sqrt(lower[0, j])
            lower[spot, j] /= lower[0, j]
            x = lower[spot, j]
            if not np.all(np.isfinite(x)):
                warn('NaN found in cholesky_band.', PydlutilsUserWarning)
                return j, ndl
    #
    # Restore padding.
    #
    ll = np.zeros(ndl.shape, dtype=ndl.dtype)
    ll[:, 0:n] = lower
    return -1, ll


def cholesky_solve(a, bb):
    """Solve the equation :math:`A x = b` where `a` is a *lower*
    Cholesky-banded matrix.

    In the :class:`~pydl.pydlutils.bspline.Bspline` machinery, `a` needs to
    be padded.  This function should only used with the output of
    :func:`~pydl.pydlutils.bspline.cholesky_band`, to ensure the proper
    padding on `a`.  Otherwise the computation is delegated to
    :func:`scipy.linalg.cho_solve_banded`.

    Parameters
    ----------
    a : :class:`numpy.ndarray`
        *Lower* Cholesky decomposition of :math:`A` in :math:`A x = b`.
    bb : :class:`numpy.ndarray`
        :math:`b` in :math:`A x = b`.

    Returns
    -------
    :class:`numpy.ndarray`
        The solution, padded to be the same shape as `bb`.
    """
    bw = a.shape[0]
    n = bb.shape[0] - bw
    x = np.zeros(bb.shape, dtype=bb.dtype)
    x[0:n] = cho_solve_banded((a[:, 0:n], True), bb[0:n])
    return x


def iterfit(xdata, ydata, invvar=None, upper=5, lower=5, x2=None,
            maxiter=10, nord=4, bkpt=None, fullbkpt=None,
            kwargs_bspline={}, kwargs_reject={}):
    """Iteratively fit a B-spline set to data, with rejection.

    Parameters
    ----------
    xdata : :class:`numpy.ndarray`
        Independent variable.
    ydata : :class:`numpy.ndarray`
        Dependent variable.
    invvar : :class:`numpy.ndarray`, optional
        Inverse variance of `ydata`.  If not set, it will be calculated based
        on the standard deviation.
    upper : :class:`int` or :class:`float`, optional
        Upper rejection threshold in units of sigma, defaults to 5 sigma.
    lower : :class:`int` or :class:`float`, optional
        Lower rejection threshold in units of sigma, defaults to 5 sigma.
    x2 : :class:`numpy.ndarray`, optional
        Orthogonal dependent variable for 2d fits.
    maxiter : :class:`int`, optional
        Maximum number of rejection iterations, default 10.  Set this to
        zero to disable rejection.
    nord : :class:`int`
    bkpt : :None:
    fullbkpt : :class:`numpy.ndarray`
    kwargs_bspline : :class:`list`
    kwargs_reject : :class:`list`

    Returns
    -------
    :func:`tuple`
        A tuple containing the fitted bspline object and an output mask.
    """
    nx = xdata.size
    if ydata.size != nx:
        raise ValueError('Dimensions of xdata and ydata do not agree.')
    if invvar is not None:
        if invvar.size != nx:
            raise ValueError('Dimensions of xdata and invvar do not agree.')
    else:
        #
        # This correction to the variance makes it the same
        # as IDL's variance()
        #
        var = ydata.var()*(float(nx)/float(nx-1))
        if var == 0:
            var = 1.0
        invvar = np.ones(ydata.shape, dtype=ydata.dtype)/var
    if x2 is not None:
        if x2.size != nx:
            raise ValueError('Dimensions of xdata and x2 do not agree.')
    yfit = np.zeros(ydata.shape, dtype=ydata.dtype)
    if invvar.size == 1:
        outmask = True
    else:
        outmask = np.ones(invvar.shape, dtype='bool')
    xsort = xdata.argsort()
    maskwork = (outmask & (invvar > 0))[xsort]
    if 'oldset' in kwargs_bspline:
        sset = kwargs_bspline['oldset']
        sset.mask = True
        sset.coeff = 0
    else:
        if not maskwork.any():
            raise ValueError('No valid data points.')
#        if 'fullbkpt' in kwargs:
#            fullbkpt = kwargs['fullbkpt']
        else:
            sset = Bspline(xdata[xsort[maskwork]], nord=nord, bkpt=bkpt,
                           fullbkpt=fullbkpt, **kwargs_bspline)
            if maskwork.sum() < sset.nord:
                warn('Number of good data points fewer than nord.',
                     PydlutilsUserWarning)
                return sset, outmask
            if x2 is not None:
                if 'xmin' in kwargs_bspline:
                    xmin = kwargs_bspline['xmin']
                else:
                    xmin = x2.min()
                if 'xmax' in kwargs_bspline:
                    xmax = kwargs_bspline['xmax']
                else:
                    xmax = x2.max()
                if xmin == xmax:
                    xmax = xmin + 1
                sset.xmin = xmin
                sset.xmax = xmax
                if 'funcname' in kwargs_bspline:
                    sset.funcname = kwargs_bspline['funcname']
    xwork = xdata[xsort]
    ywork = ydata[xsort]
    invwork = invvar[xsort]
    if x2 is not None:
        x2work = x2[xsort]
    else:
        x2work = None
    iiter = 0
    error = 0
    qdone = False
    while (error != 0 or qdone is False) and iiter <= maxiter:
        goodbk = sset.mask.nonzero()[0]
        if maskwork.sum() <= 1 or not sset.mask.any():
            sset.coeff = 0
            iiter = maxiter + 1  # End iterations
        else:
            if 'requiren' in kwargs_bspline:
                i = 0
                while xwork[i] < sset.breakpoints[goodbk[sset.nord]] and i < nx-1:
                    i += 1
                ct = 0
                for ileft in range(sset.nord, sset.mask.sum()-sset.nord+1):
                    while (sset.breakpoints[goodbk[ileft]] <= xwork[i]
                           < sset.breakpoints[goodbk[ileft + 1]] and
                           i < nx - 1):
                        ct += invwork[i]*maskwork[i] > 0
                        i += 1
                    if ct >= kwargs_bspline['requiren']:
                        ct = 0
                    else:
                        sset.mask[goodbk[ileft]] = False
            error, yfit = sset.fit(xwork, ywork, invwork*maskwork,
                                   x2=x2work)
        iiter += 1
        inmask = maskwork
        if error == -2:

            return sset, outmask
        elif error == 0:
            maskwork, qdone = djs_reject(ywork, yfit, invvar=invwork,
                                         inmask=inmask, outmask=maskwork,
                                         upper=upper, lower=lower,
                                         **kwargs_reject)
        else:
            pass
    outmask[xsort] = maskwork
    temp = yfit
    yfit[xsort] = temp
    return sset, outmask
