"""
Asymmetric Polynomial Transform Code

Shamelessly pilfered from scikit image _geometric, which contains all the
transforms in the standard repertoir, but which did not include a 2D polynomial
transform that could have different orders for X and Y.

The asymmetric polynomial algorithm was ported to Python from the IDL routine,
mod_polywarp.pro, written by M. Matuszewski.  This routine was part of the
original IDL KCWI pipeline, which can be found at:

    * https://github.com/Keck-DataReductionPipelines/KcwiDRP

The resulting transform is compatible with the scikit image transform.warp
function.

"""
import math
import numpy as np
from scipy import linalg


class GeometricTransform(object):
    """Base class for geometric transformations.

    """
    def __call__(self, coords):
        """Apply forward transformation.

        Parameters
        ----------
        coords : (N, 2) array
            Source coordinates.

        Returns
        -------
        coords : (N, 2) array
            Destination coordinates.

        """
        raise NotImplementedError()

    def inverse(self, coords):
        """Apply inverse transformation.

        Parameters
        ----------
        coords : (N, 2) array
            Destination coordinates.

        Returns
        -------
        coords : (N, 2) array
            Source coordinates.

        """
        raise NotImplementedError()

    def residuals(self, src, dst):
        """Determine residuals of transformed destination coordinates.

        For each transformed source coordinate the euclidean distance to the
        respective destination coordinate is determined.

        Parameters
        ----------
        src : (N, 2) array
            Source coordinates.
        dst : (N, 2) array
            Destination coordinates.

        Returns
        -------
        residuals : (N, ) array
            Residual for coordinate.

        """
        return np.sqrt(np.sum((self(src) - dst)**2, axis=1))

    def __add__(self, other):
        """Combine this transformation with another.

        """
        raise NotImplementedError()


class AsymmetricPolynomialTransform(GeometricTransform):
    """2D polynomial transformation.

    Has the following form::

        X = sum[j=0:order]( sum[i=0:j]( a_ji * x**(j - i) * y**i ))
        Y = sum[j=0:order]( sum[i=0:j]( b_ji * x**(j - i) * y**i ))

    Parameters
    ----------
    params : (2, N) array, optional
        Polynomial coefficients where `N * 2 = (order + 1) * (order + 2)`. So,
        a_ji is defined in `params[0, :]` and b_ji in `params[1, :]`.

    Attributes
    ----------
    params : (2, N) array
        Polynomial coefficients where `N * 2 = (order + 1) * (order + 2)`. So,
        a_ji is defined in `params[0, :]` and b_ji in `params[1, :]`.

    """

    def __init__(self, params=None):
        if params is None:
            # default to transformation which preserves original coordinates
            params = np.array([[0, 1, 0], [0, 0, 1]])
        if params.shape[0] != 2:
            raise ValueError("invalid shape of transformation parameters")
        self.params = params

    def estimate(self, src, dst, order=(2, 2)):
        """Generate a pair of warping polynomial kernels in two
           dimensions (Kx, Ky) taking  starting coordinates (Xo,Yo) to
           ending coordinates (Xi,Yi). The polynomial degrees in the
           two dimensions are allowed to differ:

           Xi = Sum (i,0,Degree[1]) (j,0,Degree[0]) Kx[i,j] Xo^j Yo^i
           Yi = Sum (i,0,Degree[1]) (j,0,Degree[0]) Ky[i,j] Xo^j Yo^i

           The degrees need not be equal, but the number of points provided
           must be greater than or equal to

           number of points = (Degree[0]+1)*(Degree[1]+1)

        Parameters
        ----------
        src : (N, 2) array
            Source coordinates.
        dst : (N, 2) array
            Destination coordinates.
        order : (2) array, int
            Polynomial order (number of coefficients is order + 1).

        Returns
        -------
        success : bool
            True, if model estimation succeeds.

        """
        if isinstance(order, int):
            order = (order, order)

        m = order[0] + 1
        n = order[1] + 1
        efford = order[0] + order[1] + 1

        plen = int(efford * (efford + 1) / 2)

        xi = dst[:, 0]
        yi = dst[:, 1]
        xo = src[:, 0]
        yo = src[:, 1]
        npts = src.shape[0]

        W = np.zeros((npts, m*n))
        X = np.zeros((npts, 2))

        medxo = np.median(abs(xo))
        medyo = np.median(abs(yo))

        X[:, 0] = xi/medxo
        X[:, 1] = yi/medyo

        xx = np.zeros(m)
        xx[0] = 1.0
        for ix in range(1, m):
            xx[ix] = xx[ix-1]/medxo

        yy = np.zeros(n)
        yy[0] = 1.0
        for iy in range(1, n):
            yy[iy] = yy[iy-1]/medyo

        zz = np.outer(yy, xx)

        U = np.zeros((npts, 2))
        U[:, 0] = xo / medxo
        U[:, 1] = yo / medyo

        for iy in range(0, n):
            for ix in range(0, m):
                row = ix * n + iy
                W[:, row] = U[:, 1]**iy * U[:, 0]**ix

        # W = W.T

        WW = np.matmul(W.T, W)

        MM = linalg.inv(WW)

        MMM = np.matmul(MM.T, W.T)

        # print('MMM shape', MMM.shape)
        # print('X shape', X.shape)
        # print('zz shape', zz.shape)
        # print('n, m', n, m)

        mmx = np.matmul(MMM, X[:, 0]).reshape((m, n))
        mmy = np.matmul(MMM, X[:, 1]).reshape((m, n))

        # print('mmx shape', mmx.shape)
        # print('mmy shape', mmy.shape)

        kx = zz * mmx.T * medxo
        ky = zz * mmy.T * medyo

        # print("kx", kx)
        # print("ky", ky)

        params = np.zeros((2, plen))

        pidx = 0
        for j in range(efford):
            for ix in range(j + 1):
                iy = j - ix
                if ix <= order[1] and iy <= order[0]:
                    params[0, pidx] = kx[ix, iy]
                    params[1, pidx] = ky[ix, iy]
                else:
                    params[0, pidx] = 0.
                    params[1, pidx] = 0.
                pidx += 1

        # print("params", params)

        self.params = params

        return True

    def __call__(self, coords):
        """Apply forward transformation.

        Parameters
        ----------
        coords : (N, 2) array
            source coordinates

        Returns
        -------
        coords : (N, 2) array
            Transformed coordinates.

        """
        x = coords[:, 0]
        y = coords[:, 1]
        u = len(self.params.ravel())
        # number of coefficients -> u = (order + 1) * (order + 2)
        order = int((- 3 + math.sqrt(9 - 4 * (2 - u))) / 2)
        dst = np.zeros(coords.shape)

        pidx = 0
        for j in range(order + 1):
            for i in range(j + 1):
                dst[:, 0] += self.params[0, pidx] * x ** (j - i) * y ** i
                dst[:, 1] += self.params[1, pidx] * x ** (j - i) * y ** i
                pidx += 1

        return dst

    def inverse(self, coords):
        raise Exception(
            'There is no explicit way to do the inverse polynomial '
            'transformation. Instead, estimate the inverse transformation '
            'parameters by exchanging source and destination coordinates,'
            'then apply the forward transformation.')


TRANSFORMS = {
    'asympolynomial': AsymmetricPolynomialTransform,
}


def estimate_transform(ttype, src, dst, **kwargs):
    """Estimate 2D geometric transformation parameters.

    You can determine the over-, well- and under-determined parameters
    with the total least-squares method.

    Number of source and destination coordinates must match.

    Parameters
    ----------
    ttype : {'asympolynomial'}
        Type of transform.
    src : (N, 2) array
        Source coordinates.
    dst : (N, 2) array
        Destination coordinates.
    kwargs : array or int
        Function parameters (src, dst, n, angle)::

            NAME / TTYPE        FUNCTION PARAMETERS
            'asympolynomial'    `src, `dst`, (`order`, `order`)
                                             (x, y polynomial orders, default
                                             order is 2)

        Also see examples below.

    Returns
    -------
    tform : :class:`GeometricTransform`
        Transform object containing the transformation parameters and providing
        access to forward and inverse transformation functions.

    Examples
    --------
    >>> import numpy as np
    >>> import geometric as tf

    >>> # estimate transformation parameters
    >>> src = np.array([0, 0, 10, 10]).reshape((2, 2))
    >>> dst = np.array([12, 14, 1, -20]).reshape((2, 2))

    >>> tform = tf.estimate_transform('asympolynomial', src, dst)

    >>> np.allclose(tform.inverse(tform(src)), src)
    True

    >>> # warp image using the estimated transformation
    >>> from skimage import data
    >>> image = data.camera()

    >>> warp(image, inverse_map=tform.inverse) # doctest: +SKIP

    >>> # create transformation with explicit parameters
    >>> tform2 = tf.AsymmetriPolynomialTransform(src, dst, order=(2, 3))
    >>> np.allclose(tform2(tform(src)))
    True

    """
    ttype = ttype.lower()
    if ttype not in TRANSFORMS:
        raise ValueError('the transformation type \'%s\' is not'
                         'implemented' % ttype)

    tform = TRANSFORMS[ttype]()
    tform.estimate(src, dst, **kwargs)

    return tform
