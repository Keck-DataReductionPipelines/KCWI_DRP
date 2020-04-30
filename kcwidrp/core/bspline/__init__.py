import astropy.utils.exceptions as aue


class PydlutilsException(Exception):
    """Base class for exceptions raised in PyDL functions.
    """
    pass


class PydlutilsUserWarning(aue.AstropyUserWarning):
    """Class for warnings issued by :mod:`pydl.pydlutils`.
    """
    pass


__all__ = ['PydlutilsException', 'PydlutilsUserWarning']