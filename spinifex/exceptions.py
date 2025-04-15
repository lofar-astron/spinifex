from __future__ import annotations


class SpinifexError(Exception):
    """Base class for errors in spinifex"""


class IonexError(SpinifexError):
    """Error in IONEX files."""


class TomionError(SpinifexError):
    """Error in TomIon files."""


class TimeResolutionError(IonexError):
    """Error in IONEX resolution."""


class FITSHeaderError(SpinifexError):
    """Error in FITS header"""


class ArrayShapeError(SpinifexError):
    """Error in array shapes"""
