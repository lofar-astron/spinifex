from __future__ import annotations


class IonexError(Exception):
    """Error in IONEX files."""


class TimeResolutionError(IonexError):
    """Error in IONEX resolution."""
