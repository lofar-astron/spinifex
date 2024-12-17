"""Time utilities for ionospheric analysis."""

from __future__ import annotations

import astropy.units as u
import numpy as np
from astropy.time import Time
from numpy.typing import ArrayLike


def get_gps_week(time: Time) -> ArrayLike:
    """Get the GPS week from a time.

    Parameters
    ----------
    time : Time
        Time(s) to get the GPS week from

    Returns
    -------
    ArrayLike
        GPS week(s)
    """
    return np.floor((time.gps * u.s).to(u.week).value).astype(int)


def get_unique_days(times: Time) -> Time:
    """Get the unique days from a list of times.

    Parameters
    ----------
    times : Time
        Times to get the unique days from.

    Returns
    -------
    Time
        Unique days
    """
    return Time(np.sort(np.unique(np.floor(times.mjd))), format="mjd")


def get_unique_days_index(times: Time) -> ArrayLike:
    """Get the unique days index from a list of times.

    Parameters
    ----------
    times : Time
        Times to get the unique days index from.

    Returns
    -------
    ArrayLike
        Unique days index
    """
    return np.searchsorted(get_unique_days(times).mjd, np.floor(times.mjd))
