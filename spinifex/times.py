"""Time utilities for ionospheric analysis."""

from __future__ import annotations

import astropy.units as u
import numpy as np
from astropy.time import Time
from numpy.typing import NDArray


def get_gps_week(time: Time) -> NDArray[np.int32]:
    """Get the GPS week from a time.

    Parameters
    ----------
    time : Time
        Time(s) to get the GPS week from

    Returns
    -------
    NDArray
        GPS week(s)
    """
    return np.array(np.floor((time.gps * u.s).to(u.week).value), dtype=np.int32)


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


def get_unique_days_index(times: Time) -> NDArray[np.int32]:
    """Get the unique days index from a list of times.

    Parameters
    ----------
    times : Time
        Times to get the unique days index from.

    Returns
    -------
    NDArray
        Unique days index
    """
    return np.array(
        np.searchsorted(get_unique_days(times).mjd, np.floor(times.mjd)), dtype=np.int32
    )


# TODO: Add docs
def get_indexlist_unique_days(unique_days: Time, times: Time) -> NDArray[np.bool_]:
    return np.array(
        np.floor(times.mjd)[np.newaxis] == unique_days.mjd[:, np.newaxis],
        dtype=np.bool_,
    )


# TODO: Add docs
def get_consecutive_days(unique_days: Time) -> NDArray[np.int32]:
    return np.cumsum(np.diff(unique_days.mjd, prepend=unique_days.mjd[0]) > 1)
