#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module for getting the Ionospheric Piercepoints"""

from __future__ import annotations

from typing import NamedTuple

import astropy.units as u
import numpy as np
from astropy.coordinates import ITRS, AltAz, EarthLocation, SkyCoord
from astropy.time import Time

r_earth = 6364.62 * u.km


class IPP(NamedTuple):
    loc: EarthLocation
    """location of the piercepoints, dimension: times x altitudes. All altitudes are assumed to be equal"""
    times: Time
    """array of times"""
    los: SkyCoord
    """Line of sight direction in ITRS coordinates"""


def get_ipp_from_skycoord(
    loc: EarthLocation, times: Time, source: SkyCoord, height_array: u.Quantity
) -> IPP:
    """Get ionospheric piercepoints

    Parameters
    ----------
    loc : EarthLocation
        observer location
    times : Time
        observation times
    source : SkyCoord
        source location
    height_array : u.Quantity
        array of altitudes

    Returns
    -------
    IPP
        Ionospheric piercepoints
    """
    aa = AltAz(location=loc, obstime=times)
    altaz = source.transform_to(aa)
    return get_ipp_from_altaz(loc, altaz, height_array)


def get_ipp_from_altaz(
    loc: EarthLocation, altaz: AltAz, height_array: u.Quantity
) -> IPP:
    """get ionospheric piercepoints from azimuth elevations

    Parameters
    ----------
    loc : EarthLocation
        observer location
    altaz : AltAz
        azimuth and elevations for all times
    height_array : u.Quantity
        array of altitudes

    Returns
    -------
    IPP
        _description_
    """
    los_dir = altaz.transform_to(ITRS)
    ipp = _get_ipp_simple(height_array=height_array, loc=loc, los_dir=los_dir)
    return IPP(
        loc=EarthLocation.from_geocentric(*ipp), times=altaz.obstime, los=los_dir
    )


def _get_ipp_simple(
    height_array: u.Quantity, loc: EarthLocation, los_dir: SkyCoord
) -> list[u.Quantity]:
    """helper function to calculate ionospheric piercepoints using a simple spherical earth model

    Parameters
    ----------
    height_array : u.Quantity
        array of altitudes
    loc : EarthLocation
        observer location
    los_dir : SkyCoord
        line of sight (ITRS)

    Returns
    -------
    list[u.Quantity]
        ipp.x, ipp.y, ipp.z positions
    """
    c_value = r_earth**2 - (r_earth + height_array) ** 2
    ipp_vector = los_dir.cartesian[np.newaxis] * height_array[:, np.newaxis]
    b_value = np.sqrt(
        loc.x * ipp_vector.x + loc.y * ipp_vector.y + loc.z * ipp_vector.z
    )  # inproduct
    alphas = -b_value + np.sqrt(b_value**2 - c_value[:, np.newaxis])
    return [
        loc.x + alphas * los_dir.x,
        loc.y + alphas * los_dir.y,
        loc.z + alphas * los_dir.z,
    ]
