#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module for getting the Ionospheric Piercepoints"""

from __future__ import annotations

from typing import NamedTuple

import astropy.units as u
import numpy as np
from astropy.constants import R_earth
from astropy.coordinates import ITRS, AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from numpy.typing import ArrayLike


class IPP(NamedTuple):
    """Ionospheric Piercepoints"""

    loc: EarthLocation
    """location of the piercepoints, dimension: times x altitudes. All altitudes are assumed to be equal"""
    times: Time
    """array of times"""
    los: SkyCoord
    """Line of sight direction in ITRS coordinates"""
    airmass: ArrayLike
    """airmass factor to convert to slant values"""
    altaz: AltAz
    """azimuth elevation"""


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
    # Note: at the moment we calculate ipp per station. I think this is ok,
    # but maybe we need to include a many stations option
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
        ionospheric piercepoints
    """
    los_dir = altaz.transform_to(ITRS)
    ipp, airmass = _get_ipp_simple(height_array=height_array, loc=loc, los_dir=los_dir)
    return IPP(
        loc=EarthLocation.from_geocentric(*ipp),
        times=altaz.obstime,
        los=los_dir,
        airmass=airmass,
        altaz=altaz,
    )


def _get_ipp_simple(
    height_array: u.Quantity, loc: EarthLocation, los_dir: SkyCoord
) -> tuple:
    """helper function to calculate ionospheric piercepoints using a simple spherical earth model
    |loc + alphas * los_dir| = R_earth + height_array, solve for alphas using abc formula
    Parameters
    ----------
    height_array : u.Quantity
        array of altitudes
    loc : EarthLocation
        observer location
    los_dir : ITRS
        line of sight, unit vector

    Returns
    -------
    tuple(list[u.Quantity], ArrayLike)
        ipp.x, ipp.y, ipp.z positions, airmass
    """
    c_value = R_earth**2 - (R_earth + height_array) ** 2
    b_value = u.Quantity(loc.geocentric) @ los_dir.cartesian.xyz
    alphas = -b_value + np.sqrt(b_value**2 - c_value[:, np.newaxis])
    ipp = (
        u.Quantity(loc.geocentric)[:, np.newaxis, np.newaxis]
        + alphas[np.newaxis] * los_dir.cartesian.xyz[:, np.newaxis]
    )
    inv_airmass = np.einsum("ijk,ik->jk", ipp, los_dir.cartesian.xyz)
    inv_airmass /= R_earth + height_array[:, np.newaxis]  # normalized
    airmass = 1.0 / inv_airmass.value
    return ipp, airmass
