#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module for getting the Ionospheric Piercepoints"""

from __future__ import annotations

from typing import NamedTuple

import astropy.units as u
import numpy as np
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


def get_ipp_from_skycoord(
    loc: EarthLocation, times: Time, source: SkyCoord, height_array: ArrayLike
) -> IPP:
    """Get the ionospheric piercepoints as EarhtLocations for a given EarthLocation, time, SkyCoord"""
    # height_array = np.arange(100, 1000, 10) * u.km
    aa = AltAz(location=loc, obstime=times)
    altaz = source.transform_to(aa)
    return get_ipp_from_altaz(loc, altaz, height_array, times)


def get_ipp_from_altaz(
    loc: EarthLocation, altaz: AltAz, height_array: ArrayLike, times: Time
) -> IPP:
    altaz_dir = altaz.transform_to(ITRS)
    ipp_vector = altaz_dir.cartesian[np.newaxis] * height_array[:, np.newaxis]
    ipp = [
        loc.to(u.m).x + ipp_vector.x,
        loc.to(u.m).y + ipp_vector.y,
        loc.to(u.m).z + ipp_vector.z,
    ]
    return IPP(loc=EarthLocation.from_geocentric(*ipp), times=times, los=altaz_dir)
