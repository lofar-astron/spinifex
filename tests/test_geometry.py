from __future__ import annotations

from astropy.utils import iers

iers.conf.auto_download = False

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from spinifex.geometry.get_ipp import IPP, get_ipp_from_skycoord


def get_test_ipp() -> IPP:
    cas_a = SkyCoord(ra=350.85 * u.deg, dec=58.815 * u.deg)
    lon = 6.367 * u.deg
    lat = 52.833 * u.deg
    heights = np.arange(100, 2000, 100) * u.km
    dwingeloo = EarthLocation(lon=lon, lat=lat, height=0 * u.km)
    times = Time("2020-01-20T01:00:00") + np.arange(0, 10) * 15 * u.min
    return get_ipp_from_skycoord(
        loc=dwingeloo, times=times, source=cas_a, height_array=heights
    )


def test_geometry():
    """Test the geometry module"""
    ipp = get_test_ipp()
    assert isinstance(ipp, IPP)
    assert isinstance(ipp.airmass, np.ndarray)
