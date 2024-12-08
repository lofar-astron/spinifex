from __future__ import annotations

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from spinifex.geometry import get_ipp


def test_geometry():
    """Test the geometry module"""
    source = SkyCoord.from_name("Cas A")
    lon = 6.367 * u.deg
    lat = 52.833 * u.deg
    heights = np.arange(100, 2000, 100) * u.km
    dwingeloo = EarthLocation(lon=lon, lat=lat, height=0 * u.km)
    times = Time("2020-03-21T01:00:00") + np.arange(0, 10) * 15 * u.min
    ipp = get_ipp.get_ipp_from_skycoord(
        loc=dwingeloo, times=times, source=source, height_array=heights
    )
    assert isinstance(ipp, get_ipp.IPP)
    assert isinstance(ipp.airmass, np.ndarray)
