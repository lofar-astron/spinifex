#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0
# pylint: disable=duplicate-code
"""Testing of the main module"""

from __future__ import annotations

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from spinifex import get_rm


def test_get_rm():
    """Test that get_rm does not crash"""
    times = Time("2020-01-08T01:00:00") + np.arange(10) * 25 * u.min
    source = SkyCoord.from_name("Cas A")
    lon = 6.367 * u.deg
    lat = 52.833 * u.deg
    dwingeloo = EarthLocation(lon=lon, lat=lat, height=0 * u.km)
    rm = get_rm.get_rm_from_skycoord(loc=dwingeloo, times=times, source=source)
    assert isinstance(rm.rm, np.ndarray)
    assert rm.rm.shape == times.shape
    assert np.isclose(rm.rm[0], 68.7, 0.1)
