#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Testing of magnetic"""

from __future__ import annotations

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from spinifex.geometry import get_ipp
from spinifex.magnetic import models


def is_convertible_to_unit(quantity: u.Quantity, unit: u.Unit) -> bool:
    """Test if unit is convertible to a given unit"""
    try:
        _ = quantity.to(unit)
        return True
    except u.UnitsError:
        return False


def test_get_magnetic_field():
    """Test that get_magnetic does not crash"""
    source = SkyCoord.from_name("Cas A")
    lon = 6.367 * u.deg
    lat = 52.833 * u.deg
    heights = np.arange(100, 2000, 100) * u.km
    dwingeloo = EarthLocation(lon=lon, lat=lat, height=0 * u.km)
    times = Time("2020-03-21T01:00:00") + np.arange(0, 10) * 15 * u.min
    ipp = get_ipp.get_ipp_from_skycoord(
        loc=dwingeloo, times=times, source=source, height_array=heights
    )
    field = models.ppigrf(ipp)
    assert is_convertible_to_unit(field, u.tesla)
