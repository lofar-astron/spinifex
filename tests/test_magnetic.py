#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Testing of magnetic"""

from __future__ import annotations

import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
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
    above_dwingeloo = EarthLocation(
        lon=6.367 * u.deg, lat=52.833 * u.deg, height=100 * u.km
    )
    time = Time("2024-01-01T13:42")
    field = models.ppigrf(above_dwingeloo, time)
    assert is_convertible_to_unit(field, u.tesla)
