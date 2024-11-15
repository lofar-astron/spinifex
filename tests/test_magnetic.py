#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Testing of magnetic"""
from unittest import TestCase

import astropy.units as u
from astropy.coordinates import EarthLocation

from spinifex.magnetic import get_magnetic_field


def is_convertible_to_unit(quantity: u.Quantity, unit: u.Unit) -> bool:
    """Test if unit is convertible to a given unit"""
    try:
        _ = quantity.to(unit)
        return True
    except u.UnitsError:
        return False


class TestMagnetic(TestCase):
    """Testing of magnetic"""

    def test_get_magnetic_field(self):
        """Test that get_magnetic does not crash"""
        above_dwingeloo = EarthLocation(
            lon=6.367 * u.deg, lat=52.833 * u.deg, height=100 * u.km
        )
        field = get_magnetic_field(above_dwingeloo)
        assert is_convertible_to_unit(field, u.tesla)
