"""Testing of magnetic"""

from __future__ import annotations

from astropy.utils import iers

iers.conf.auto_download = False

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from spinifex.geometry.get_ipp import get_ipp_from_skycoord
from spinifex.magnetic import magnetic_models


def is_convertible_to_unit(quantity: u.Quantity, unit: u.Unit) -> bool:
    """Test if unit is convertible to a given unit"""
    try:
        _ = quantity.to(unit)
        return True
    except u.UnitsError:
        return False


def test_get_magnetic_field():
    """Test that get_magnetic does not crash"""
    cas_a = SkyCoord(ra=350.85 * u.deg, dec=58.815 * u.deg)
    lon = 6.367 * u.deg
    lat = 52.833 * u.deg
    heights = np.arange(100, 2000, 100) * u.km
    dwingeloo = EarthLocation(lon=lon, lat=lat, height=0 * u.km)
    times = Time("2020-01-20T01:00:00") + np.arange(0, 10) * 15 * u.min
    ipp = get_ipp_from_skycoord(
        loc=dwingeloo, times=times, source=cas_a, height_array=heights
    )
    field = magnetic_models.ppigrf(ipp)
    assert field.magnetic_field.shape == times.shape + heights.shape
    assert field.magnetic_field_error.shape == times.shape + heights.shape
    assert is_convertible_to_unit(field.magnetic_field, u.tesla)
    assert np.isclose(field.magnetic_field[0, 0].value, -5335, 0.5)
