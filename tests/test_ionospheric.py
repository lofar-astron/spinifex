#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0
# pylint: disable=duplicate-code
"""Testing of ionospheric"""

from __future__ import annotations

from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation
from astropy.time import Time
from spinifex.ionospheric import ionospheric_models
from spinifex.ionospheric.ionex_parser import read_ionex


def test_get_ionosphere():
    """Test that get_ionosphere does not crash"""
    datapath = Path(__file__).parent / "data"

    testdata = datapath / "CODG0080.20I"

    ionex = read_ionex(testdata)
    assert ionex.tec.shape == (25, 73, 71)
    times = Time(ionex.times[:2] + 30 * u.min)
    height = np.linspace(100, 1000, 2) * u.km
    lon = np.array([[6.367, 6.5], [3.367, 7.5]]) * u.deg
    lat = np.array([[52.833, 60.0], [52.833, 60.0]]) * u.deg  # h x times
    above_dwingeloo = EarthLocation(lon=lon, lat=lat, height=height[:, np.newaxis])
    tec = ionospheric_models.ionex(loc=above_dwingeloo, times=times)
    assert tec.shape == (2,)
    tec = ionospheric_models.ionex(loc=above_dwingeloo[:, :1], times=times[:1])
    assert tec.shape == (1,)
    tec = ionospheric_models.ionex_iri(loc=above_dwingeloo[:, :], times=times[:])
    assert tec.shape == (2, 2)
