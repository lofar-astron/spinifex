#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0
# pylint: disable=duplicate-code
"""Testing of ionospheric"""

from __future__ import annotations

from importlib import resources

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from spinifex.geometry.get_ipp import get_ipp_from_skycoord
from spinifex.ionospheric import ionospheric_models
from spinifex.ionospheric.ionex_parser import (
    read_ionex,
    unique_days_from_ionex,
    unique_days_from_ionex_files,
)


def test_get_ionosphere():
    """Test that get_ionosphere does not crash"""
    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        testdata = datapath / "CODG0080.20I.gz"

    ionex = read_ionex(testdata)
    assert ionex.tec.shape == (25, 73, 71)
    source = SkyCoord.from_name("Cas A")
    lon = 6.367 * u.deg
    lat = 52.833 * u.deg
    heights = np.arange(100, 2000, 100) * u.km
    dwingeloo = EarthLocation(lon=lon, lat=lat, height=0 * u.km)
    times = Time("2020-01-08T01:00:00") + np.arange(0, 10) * 15 * u.min
    ipp = get_ipp_from_skycoord(
        loc=dwingeloo, times=times, source=source, height_array=heights
    )
    tec = ionospheric_models.ionex(ipp)
    assert tec.shape == ipp.loc.shape
    tec = ionospheric_models.ionex_iri(ipp)
    assert tec.shape == ipp.loc.shape


def test_read_zcompressed():
    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        testdata = datapath / "casg0010.99i.Z"
    ionex = read_ionex(testdata)
    assert ionex.tec.shape == (12, 73, 71)


def test_unique_days():
    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        testdata = datapath / "CODG0080.20I.gz"
        other_data = datapath / "casg0010.99i.Z"

    unique_days = unique_days_from_ionex_files(testdata)

    assert unique_days.value.shape == (1,)
    assert unique_days.value[0] == 58856.0

    ionex = read_ionex(testdata)
    unique_days = unique_days_from_ionex(ionex)
    assert unique_days.value.shape == (1,)
    assert unique_days.value[0] == 58856.0

    unique_days = unique_days_from_ionex_files([testdata, testdata])
    assert unique_days.value.shape == (1,)
    assert unique_days.value[0] == 58856.0

    ionex = read_ionex(testdata)
    unique_days = unique_days_from_ionex([ionex, ionex])
    assert unique_days.value.shape == (1,)
    assert unique_days.value[0] == 58856.0

    unique_days = unique_days_from_ionex_files([testdata, other_data])
    assert unique_days.value.shape == (2,)
    assert unique_days.value[0] == 51179.0
    assert unique_days.value[1] == 58856.0
