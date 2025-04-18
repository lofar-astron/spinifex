"""Testing of ionospheric"""

from __future__ import annotations

from importlib import resources
from typing import Any

from astropy.utils import iers
from spinifex.ionospheric.models import parse_iono_kwargs

iers.conf.auto_download = False

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from spinifex.geometry.get_ipp import IPP, get_ipp_from_skycoord
from spinifex.ionospheric import ionospheric_models
from spinifex.ionospheric.ionex_parser import (
    read_ionex,
    unique_days_from_ionex,
    unique_days_from_ionex_files,
)


@pytest.fixture
def ipp() -> IPP:
    cas_a = SkyCoord(ra=350.85 * u.deg, dec=58.815 * u.deg)
    lon = 6.367 * u.deg
    lat = 52.833 * u.deg
    heights = np.arange(100, 2000, 100) * u.km
    dwingeloo = EarthLocation(lon=lon, lat=lat, height=0 * u.km)
    times = Time("2020-01-08T01:00:00") + np.arange(0, 10) * 15 * u.min
    return get_ipp_from_skycoord(
        loc=dwingeloo, times=times, source=cas_a, height_array=heights
    )


@pytest.fixture
def ipp2() -> IPP:
    cas_a = SkyCoord(ra=350.85 * u.deg, dec=58.815 * u.deg)
    lon = 6.367 * u.deg
    lat = 52.833 * u.deg
    heights = np.arange(100, 2000, 100) * u.km
    dwingeloo = EarthLocation(lon=lon, lat=lat, height=0 * u.km)
    times = Time("2020-01-08T01:00:00") + np.arange(0, 10) * 3.5 * u.hr
    return get_ipp_from_skycoord(
        loc=dwingeloo, times=times, source=cas_a, height_array=heights
    )


def test_get_ionosphere(ipp):
    """Test that get_ionosphere does not crash"""
    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        testdata = datapath / "codg0080.20i.Z"

    ionex = read_ionex(testdata)
    assert ionex.tec.shape == (25, 73, 71)

    tec = ionospheric_models.ionex(ipp)
    assert tec.electron_density.shape == ipp.loc.shape
    tec = ionospheric_models.ionex_iri(ipp)
    assert tec.electron_density.shape == ipp.loc.shape


def test_read_zcompressed():
    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        testdata = datapath / "casg0010.99i.Z"
    ionex = read_ionex(testdata)
    assert ionex.tec.shape == (12, 73, 71)


def test_ionosphere_ionex(ipp):
    iono_kwargs: dict[str, Any] = {}
    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        iono_kwargs["output_directory"] = datapath
        iono_kwargs["prefix"] = "esa"
        iono_kwargs["server"] = "cddis"
        options = parse_iono_kwargs(ionospheric_models.ionex, **iono_kwargs)
        tec = ionospheric_models.ionex(ipp, options=options)
        assert tec.electron_density.shape == ipp.loc.shape
        assert tec.electron_density_error.shape[0] == ipp.loc.shape[0]

        # Test bad arguments
        with pytest.raises(TypeError):
            options = parse_iono_kwargs(ionospheric_models.ionex, bad_arg="bad")


def test_ionosphere_ionex_multiple_days(ipp2):
    iono_kwargs: dict[str, Any] = {}
    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        iono_kwargs["output_directory"] = datapath
        iono_kwargs["prefix"] = "esa"
        iono_kwargs["server"] = "cddis"
        options = parse_iono_kwargs(ionospheric_models.ionex, **iono_kwargs)
        tec = ionospheric_models.ionex(ipp2, options=options)
        assert tec.electron_density.shape == ipp2.loc.shape
        assert tec.electron_density_error.shape[0] == ipp2.loc.shape[0]


def test_unique_days():
    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        testdata = datapath / "codg0080.20i.Z"
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
