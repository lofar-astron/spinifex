"""Testing of tomion ionospheric model"""

from __future__ import annotations

from importlib import resources

from astropy.utils import iers
from spinifex.ionospheric.models import parse_iono_kwargs
from spinifex.ionospheric.tomion_parser import TOMOION_FORMAT_DICT

iers.conf.auto_download = False

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from spinifex.geometry.get_ipp import IPP, get_ipp_from_skycoord
from spinifex.ionospheric import ionospheric_models


@pytest.fixture
def ipp() -> IPP:
    cas_a = SkyCoord(ra=350.85 * u.deg, dec=58.815 * u.deg)
    lon = 6.367 * u.deg
    lat = 52.833 * u.deg
    heights = np.arange(100, 2000, 100) * u.km
    dwingeloo = EarthLocation(lon=lon, lat=lat, height=0 * u.km)
    times = Time("2023-01-08T01:00:00") + np.arange(0, 10) * 15 * u.min
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
    times = Time("2023-01-08T01:00:00") + np.arange(0, 10) * 3.5 * u.hr
    return get_ipp_from_skycoord(
        loc=dwingeloo, times=times, source=cas_a, height_array=heights
    )


def test_ionosphere_tomion(ipp):
    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        options = parse_iono_kwargs(
            ionospheric_models.tomion,
            output_directory=datapath,
        )
        tec = ionospheric_models.tomion(ipp, options=options)
        assert tec.electron_density.shape == ipp.loc.shape
        assert tec.electron_density_error.shape == ipp.loc.shape

        # Test bad arguments
        with pytest.raises(TypeError):
            options = parse_iono_kwargs(ionospheric_models.tomion, bad_arg="bad")


def test_ionosphere_tomionmultiple_days(ipp2):
    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        options = parse_iono_kwargs(
            ionospheric_models.tomion,
            output_directory=datapath,
        )
        tec = ionospheric_models.tomion(ipp2, options=options)
        assert tec.electron_density.shape == ipp2.loc.shape
        assert tec.electron_density_error.shape == ipp2.loc.shape


def test_constants():
    # Previous constants from first implementation
    tomion_format = [
        "mjd",
        "index",
        "value",
        "stddev",
        "type",
        "number_of_observations",
        "height",
        "ra",
        "dec",
        "i",
        "j",
        "k",
        "label",
        "longitude",
        "lst",
        "year",
        "doy",
        "month",
        "dom",
    ]

    data_types = [
        float,
        int,
        float,
        float,
        str,
        int,
        float,
        float,
        float,
        int,
        int,
        int,
        str,
        float,
        float,
        int,
        int,
        int,
        int,
    ]

    assert list(TOMOION_FORMAT_DICT.keys()) == tomion_format
    assert list(TOMOION_FORMAT_DICT.values()) == data_types
