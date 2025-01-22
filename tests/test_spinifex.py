#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0
# pylint: disable=duplicate-code
"""Testing of the main module"""

from __future__ import annotations

from importlib import resources
from typing import Any

from astropy.utils import iers  # pylint: disable=ungrouped-imports

iers.conf.auto_download = False

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from spinifex import get_rm
from spinifex.ionospheric import ionospheric_models


def test_get_rm():
    """Test that get_rm does not crash"""
    times = Time("2020-01-08T01:00:00") + np.arange(10) * 25 * u.min
    cas_a = SkyCoord(ra=350.85 * u.deg, dec=58.815 * u.deg)
    lon = 6.367 * u.deg
    lat = 52.833 * u.deg
    dwingeloo = EarthLocation(lon=lon, lat=lat, height=0 * u.km)
    iono_kwargs: dict[str, Any] = {}
    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        iono_kwargs["output_directory"] = datapath
        iono_kwargs["prefix"] = "cod"
        iono_kwargs["server"] = "cddis"
        rm = get_rm.get_rm_from_skycoord(
            loc=dwingeloo,
            times=times,
            source=cas_a,
            iono_model=ionospheric_models.ionex,
            iono_kwargs=iono_kwargs,
        )
        assert isinstance(rm.rm, np.ndarray)
        assert rm.rm.shape == times.shape
        assert np.isclose(rm.rm[0], 0.0687, 0.001)
