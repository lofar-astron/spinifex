#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module for getting the Earth magnetic field"""

from __future__ import annotations

import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from ppigrf import igrf


def get_ppigrf_magnetic_field(loc: EarthLocation, date: Time) -> u.Quantity:
    """Get the magnetic field at a given EarthLocation"""
    lon_deg = loc.lon.to(u.deg).value
    lat_deg = loc.lat.to(u.deg).value
    height_km = loc.height.to(u.km).value
    be, bn, bu = igrf(lon=lon_deg, lat=lat_deg, h=height_km, date=date.to_datetime())
    return u.Quantity((be, bn, bu) * u.nanotesla)
