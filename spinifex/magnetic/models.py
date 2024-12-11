#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module for getting the Earth magnetic field"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time
from ppigrf import igrf


class MagneticFieldFunction(Protocol):
    """Magnetic field callable"""

    def __call__(self, loc: EarthLocation, date: Time) -> u.Quantity: ...


@dataclass
class MagneticModels:
    """Supported magnetic field models"""

    ppigrf: MagneticFieldFunction


def get_ppigrf_magnetic_field(loc: EarthLocation, date: Time) -> u.Quantity:
    """Get the magnetic field at a given EarthLocation"""
    lon_deg = loc.lon.to(u.deg).value
    lat_deg = loc.lat.to(u.deg).value
    height_km = loc.height.to(u.km).value
    b_field_east, b_field_north, b_field_up = igrf(
        lon=lon_deg, lat=lat_deg, h=height_km, date=date.to_datetime()
    )
    return u.Quantity((b_field_east, b_field_north, b_field_up) * u.nanotesla)


magnetic_models = MagneticModels(ppigrf=get_ppigrf_magnetic_field)
