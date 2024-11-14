#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

""" Module for getting the Earth magnetic field """

__all__ = ["get_magnetic_field"]

from astropy.coordinates import EarthLocation
import astropy.units as u


def get_magnetic_field(loc: EarthLocation) -> u.Quantity:
    """Get the magnetic field at a given EarthLocation"""
    lon_deg = loc.longitude.to(u.deg).value
    lat_deg = loc.latitude.to(u.deg).value
    height_m = loc.height.to(u.m).value

    return u.Quantity((3, 5, 4) * u.tesla)
