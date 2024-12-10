#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module for getting the Earth magnetic field"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import astropy.units as u
import numpy as np
from astropy.coordinates import ITRS, AltAz
from ppigrf import igrf

from spinifex.geometry.get_ipp import IPP


class MagneticFieldFunction(Protocol):
    """Magnetic field callable"""

    def __call__(self, IPP) -> u.Quantity: ...


@dataclass
class MagneticModels:
    """Supported magnetic field models"""

    ppigrf: MagneticFieldFunction


def get_ppigrf_magnetic_field(ipp: IPP) -> u.Quantity:
    """Get the magnetic field at a given EarthLocation"""
    loc = ipp.loc
    date = ipp.times[0]
    b_e, b_n, b_u = igrf(
        lon=loc.lon.deg,
        lat=loc.lat.deg,
        h=loc.height.to(u.km).value,
        date=date.to_datetime(),
    )
    # ppigrf adds an extra axis for time, we remove it by taking the first element
    b_magn = np.sqrt(b_e**2 + b_n**2 + b_u**2)[0]
    b_az = np.arctan2(b_e, b_n)
    b_el = np.arctan2(b_u, np.sqrt(b_n**2 + b_e**2))
    b_altaz = AltAz(az=b_az[0] * u.rad, alt=b_el[0] * u.rad, location=loc)
    b_itrs = b_altaz.transform_to(ITRS())

    # project to LOS
    los = ipp.los
    b_par = los.x * b_itrs.x + los.y * b_itrs.y + los.z * b_itrs.z
    b_par = b_par * b_magn
    # magnitude along LOS,
    return u.Quantity(b_par * u.nanotesla)


magnetic_models = MagneticModels(ppigrf=get_ppigrf_magnetic_field)
