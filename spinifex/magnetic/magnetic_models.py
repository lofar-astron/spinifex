#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module for getting the Earth magnetic field"""

from __future__ import annotations

import astropy.units as u
import numpy as np
from ppigrf import igrf_gc

from spinifex.geometry.get_ipp import IPP


def get_ppigrf_magnetic_field(ipp: IPP) -> u.Quantity:
    """Get the magnetic field at a given EarthLocation"""
    loc = ipp.loc
    phi = np.arctan2(loc.x, loc.y)
    theta = np.arctan2(loc.z, np.sqrt(loc.x**2 + loc.y**2))
    r = np.sqrt(loc.x**2 + loc.y**2 + loc.z**2)
    # possibly these conversions are easier in astropy
    date = ipp.times[0]
    b_r, b_theta, b_phi = igrf_gc(
        r=r.to(u.km).value,
        theta=theta.to(u.deg).value,
        phi=phi.to(u.deg).value,
        date=date.to_datetime(),
    )
    costheta = np.cos(np.radians(b_theta))
    sintheta = np.sin(np.radians(b_theta))
    cosphi = np.cos(np.radians(b_phi))
    sinphi = np.sin(np.radians(b_phi))
    b_xyz = np.array([b_r * sinphi * costheta, b_r * cosphi * costheta, b_r * sintheta])

    # project to LOS
    los = ipp.los
    b_par = (
        los.x * b_xyz[0] + los.y * b_xyz[1] + los.z * b_xyz[2]
    )  # magnitude along LOS
    return u.Quantity(b_par * u.nanotesla)
