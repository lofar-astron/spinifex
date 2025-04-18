"""Module for getting the Earth magnetic field"""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple, Protocol

import astropy.units as u
import numpy as np
from astropy.coordinates import ITRS, AltAz
from ppigrf import igrf

from spinifex.geometry.get_ipp import IPP
from spinifex.times import get_unique_days


class MagneticProfile(NamedTuple):
    """Data object to hold Magnetic field profile and uncertainties"""

    magnetic_field: u.Quantity
    magnetic_field_error: u.Quantity


class MagneticFieldFunction(Protocol):
    """Magnetic field callable"""

    def __call__(self, ipp: IPP) -> MagneticProfile: ...


@dataclass
class MagneticModels:
    """Supported magnetic field models"""

    ppigrf: MagneticFieldFunction


def get_ppigrf_magnetic_field(ipp: IPP) -> u.Quantity:
    """Get the magnetic field at a given EarthLocation"""

    RMS_E = 87
    RMS_N = 73
    RMS_U = 114

    # constants from https://geomag.bgs.ac.uk/research/modelling/IGRF.html

    unique_days = get_unique_days(ipp.times)
    b_par = np.zeros(ipp.loc.shape, dtype=float)
    relative_uncertainty = np.zeros_like(b_par)
    for u_day in unique_days:
        indices = np.floor(ipp.times.mjd) == np.floor(u_day.mjd)
        loc = ipp.loc[indices]
        b_e, b_n, b_u = igrf(
            lon=loc.lon.deg,
            lat=loc.lat.deg,
            h=loc.height.to(u.km).value,
            date=u_day.to_datetime(),
        )
        # ppigrf adds an extra axis for time, we remove it by taking the first element
        b_magn = np.sqrt(b_e**2 + b_n**2 + b_u**2)[0]
        # relative uncertainty is 1/2 of the relative uncertainty (rms / b_magn**2) of the
        # sum of individual uncertainties of the squares (2 * rms_<enu> * b_<enu>)
        # multiply by b_magn to get absolute value
        rms = (RMS_E * b_e + RMS_N * b_n + RMS_U * b_u) / (b_magn)
        relative_uncertainty[indices] = rms / b_magn
        b_az = np.arctan2(b_e, b_n)
        b_el = np.arctan2(b_u, np.sqrt(b_n**2 + b_e**2))
        b_altaz = AltAz(az=b_az[0] * u.rad, alt=b_el[0] * u.rad, location=loc)
        b_itrs = b_altaz.transform_to(ITRS())

        # project to LOS
        los = ipp.los[indices][:, np.newaxis]
        b_par[indices] = los.x * b_itrs.x + los.y * b_itrs.y + los.z * b_itrs.z
        b_par[indices] *= b_magn
        # magnitude along LOS,

    return MagneticProfile(
        magnetic_field=u.Quantity(b_par * u.nanotesla),
        magnetic_field_error=u.Quantity(
            np.abs(b_par * relative_uncertainty) * u.nanotesla
        ),
    )


magnetic_models = MagneticModels(ppigrf=get_ppigrf_magnetic_field)


def parse_magnetic_model(magnetic_model_name: str) -> MagneticFieldFunction:
    """parse magnetic model name

    Parameters
    ----------
    magnetic_model_name : str
        name of the magnetic model

    Returns
    -------
    MagneticFieldFunction
        magnetic field function

    Raises
    ------
    TypeError
        if the magnetic model is not known

    """

    try:
        return getattr(magnetic_models, magnetic_model_name)  # type: ignore[no-any-return]
    except AttributeError as e:
        msg = f"Unknown magnetic model {magnetic_model_name}. Supported models are {list(magnetic_models.__annotations__.keys())}"
        raise TypeError(msg) from e
