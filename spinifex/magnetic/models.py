"""Module for getting the Earth magnetic field"""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple, Protocol

import astropy.units as u
import numpy as np
from ppigrf import igrf_gc

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


def _transform_b_gc_to_itrs(
    b_r: float, b_theta: float, b_phi: float, lon: float, lat: float
) -> tuple[float, float, float]:
    """Helper function to transform local b_field in geocentric angles to itrs x y z

    Parameters
    ----------
    b_r : float
        radial b
    b_theta : float
        latitude b
    b_phi : float
        longitude b
    lon : float
        longitude of ipp (radians)
    lat : float
        latitude of ipp (radians)

    Returns
    -------
    tuple[float,float, float]
        x, y, z direction (not normalized) of b_field
    """
    theta_rad = np.pi / 2 - lat  # colatitude
    sin_theta = np.sin(theta_rad)
    cos_theta = np.cos(theta_rad)
    sin_phi = np.sin(lon)
    cos_phi = np.cos(lon)

    b_x = sin_theta * cos_phi * b_r + cos_theta * cos_phi * b_theta - sin_phi * b_phi

    b_y = sin_theta * sin_phi * b_r + cos_theta * sin_phi * b_theta + cos_phi * b_phi

    b_z = cos_theta * b_r - sin_theta * b_theta

    return b_x, b_y, b_z


def get_ppigrf_magnetic_field(ipp: IPP) -> MagneticProfile:
    """Get the magnetic field at a given EarthLocation"""

    RMS_E = 87
    RMS_N = 73
    RMS_U = 114

    # constants from https://geomag.bgs.ac.uk/research/modelling/IGRF.html

    unique_days = get_unique_days(ipp.times)
    b_par = np.zeros(ipp.lon.shape, dtype=float)
    relative_uncertainty = np.zeros_like(b_par)

    # ppigrf uses proper geodetic coordinates,
    # use Br, Btheta, Bphi = ppigrf.igrf_gc(r, theta, phi, date)

    for u_day in unique_days:
        indices = np.floor(ipp.times.mjd) == np.floor(u_day.mjd)
        # loc = ipp.loc[indices]
        b_r, b_theta, b_phi = igrf_gc(
            theta=ipp.lat.to(u.deg).value,
            phi=ipp.lon.to(u.deg).value,
            h=ipp.height.to(u.km).value,  # geocentric radius in km
            date=u_day.to_datetime(),
        )
        # ppigrf adds an extra axis for time, we remove it by taking the first element
        b_magn = np.sqrt(b_r**2 + b_theta**2 + b_phi**2)[0]
        # relative uncertainty is 1/2 of the relative uncertainty (rms / b_magn**2) of the
        # sum of individual uncertainties of the squares (2 * rms_<enu> * b_<enu>)
        # multiply by b_par to get absolute value
        rms = (RMS_E * b_phi + RMS_N * b_theta + RMS_U * b_r) / (b_magn)
        relative_uncertainty[indices] = rms / b_magn

        b_x, b_y, b_z = _transform_b_gc_to_itrs(
            b_r, b_theta, b_phi, ipp.lon.rad, ipp.lat.rad
        )

        # project to LOS
        los = ipp.los[indices][:, np.newaxis]
        b_par[indices] = los[:, :, 0] * b_x + los[:, :, 1] * b_y + los[:, :, 2] * b_z
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
