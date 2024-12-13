from __future__ import annotations

from typing import NamedTuple

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from numpy.typing import ArrayLike

from spinifex.geometry.get_ipp import IPP, get_ipp_from_altaz, get_ipp_from_skycoord
from spinifex.ionospheric.models import ModelDensityFunction, ionospheric_models
from spinifex.magnetic.models import MagneticFieldFunction, magnetic_models

default_height = np.array([450.0]) * u.km


class RM(NamedTuple):
    """object with all rotation measures"""

    rm: ArrayLike
    """rotation measures"""
    times: Time
    """time axis"""
    b_parallel: ArrayLike
    """parallel magnetic field"""
    electron_density: ArrayLike
    """electron content"""
    height: ArrayLike
    """array of altitudes (km)"""
    azimuth: ArrayLike
    """array of azimuths (degrees)"""
    elevation: ArrayLike
    """array of elevation (degrees)"""


def _get_rm(
    ipp: IPP,
    iono_model: ModelDensityFunction = ionospheric_models.ionex,
    magnetic_model: MagneticFieldFunction = magnetic_models.ppigrf,
    iono_kwargs: dict | None = None,
) -> RM:
    """Get the rotation measures for a given set of ionospheric piercepoints

    Parameters
    ----------
    ipp : IPP
        ionospheric piercepoints
    iono_model : ModelDensityFunction, optional
        ionospheric model, by default ionospheric_models.ionex
    magnetic_model : MagneticFieldFunction, optional
        geomagnetic model, by default magnetic_models.ppigrf
    iono_kwargs : dict
        options for the ionospheric model, by default {}

    Returns
    -------
    RM
        rotation measures object
    """
    iono_kwargs = iono_kwargs or {}
    density_profile = iono_model(ipp=ipp, **iono_kwargs)
    magnetic_profile = magnetic_model(ipp=ipp)
    rm = np.sum(
        -2.62e-6 * density_profile * magnetic_profile.to(u.nT).value * ipp.airmass,
        axis=0,
    )
    return RM(
        rm=rm,
        times=ipp.times,
        b_parallel=magnetic_profile,
        electron_density=density_profile,
        height=ipp.loc[:, 0].height.to(u.km).value,
        azimuth=ipp.altaz.az.deg,
        elevation=ipp.altaz.alt.deg,
    )


def get_rm_from_altaz(
    loc: EarthLocation,
    altaz: AltAz,
    height_array: ArrayLike = default_height,
    iono_model: ModelDensityFunction = ionospheric_models.ionex,
    magnetic_model: MagneticFieldFunction = magnetic_models.ppigrf,
    iono_kwargs: dict | None = None,
) -> RM:
    """get rotation measures for user defined altaz coordinates

    Parameters
    ----------
    loc : EarthLocation
        observer location
    altaz : AltAz
        altaz coordinates
    height_array : ArrayLike, optional
        altitudes, by default default_height
    iono_model : ModelDensityFunction, optional
        ionospheric model, by default ionospheric_models.ionex
    magnetic_model : MagneticFieldFunction, optional
        geomagnetic model, by default magnetic_models.ppigrf
    iono_kwargs : dict
        options for the ionospheric model, by default {}
    Returns
    -------
    RM
        rotation measure object
    """
    ipp = get_ipp_from_altaz(loc=loc, altaz=altaz, height_array=height_array)
    return _get_rm(
        ipp=ipp,
        iono_model=iono_model,
        magnetic_model=magnetic_model,
        iono_kwargs=iono_kwargs,
    )


def get_rm_from_skycoord(
    loc: EarthLocation,
    times: Time,
    source: SkyCoord,
    height_array=default_height,
    iono_model=ionospheric_models.ionex,
    magnetic_model=magnetic_models.ppigrf,
    iono_kwargs: dict | None = None,
):
    """get rotation measures for user defined times and source coordinate

    Parameters
    ----------
    loc : EarthLocation
        observer location
    times : Time
        times
    source : SkyCoord
        coordinates of the source
    height_array : ArrayLike, optional
        altitudes, by default default_height
    iono_model : ModelDensityFunction, optional
        ionospheric model, by default ionospheric_models.ionex
    magnetic_model : MagneticFieldFunction, optional
        geomagnetic model, by default magnetic_models.ppigrf
    iono_kwargs : dict
        options for the ionospheric model, by default {}


    Returns
    -------
    RM
        rotation measure object
    """

    ipp = get_ipp_from_skycoord(
        loc=loc, times=times, source=source, height_array=height_array
    )
    return _get_rm(
        ipp=ipp,
        iono_model=iono_model,
        magnetic_model=magnetic_model,
        iono_kwargs=iono_kwargs,
    )
