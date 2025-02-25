"""Module to calculate rotation measures"""

from __future__ import annotations

from typing import Any, NamedTuple

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from numpy.typing import NDArray

from spinifex.geometry import IPP, get_ipp_from_altaz, get_ipp_from_skycoord
from spinifex.ionospheric import ModelDensityFunction, ionospheric_models
from spinifex.logger import logger
from spinifex.magnetic import MagneticFieldFunction, magnetic_models

DEFAULT_IONO_HEIGHT = np.array([450.0]) * u.km


class RM(NamedTuple):
    """object with all rotation measures"""

    rm: NDArray[Any]
    """rotation measures"""
    times: Time
    """time axis"""
    b_parallel: NDArray[Any]
    """parallel magnetic field"""
    electron_density: NDArray[Any]
    """electron content"""
    height: NDArray[Any]
    """array of altitudes (km)"""
    azimuth: NDArray[Any]
    """array of azimuths (degrees)"""
    elevation: NDArray[Any]
    """array of elevation (degrees)"""
    loc: EarthLocation
    """observer location"""


def _get_rm(
    ipp: IPP,
    iono_model: ModelDensityFunction = ionospheric_models.ionex,
    magnetic_model: MagneticFieldFunction = magnetic_models.ppigrf,
    iono_kwargs: dict[str, Any] | None = None,
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
    logger.info("Calculating rotation measure")
    iono_kwargs = iono_kwargs or {}
    density_profile = iono_model(ipp=ipp, iono_kwargs=iono_kwargs)
    magnetic_profile = magnetic_model(ipp=ipp)
    b_field_to_rm = -2.62e-6  # TODO: What are the units of this constant?
    rm = np.sum(
        b_field_to_rm * density_profile * magnetic_profile.to(u.nT).value * ipp.airmass,
        axis=1,
    )
    return RM(
        rm=rm,
        times=ipp.times,
        b_parallel=magnetic_profile,
        electron_density=density_profile,
        height=ipp.loc.height.to(u.km).value,
        azimuth=ipp.altaz.az.deg,
        elevation=ipp.altaz.alt.deg,
        loc=ipp.station_loc,
    )


def get_average_rm(rm: RM) -> RM:
    profile_weights = np.sum(rm.electron_density, axis=1, keepdims=True)
    return RM(
        rm=rm.rm.mean(),
        times=rm.times.mean(),
        b_parallel=np.sum(
            rm.b_parallel * rm.electron_density / profile_weights, axis=1
        ).mean(),
        electron_density=profile_weights.mean(),
        height=np.sum(rm.height * rm.electron_density / profile_weights, axis=1).mean(),
        azimuth=np.degrees(np.angle(np.sum(np.exp(1.0j * np.radians(rm.azimuth))))),
        elevation=rm.elevation.mean(),
    )


def get_rm_from_altaz(
    loc: EarthLocation,
    altaz: AltAz,
    height_array: u.Quantity = DEFAULT_IONO_HEIGHT,
    iono_model: ModelDensityFunction = ionospheric_models.ionex,
    magnetic_model: MagneticFieldFunction = magnetic_models.ppigrf,
    iono_kwargs: dict[str, Any] | None = None,
) -> RM:
    """get rotation measures for user defined altaz coordinates

    Parameters
    ----------
    loc : EarthLocation
        observer location
    altaz : AltAz
        altaz coordinates
    height_array : u.Quantity, optional
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
    height_array: u.Quantity = DEFAULT_IONO_HEIGHT,
    iono_model: ModelDensityFunction = ionospheric_models.ionex,
    magnetic_model: MagneticFieldFunction = magnetic_models.ppigrf,
    iono_kwargs: dict[str, Any] | None = None,
) -> RM:
    """get rotation measures for user defined times and source coordinate

    Parameters
    ----------
    loc : EarthLocation
        observer location
    times : Time
        times
    source : SkyCoord
        coordinates of the source
    height_array : NDArray, optional
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
