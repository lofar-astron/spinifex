"""Module to calculate electron densities"""

from __future__ import annotations

from typing import Any, NamedTuple

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from numpy.typing import NDArray

from spinifex.geometry import IPP, get_ipp_from_altaz, get_ipp_from_skycoord
from spinifex.get_rm import DEFAULT_IONO_HEIGHT
from spinifex.ionospheric import ModelDensityFunction, ionospheric_models
from spinifex.logger import logger


class DTEC(NamedTuple):
    """object with all electron densities"""

    times: Time
    """time axis"""
    electron_density: NDArray[Any]
    """electron content"""
    airmass: NDArray[Any]
    """conversion from vertical to slant TEC"""
    height: NDArray[Any]
    """array of altitudes (km)"""
    loc: EarthLocation
    """observer location"""


def _get_dtec(
    ipp: IPP,
    iono_model: ModelDensityFunction = ionospheric_models.ionex,
    iono_kwargs: dict[str, Any] | None = None,
) -> DTEC:
    """Get the electron densities for a given set of ionospheric piercepoints

    Parameters
    ----------
    ipp : IPP
        ionospheric piercepoints
    iono_model : ModelDensityFunction, optional
        ionospheric model, by default ionospheric_models.ionex
    iono_kwargs : dict
        options for the ionospheric model, by default {}

    Returns
    -------
    DTEC
        electron densities object
    """
    logger.info("Calculating electron density")
    iono_kwargs = iono_kwargs or {}
    density_profile = iono_model(ipp=ipp, iono_kwargs=iono_kwargs)
    return DTEC(
        times=ipp.times,
        electron_density=density_profile,
        airmass=ipp.airmass,
        height=ipp.loc.height.to(u.km).value,
        loc=ipp.station_loc,
    )


def get_dtec_from_altaz(
    loc: EarthLocation,
    altaz: AltAz,
    height_array: u.Quantity = DEFAULT_IONO_HEIGHT,
    iono_model: ModelDensityFunction = ionospheric_models.ionex,
    iono_kwargs: dict[str, Any] | None = None,
) -> DTEC:
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
    iono_kwargs : dict
        options for the ionospheric model, by default {}

    Returns
    -------
    DTEC
        electron density object
    """
    ipp = get_ipp_from_altaz(loc=loc, altaz=altaz, height_array=height_array)
    return _get_dtec(
        ipp=ipp,
        iono_model=iono_model,
        iono_kwargs=iono_kwargs,
    )


def get_dtec_from_skycoord(
    loc: EarthLocation,
    times: Time,
    source: SkyCoord,
    height_array: u.Quantity = DEFAULT_IONO_HEIGHT,
    iono_model: ModelDensityFunction = ionospheric_models.ionex,
    iono_kwargs: dict[str, Any] | None = None,
) -> DTEC:
    """get electron densities for user defined times and source coordinate

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
    iono_kwargs : dict
        options for the ionospheric model, by default {}


    Returns
    -------
    DTEC
        relectron density object
    """

    ipp = get_ipp_from_skycoord(
        loc=loc, times=times, source=source, height_array=height_array
    )
    return _get_dtec(
        ipp=ipp,
        iono_model=iono_model,
        iono_kwargs=iono_kwargs,
    )
