"""Module to calculate electron densities"""

from __future__ import annotations

from typing import Any, NamedTuple

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from numpy.typing import NDArray

from spinifex.geometry import IPP, get_ipp_from_altaz, get_ipp_from_skycoord
from spinifex.get_rm import DEFAULT_IONO_HEIGHT
from spinifex.ionospheric import ModelDensityFunction
from spinifex.ionospheric.iri_density import IRI_HEIGHTS
from spinifex.ionospheric.models import (
    O,
    parse_iono_kwargs,
    parse_iono_model,
)
from spinifex.ionospheric.tomion_parser import TOMION_HEIGHTS
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
    iono_model: ModelDensityFunction[O],
    iono_options: O | None = None,
) -> DTEC:
    """Get the electron densities for a given set of ionospheric piercepoints

    Parameters
    ----------
    ipp : IPP
        ionospheric piercepoints
    iono_model : ModelDensityFunction, optional
        ionospheric model, by default ionospheric_models.ionex
    iono_options : OptionType, optional
        options for the ionospheric model, by default None

    Returns
    -------
    DTEC
        electron densities object
    """
    logger.info("Calculating electron density")
    density_profile = iono_model(ipp=ipp, options=iono_options)
    return DTEC(
        times=ipp.times,
        electron_density=density_profile.electron_density,
        airmass=ipp.airmass,
        height=ipp.loc.height.to(u.km).value,
        loc=ipp.station_loc,
    )


def _get_dtec_from_altaz(
    loc: EarthLocation,
    altaz: AltAz,
    iono_model: ModelDensityFunction[O],
    height_array: u.Quantity = DEFAULT_IONO_HEIGHT,
    iono_options: O | None = None,
) -> DTEC:
    """get electron densities for user defined altaz coordinates

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
    iono_options : dict
        keyword arguments for the ionospheric model


    Returns
    -------
    DTEC
        electron density object
    """
    ipp = get_ipp_from_altaz(loc=loc, altaz=altaz, height_array=height_array)
    return _get_dtec(
        ipp=ipp,
        iono_model=iono_model,
        iono_options=iono_options,
    )


def get_dtec_from_altaz(
    loc: EarthLocation,
    altaz: AltAz,
    iono_model_name: str = "ionex",
    **iono_kwargs: Any,
) -> DTEC:
    """get rotation measures for user defined altaz coordinates

    Parameters
    ----------
    loc : EarthLocation
        observer location
    altaz : AltAz
        altaz coordinates
    iono_model_name : str, optional
        ionospheric model name, by default "ionex". Must be a supported ionospheric model.
    iono_options : dict
        keyword arguments for the ionospheric model


    Returns
    -------
    DTEC
        electron density object
    """
    iono_model = parse_iono_model(iono_model_name)
    height_array = DEFAULT_IONO_HEIGHT
    if iono_model_name == "tomion":
        height_array = TOMION_HEIGHTS
    elif iono_model_name == "ionex_iri":
        height_array = IRI_HEIGHTS

    iono_options = parse_iono_kwargs(iono_model=iono_model, **iono_kwargs)
    return _get_dtec_from_altaz(
        loc=loc,
        altaz=altaz,
        height_array=height_array,
        iono_model=iono_model,
        iono_options=iono_options,
    )


def _get_dtec_from_skycoord(
    loc: EarthLocation,
    times: Time,
    source: SkyCoord,
    iono_model: ModelDensityFunction[O],
    height_array: u.Quantity = DEFAULT_IONO_HEIGHT,
    iono_options: O | None = None,
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
        keyword arguments for the ionospheric model


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
        iono_options=iono_options,
    )


def get_dtec_from_skycoord(
    loc: EarthLocation,
    times: Time,
    source: SkyCoord,
    iono_model_name: str = "ionex",
    # iono_model: ModelDensityFunction[IonoOptions] = ionospheric_models.ionex,
    **iono_kwargs: Any,
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
    iono_model_name : str, optional
        ionospheric model name, by default "ionex". Must be a supported ionospheric model.
    iono_kwargs : dict
        keyword arguments for the ionospheric model


    Returns
    -------
    DTEC
        relectron density object
    """
    iono_model = parse_iono_model(iono_model_name)
    height_array = DEFAULT_IONO_HEIGHT
    if iono_model_name == "tomion":
        height_array = TOMION_HEIGHTS
    elif iono_model_name == "ionex_iri":
        height_array = IRI_HEIGHTS
    iono_options = parse_iono_kwargs(iono_model=iono_model, **iono_kwargs)
    return _get_dtec_from_skycoord(
        loc=loc,
        times=times,
        source=source,
        height_array=height_array,
        iono_model=iono_model,
        iono_options=iono_options,
    )
