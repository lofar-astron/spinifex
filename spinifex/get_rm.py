"""Module to calculate rotation measures"""

from __future__ import annotations

from typing import Any, NamedTuple

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from numpy.typing import NDArray

from spinifex.geometry import IPP, get_ipp_from_altaz, get_ipp_from_skycoord
from spinifex.ionospheric import ModelDensityFunction
from spinifex.ionospheric.iri_density import IRI_HEIGHTS
from spinifex.ionospheric.models import (
    O,
    parse_iono_kwargs,
    parse_iono_model,
)
from spinifex.ionospheric.tomion_parser import TOMION_HEIGHTS
from spinifex.logger import logger
from spinifex.magnetic import MagneticFieldFunction
from spinifex.magnetic.models import parse_magnetic_model

DEFAULT_IONO_HEIGHT = np.array([450.0]) * u.km


class RM(NamedTuple):
    """object with all rotation measures"""

    rm: NDArray[Any]
    """rotation measures"""
    rm_error: NDArray[Any]
    """error on rotation measures"""
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
    iono_model: ModelDensityFunction[O],
    magnetic_model: MagneticFieldFunction,
    iono_options: O | None = None,
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
    iono_options : OptionType, optional
        options for the ionospheric model, by default None

    Returns
    -------
    RM
        rotation measures object
    """
    logger.info("Calculating rotation measure")
    density_profile = iono_model(ipp=ipp, options=iono_options)

    magnetic_profile = magnetic_model(ipp=ipp)
    b_field_to_rm = -2.62e-6  # TODO: What are the units of this constant?

    rm = np.sum(
        b_field_to_rm
        * density_profile.electron_density
        * magnetic_profile.magnetic_field.to(u.nT).value
        * ipp.airmass,
        axis=1,
    )
    relative_uncertainty = np.abs(
        np.sum(density_profile.electron_density_error, axis=1)
        / np.sum(density_profile.electron_density, axis=1)
    ) + np.abs(
        np.sum(magnetic_profile.magnetic_field_error.to(u.nT).value, axis=1)
        / np.sum(magnetic_profile.magnetic_field.to(u.nT).value, axis=1)
    )
    rm_error = relative_uncertainty * np.abs(rm)
    return RM(
        rm=rm,
        rm_error=rm_error,
        times=ipp.times,
        b_parallel=magnetic_profile.magnetic_field,
        electron_density=density_profile.electron_density,
        height=ipp.loc.height.to(u.km).value,
        azimuth=ipp.altaz.az.deg,
        elevation=ipp.altaz.alt.deg,
        loc=ipp.station_loc,
    )


def get_average_rm(rm: RM) -> RM:
    profile_weights = np.sum(rm.electron_density, axis=1, keepdims=True)
    return RM(
        rm=rm.rm.mean(),
        rm_error=np.sqrt(np.mean(rm.rm_error**2)),
        times=rm.times.mean(),
        b_parallel=np.sum(
            rm.b_parallel * rm.electron_density / profile_weights,
            axis=1,
        ).mean(),
        electron_density=profile_weights.mean(),
        height=np.sum(rm.height * rm.electron_density / profile_weights, axis=1).mean(),
        azimuth=np.degrees(np.angle(np.sum(np.exp(1.0j * np.radians(rm.azimuth))))),
        elevation=rm.elevation.mean(),
        loc=rm.loc,
    )


def _get_rm_from_altaz(
    loc: EarthLocation,
    altaz: AltAz,
    iono_model: ModelDensityFunction[O],
    magnetic_model: MagneticFieldFunction,
    height_array: u.Quantity = DEFAULT_IONO_HEIGHT,
    iono_options: O | None = None,
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
        keyword arguments for the ionospheric model

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
        iono_options=iono_options,
    )


def get_rm_from_altaz(
    loc: EarthLocation,
    altaz: AltAz,
    iono_model_name: str = "ionex",
    magnetic_model_name: str = "ppigrf",
    **iono_kwargs: Any,
) -> RM:
    """get rotation measures for user defined altaz coordinates

    Parameters
    ----------
    loc : EarthLocation
        observer location
    altaz : AltAz
        altaz coordinates
    iono_model_name : str, optional
        ionospheric model name, by default "ionex". Must be a supported ionospheric model.
    magnetic_model_name : str, optional
        geomagnetic model name, by default "ppigrf". Must be a supported geomagnetic model.
    iono_kwargs : dict
        keyword arguments for the ionospheric model

    Returns
    -------
    RM
        rotation measure object
    """
    iono_model = parse_iono_model(iono_model_name)
    height_array = DEFAULT_IONO_HEIGHT
    if iono_model_name == "tomion":
        height_array = TOMION_HEIGHTS
    elif iono_model_name == "ionex_iri":
        height_array = IRI_HEIGHTS

    iono_options = parse_iono_kwargs(iono_model=iono_model, **iono_kwargs)
    magnetic_model = parse_magnetic_model(magnetic_model_name)
    return _get_rm_from_altaz(
        loc=loc,
        altaz=altaz,
        height_array=height_array,
        iono_model=iono_model,
        magnetic_model=magnetic_model,
        iono_options=iono_options,
    )


def _get_rm_from_skycoord(
    loc: EarthLocation,
    times: Time,
    source: SkyCoord,
    iono_model: ModelDensityFunction[O],
    magnetic_model: MagneticFieldFunction,
    height_array: u.Quantity = DEFAULT_IONO_HEIGHT,
    iono_options: O | None = None,
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
    iono_options : IonoOptions, optional
        options for the ionospheric model, by default None


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
        iono_options=iono_options,
    )


def get_rm_from_skycoord(
    loc: EarthLocation,
    times: Time,
    source: SkyCoord,
    iono_model_name: str = "ionex",
    magnetic_model_name: str = "ppigrf",
    **iono_kwargs: Any,
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
    iono_model_name : str, optional
        ionospheric model name, by default "ionex". Must be a supported ionospheric model.
    magnetic_model_name : str, optional
        geomagnetic model name, by default "ppigrf". Must be a supported geomagnetic model.
    iono_kwargs : dict
        keyword arguments for the ionospheric model


    Returns
    -------
    RM
        rotation measure object
    """
    iono_model = parse_iono_model(iono_model_name)
    height_array = DEFAULT_IONO_HEIGHT
    if iono_model_name == "tomion":
        height_array = TOMION_HEIGHTS
    elif iono_model_name == "ionex_iri":
        height_array = IRI_HEIGHTS
    iono_options = parse_iono_kwargs(iono_model=iono_model, **iono_kwargs)
    magnetic_model = parse_magnetic_model(magnetic_model_name)
    return _get_rm_from_skycoord(
        loc=loc,
        times=times,
        source=source,
        height_array=height_array,
        iono_model=iono_model,
        magnetic_model=magnetic_model,
        iono_options=iono_options,
    )
