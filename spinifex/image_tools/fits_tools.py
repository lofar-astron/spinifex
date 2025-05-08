"""Predict / correct ionospheric RM in FITS images"""

# Follows report by Van Eck (2021) and
# Python implementation in FRion https://github.com/CIRADA-Tools/FRion/
# FRion License: MIT Copyright (c) 2021 Cameron Van Eck

from __future__ import annotations

from pathlib import Path
from typing import Any, NamedTuple

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS

from spinifex.exceptions import FITSHeaderError
from spinifex.get_rm import RM, get_rm_from_skycoord
from spinifex.image_tools.image_tools import IntegratedRM, get_integrated_rm
from spinifex.logger import logger


class FITSMetaData(NamedTuple):
    """Metadata from a FITS file"""

    start_time: Time
    """Start time of obsrvation"""
    duration: TimeDelta
    """Duration of observation"""
    location: EarthLocation
    """Location of observation"""
    name: str
    """Name of FITS file as str"""
    source: SkyCoord
    """Target of observation"""


def get_metadata_from_fits(fits_path: Path) -> FITSMetaData:
    """Get metadata from a FITS file header

    Parameters
    ----------
    fits_path : Path
        FITS file

    Returns
    -------
    FitsMetaData
        FITS image metadata
    """

    header = fits.getheader(fits_path)
    wcs = WCS(header)
    wcs_celestial = wcs.celestial
    source = SkyCoord(wcs_celestial.pixel_to_world(*wcs_celestial.wcs.crpix))

    # Lots of try/excepts/if/elses here - I'm sorry
    # This is a bit gross, but I'm being super paranoid

    if "MJD-OBS" in header:
        start_time = Time(header["MJD-OBS"], format="mjd")

    else:
        try:
            start_time = Time(
                header["DATE-OBS"],
                format="fits",
                scale="utc",
            )
        except KeyError as key_error:
            msg = "Keys `DATE-OBS` or `MJD-OBS` not in the header - we can't get the time of the observation"
            raise FITSHeaderError(msg) from key_error
    try:
        duration = TimeDelta(header["DURATION"], format="sec")
    except KeyError as key_error:
        msg = "`DURATION` not in the header - we can't get the time of the observation"
        raise FITSHeaderError(msg) from key_error

    # Get telescope location
    if "TELESCOP" in header:
        msg = "Getting telescope location from `TELESCOP` key"
        logger.info(msg)
        location = EarthLocation.of_site(site_name=header["TELESCOP"])

    elif all(key in header for key in ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z")):
        msg = "Getting telescope location from `OBSGEO` keys"
        logger.info(msg)
        location = EarthLocation(
            x=header["OBSGEO-X"], y=header["OBSGEO-Y"], z=header["OBSGEO-Z"]
        )
    elif all(key in header for key in ("ALT-OBS", "LAT-OBS", "LONG-OBS")):
        msg = "Getting telescope location from `ALT-,LAT-,LONG-OBS` keys"
        logger.info(msg)
        location = EarthLocation(
            lon=header["LONG-LOBS"], lat=header["LAT-OBS"], height=header["ALT-OBS"]
        )
    else:
        msg = "Cannot determine telescope coordinates from FITS header"
        raise FITSHeaderError(msg)

    return FITSMetaData(
        start_time=start_time,
        duration=duration,
        location=location,
        name=fits_path.name,
        source=source,
    )


def get_rm_from_fits(
    fits_path: Path,
    timestep: u.Quantity = 15 * u.min,
    iono_model_name: str = "ionex",
    magnetic_model_name: str = "ppigrf",
    **iono_kwargs: Any,
) -> RM:
    """Get the ionospheric RM from a FITS file

    Parameters
    ----------
    fits_path : Path
        Path to FITS file
    timestep : u.Quantity, optional
        Time step to evaluate time-dependant RM, by default 15*u.min
    iono_model_name : str, optional
        ionospheric model name, by default "ionex". Must be a supported ionospheric model.
    magnetic_model_name : str, optional
        geomagnetic model name, by default "ppigrf". Must be a supported geomagnetic model.
    iono_kwargs : dict
        Keyword arguments for the ionospheric model

    Returns
    -------
    RM
        Rotation measure object
    """
    fits_metadata = get_metadata_from_fits(fits_path=fits_path)

    start_time = fits_metadata.start_time
    duration = fits_metadata.duration
    end_time = start_time + duration
    times_mjd = np.arange(
        start=start_time.mjd, stop=end_time.mjd, step=timestep.to(u.day).value
    )
    times = Time(times_mjd, format="mjd")

    return get_rm_from_skycoord(
        loc=fits_metadata.location,
        times=times,
        source=fits_metadata.source,
        iono_model_name=iono_model_name,
        magnetic_model_name=magnetic_model_name,
        **iono_kwargs,
    )


def get_freq_from_fits(fits_path: Path) -> u.Quantity:
    """Get frequency array from FITS file

    Parameters
    ----------
    fits_path : Path
        Path to FITS file

    Returns
    -------
    u.Quantity
        Frequency array

    Raises
    ------
    FITSHeaderError
        If no spectral axis is found
    """
    header = fits.getheader(fits_path)

    if "SPECSYS" not in header:
        msg = "Could not find `SPECSYS` in header, assuming `TOPOCENT`"
        logger.warning(msg)
        header["SPECSYS"] = "TOPOCENT"

    wcs = WCS(header)
    if not wcs.has_spectral:
        msg = "FITS header has no spectral axis"
        raise FITSHeaderError(msg)

    specral_wcs = wcs.spectral

    return specral_wcs.pixel_to_world(np.arange(*specral_wcs.array_shape))


def get_integrated_rm_from_fits(
    fits_path: Path,
    timestep: u.Quantity = 15 * u.min,
    iono_model_name: str = "ionex",
    magnetic_model_name: str = "ppigrf",
    **iono_kwargs: Any,
) -> IntegratedRM:
    """Computed the integrated RM effect following Van Eck (2021)

    Parameters
    ----------
    fits_path : Path
        Path to FITS file
    timestep : u.Quantity, optional
        Timestep to use for computing time-dependent RM, by default 15*u.min
    iono_model_name : str, optional
        ionospheric model name, by default "ionex". Must be a supported ionospheric model.
    magnetic_model_name : str, optional
        geomagnetic model name, by default "ppigrf". Must be a supported geomagnetic model.
    iono_kwargs : dict
        keyword arguments for the ionospheric model

    Returns
    -------
    IntegratedRM
        Integrated rotation measure object
    """

    # Get the time-dependent RM
    time_dep_rm = get_rm_from_fits(
        fits_path=fits_path,
        timestep=timestep,
        iono_model_name=iono_model_name,
        magnetic_model_name=magnetic_model_name,
        **iono_kwargs,
    )

    freq_arr = get_freq_from_fits(fits_path)

    return get_integrated_rm(
        time_dep_rm=time_dep_rm,
        freq_arr=freq_arr,
    )


def correct_fits_images(
    stokes_q_path: Path,
    stokes_u_path: Path,
    integrated_rm: IntegratedRM,
    ext: int = 0,
    suffix: str = "fr_corr",
) -> None:
    msg = "TODO: Implement correct_fits_images"
    raise NotImplementedError(msg)
