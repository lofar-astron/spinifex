"""Module for getting the Ionospheric Piercepoints"""

from __future__ import annotations

from typing import NamedTuple

import astropy.units as u
import numpy as np
from astropy.constants import R_earth
from astropy.coordinates import ITRS, AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from numpy.typing import NDArray

from spinifex.logger import logger


class IPP(NamedTuple):
    """Ionospheric Piercepoints"""

    loc: EarthLocation
    """location of the piercepoints, dimension: times x altitudes. All altitudes are assumed to be equal"""
    times: Time
    """array of times"""
    los: SkyCoord
    """Line of sight direction in ITRS coordinates"""
    airmass: NDArray[np.float64]
    """airmass factor to convert to slant values"""
    altaz: AltAz
    """azimuth elevation"""
    station_loc: EarthLocation
    """Observer Location"""


def get_ipp_from_skycoord(
    loc: EarthLocation, times: Time, source: SkyCoord, height_array: u.Quantity
) -> IPP:
    """Get ionospheric piercepoints

    Parameters
    ----------
    loc : EarthLocation
        observer location
    times : Time
        observation times
    source : SkyCoord
        source location
    height_array : u.Quantity
        array of altitudes

    Returns
    -------
    IPP
        Ionospheric piercepoints
    """
    # Note: at the moment we calculate ipp per station. I think this is ok,
    # but maybe we need to include a many stations option
    aa = AltAz(location=loc, obstime=times)
    altaz = source.transform_to(aa)
    return get_ipp_from_altaz(loc, altaz, height_array)


def get_ipp_from_altaz(
    loc: EarthLocation, altaz: AltAz, height_array: u.Quantity
) -> IPP:
    """get ionospheric piercepoints from azimuth elevations

    Parameters
    ----------
    loc : EarthLocation
        observer location
    altaz : AltAz
        azimuth and elevations for all times
    height_array : u.Quantity
        array of altitudes

    Returns
    -------
    IPP
        ionospheric piercepoints
    """
    if not altaz.obstime.shape or altaz.obstime.shape != altaz.az.shape:
        altaz = _make_dimensions_match(altaz)
    los_dir = altaz.transform_to(ITRS(obstime=altaz.obstime))
    ipp, airmass = _get_ipp_simple(height_array=height_array, loc=loc, los_dir=los_dir)
    return IPP(
        loc=EarthLocation.from_geocentric(*ipp),
        times=altaz.obstime,
        los=los_dir,
        airmass=airmass,
        altaz=altaz,
        station_loc=loc,
    )


def _make_dimensions_match(altaz: AltAz) -> AltAz:
    """Helper function to change time dimensions suchthat they correspond to the altaz dimension

    Parameters
    ----------
    altaz : AltAz
        the altaz object

    Returns
    -------
    AltAz
        altaz object with matching obstime dimension

    Raises
    ------
    NotImplementedError
        multiple times with different shape than altaz is not implemented yet
    """
    times = altaz.obstime
    az = altaz.az
    # if multiple azimuth/altitudes for one time, just increase dimensions of time
    if not times.shape:
        times = Time(times.mjd * np.ones(az.shape), format="mjd")
    if times.shape != az.shape:
        msg = (
            "Support for multiple times for azimuth/elevation grids is not implemented"
        )
        raise NotImplementedError(msg)

    return AltAz(az=altaz.az, alt=altaz.alt, obstime=times, location=altaz.location)


# TODO: Create return type for this function
def _get_ipp_simple(
    height_array: u.Quantity, loc: EarthLocation, los_dir: SkyCoord
) -> tuple[list[u.Quantity], NDArray[np.float64]]:
    r"""helper function to calculate ionospheric piercepoints using a simple spherical earth model

    .. code-block::

        |loc + alphas * los_dir| = R_earth + height_array

    solve for alphas using abc formula

    Parameters
    ----------
    height_array : u.Quantity
        array of altitudes
    loc : EarthLocation
        observer location
    los_dir : ITRS
        line of sight, unit vector

    Returns
    -------
    tuple(list[u.Quantity], NDArray)
        ipp.x, ipp.y, ipp.z positions, airmass
    """
    logger.info("Calculating ionospheric piercepoints")
    c_value = R_earth**2 - (R_earth + height_array) ** 2
    b_value = u.Quantity(loc.geocentric) @ los_dir.cartesian.xyz
    b_value = b_value[:, np.newaxis]
    alphas = -b_value + np.sqrt(b_value**2 - c_value)
    ipp = (
        u.Quantity(loc.geocentric)[:, np.newaxis, np.newaxis]
        + alphas * los_dir.cartesian.xyz[:, :, np.newaxis]
    )
    inv_airmass = np.einsum("ijk,ij->jk", ipp, los_dir.cartesian.xyz)
    inv_airmass /= R_earth + height_array  # normalized
    airmass = (
        1.0 / inv_airmass.decompose().value
    )  # if you forget the .decompose it can have airmass in (m/km)
    return ipp, airmass
