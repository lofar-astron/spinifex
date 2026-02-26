"""Module for getting the Ionospheric Piercepoints"""

from __future__ import annotations

from typing import NamedTuple

import astropy.units as u
import numpy as np
from astropy.coordinates import ITRS, AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from numpy.typing import NDArray

from spinifex.logger import logger

R_EARTH_MEAN = 6371.0 * u.km
R_EARTH_MEAN_KM = 6371.0


class IPP(NamedTuple):
    """Ionospheric Piercepoints"""

    lon: u.Quantity
    """longitude on a spherical Earth approximation"""
    lat: u.Quantity
    """latitude on a spherical Earth approximation"""
    height: u.Quantity
    """geocentric height above Earth center"""
    times: Time
    """array of times"""
    los: NDArray[np.float64]
    """Line of sight direction in ITRS coordinates, normalized"""
    airmass: NDArray[np.float64]
    """airmass factor to convert to slant values"""
    altaz: AltAz
    """azimuth elevation"""
    station_loc: EarthLocation
    """Observer Location"""


def get_ipp_from_itrs(
    loc: EarthLocation, times: Time, los_dir: ITRS, height_array: u.Quantity
) -> IPP:
    if not times.shape:
        times = Time(np.array([times.mjd]), format="mjd")
    los_vector = los_dir.cartesian.xyz.value
    los_vector /= np.linalg.norm(los_vector, axis=0)
    ipp, airmass = _get_ipp_simple(
        height_array=height_array, loc=loc, los_dir=los_vector
    )
    lon, lat, radius = _xyz_to_lonlat_spherical(*[i.to(u.m).value for i in ipp])
    # ipploc = EarthLocation.from_geodetic(lon, lat, height_array)  # TOO slow
    altaz = los_dir.transform_to(AltAz(obstime=times, location=loc))
    return IPP(
        lon=lon,
        lat=lat,
        height=radius,
        times=times,
        los=los_vector.T,
        airmass=airmass,
        altaz=altaz,
        station_loc=loc,
    )


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
    if not times.shape:
        times = Time(np.array([times.mjd]), format="mjd")

    source_itrs = source.transform_to(ITRS(obstime=times, location=loc))
    return get_ipp_from_itrs(loc, times, source_itrs, height_array)


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
    los_dir = altaz.transform_to(ITRS(obstime=altaz.obstime, location=altaz.location))
    return get_ipp_from_itrs(loc, altaz.obstime, los_dir, height_array)


def _xyz_to_lonlat_spherical(
    x_m: float, y_m: float, z_m: float
) -> tuple[u.Quantity, u.Quantity, u.Quantity]:
    """Convert geocentric XYZ to geodetic lon/lat/radius using spherical model."""
    lon_rad = np.arctan2(y_m, x_m)
    r_m = np.sqrt(x_m**2 + y_m**2 + z_m**2)
    lat_rad = np.arcsin(np.clip(z_m / r_m, -1, 1))

    return (
        np.rad2deg(lon_rad) * u.deg,
        np.rad2deg(lat_rad) * u.deg,
        r_m * u.m,
    )  # ← Geocentric radius!


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
def _get_ipp_simple_SLOW(
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
    r_station = u.Quantity(loc.geocentric)

    c_value = np.sum(r_station**2) - (R_EARTH_MEAN + height_array) ** 2
    if len(los_dir.shape) == 1:
        los_dir = los_dir[:, np.newaxis]  # make sure b_values is an array
    b_value = r_station @ los_dir
    b_value = b_value[:, np.newaxis]
    alphas = -b_value + np.sqrt(b_value**2 - c_value)
    ipp = (
        u.Quantity(loc.geocentric)[:, np.newaxis, np.newaxis]
        + alphas * los_dir[:, :, np.newaxis]
    )
    inv_airmass = np.einsum("ijk,ij->jk", ipp, los_dir)
    inv_airmass /= R_EARTH_MEAN + height_array  # normalized
    airmass = (
        1.0 / inv_airmass.decompose().value
    )  # if you forget the .decompose it can have airmass in (m/km)
    return ipp, airmass


def _get_ipp_simple_FAST(
    height_array: u.Quantity, loc: EarthLocation, los_dir: NDArray[np.float64]
) -> tuple[list[u.Quantity], NDArray[np.float64]]:
    """
    Optimized version: Convert to numpy arrays immediately, compute in km.

    Parameters
    ----------
    height_array : u.Quantity
        Array of altitudes (will be converted to km)
    loc : EarthLocation
        Observer location
    los_dir : NDArray
        Line of sight unit vectors (3, n_times) or (3,)

    Returns
    -------
    tuple
        ipp positions (as list of Quantities in meters), airmass (dimensionless)
    """
    logger.info("Calculating ionospheric piercepoints")

    # Convert everything to numpy arrays in km (do this ONCE)
    height_km = height_array.to(u.km).value  # Shape: (n_heights,)

    # Get station position in km
    loc_xyz_km = np.array(
        [loc.x.to(u.km).value, loc.y.to(u.km).value, loc.z.to(u.km).value]
    )  # meters

    # Get LOS direction (already unitless)
    if los_dir.ndim == 1:
        los_dir = los_dir[:, np.newaxis]  # Shape: (3, n_times)

    # All calculations in km (pure numpy - FAST!)
    loc_radius_squared = np.sum(loc_xyz_km**2)  # scalar
    target_radius = R_EARTH_MEAN_KM + height_km  # Shape: (n_heights,)
    c_value = loc_radius_squared - target_radius**2  # Shape: (n_heights,)

    # Dot product: loc · los (pure numpy)
    b_value = loc_xyz_km @ los_dir  # Shape: (n_times,)
    b_value = b_value[:, np.newaxis]  # Shape: (n_times, 1)

    # Solve quadratic (pure numpy - FAST!)
    # alpha = -b + sqrt(b^2 - c)
    # Broadcasting: (n_times, 1) and (n_heights,) → (n_times, n_heights)
    alphas = -b_value + np.sqrt(b_value**2 - c_value)  # Shape: (n_times, n_heights)

    # Calculate pierce points (pure numpy - FAST!)
    # ipp = loc + alpha * los
    # Broadcasting: (3, 1, 1) + (n_times, n_heights) * (3, n_times, 1)
    ipp_km = (
        loc_xyz_km[:, np.newaxis, np.newaxis]  # (3, 1, 1)
        + alphas[np.newaxis, :, :]
        * los_dir[:, :, np.newaxis]  # (3, n_times, n_heights)
    )

    # Calculate airmass (pure numpy - FAST!)
    # inv_airmass = ipp · los / |ipp|
    inv_airmass = np.einsum("ijk,ij->jk", ipp_km, los_dir)  # (n_times, n_heights)
    inv_airmass /= target_radius  # Normalize by target radius
    airmass = 1.0 / inv_airmass  # (n_times, n_heights)

    # Convert back to meters with units (do this ONCE at the end)
    ipp_m = ipp_km * 1000.0  # Convert km back to meters
    ipp_with_units = u.Quantity(ipp_m, unit=u.m)
    # [x, y, z] each with shape (n_times, n_heights)

    return ipp_with_units, airmass


_get_ipp_simple = _get_ipp_simple_FAST
