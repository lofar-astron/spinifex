"""Module of ionex manipulation tools"""

from __future__ import annotations

from importlib import resources
from pathlib import Path

import numpy as np
from astropy.time import Time
from numpy.typing import ArrayLike

from spinifex.geometry.get_ipp import IPP
from spinifex.ionospheric.index_tools import (
    _compute_index_and_weights,
    get_indices_axis,
)
from spinifex.ionospheric.ionex_download import SortedIonexPaths, download_ionex
from spinifex.ionospheric.ionex_parser import IonexData, read_ionex


def interpolate_ionex(
    ionex: IonexData,
    lons: ArrayLike,
    lats: ArrayLike,
    times: Time,
    apply_earth_rotation: float = 1,
) -> ArrayLike:
    """Interpolate ionex data to a given lon/lat/height grid.
    lons, lats, times all should have the same length
    apply_earth_rotation:
    This is assuming that the TEC maps move according to the rotation
    Earth (following method 3 of interpolation described in the IONEX
    document). Experiments with high time resolution ROB data show that
    this is not really the case, resulting in strange wavelike structures
    when applying this smart interpolation.
    Parameters
    ----------
    ionex : IonexData
        ionex object containing the information of the ionex file
    lons : ArrayLike
        longitudes (deg) of all points to interpolate to
    lats : ArrayLike
        lattitudes (deg) of all points to interpolate to
    times : Time
        times of all points to interpolate to
    apply_earth_rotation : float, optional
        specify (with a number between 0 and 1) how much of the earth rotation
        is taken in to account in the interpolation step., by default 1

    Returns
    -------
    ArrayLike
        array with interpolated tec values
    """
    timeindex = _compute_index_and_weights(ionex.times.mjd, times.mjd)
    latindex = _compute_index_and_weights(ionex.lats, lats)
    # take into account earth rotation
    if apply_earth_rotation > 0:
        rot1 = (
            (times.mjd - ionex.times.mjd[timeindex.idx1]) * 360.0
        ) * apply_earth_rotation
        rot2 = (
            (times.mjd - ionex.times.mjd[timeindex.idx2]) * 360.0
        ) * apply_earth_rotation
        lonindex1 = get_indices_axis(lons + rot1, ionex.lons, wrap_unit=360)
        lonindex2 = get_indices_axis(lons + rot2, ionex.lons, wrap_unit=360)
    else:
        lonindex1 = get_indices_axis(lons, ionex.lons, wrap_unit=360)
        lonindex2 = lonindex1
    tecdata = (
        ionex.tec[timeindex.idx1, lonindex1.idx1, latindex.idx1]
        * lonindex1.w1
        * timeindex.w1
        * latindex.w1
    )
    tecdata += (
        ionex.tec[timeindex.idx1, lonindex1.idx2, latindex.idx1]
        * lonindex1.w2
        * timeindex.w1
        * latindex.w1
    )

    tecdata += (
        ionex.tec[timeindex.idx2, lonindex2.idx1, latindex.idx1]
        * lonindex2.w1
        * timeindex.w2
        * latindex.w1
    )
    tecdata += (
        ionex.tec[timeindex.idx2, lonindex2.idx2, latindex.idx1]
        * lonindex2.w2
        * timeindex.w2
        * latindex.w1
    )

    tecdata += (
        ionex.tec[timeindex.idx1, lonindex1.idx1, latindex.idx2]
        * lonindex1.w1
        * timeindex.w1
        * latindex.w2
    )
    tecdata += (
        ionex.tec[timeindex.idx1, lonindex1.idx2, latindex.idx2]
        * lonindex1.w2
        * timeindex.w1
        * latindex.w2
    )

    tecdata += (
        ionex.tec[timeindex.idx2, lonindex2.idx1, latindex.idx2]
        * lonindex2.w1
        * timeindex.w2
        * latindex.w2
    )
    tecdata += (
        ionex.tec[timeindex.idx2, lonindex2.idx2, latindex.idx2]
        * lonindex2.w2
        * timeindex.w2
        * latindex.w2
    )
    return tecdata


def get_ionex_file(
    times: Time, server: str | None = None, prefix: str = "CODG"
) -> Path | None:
    """Find ionex file locally or online.

    Parameters
    ----------
    times : Time
        time for which to retrieve
    server : str, optional
        URL of a server to query. If not specified, local files will be searched.
    prefix : str, optional
        name of the ionex provider, by default "CODG"

    Returns
    -------
    str
        _description_
    """

    if server is None:
        doy = times[0].datetime.timetuple().tm_yday
        yy = times[0].ymdhms[0]
        yy = yy - 2000 if yy > 2000 else yy - 1990
        with resources.as_file(resources.files("spinifex.data.tests")) as test_data:
            return SortedIonexPaths(
                ionex_list=[test_data / f"{prefix}{doy:03d}0.{yy:02d}I.gz"],
                unique_days=times[0],
            )
    else:
        return download_ionex(times=times, server=server, prefix=prefix)


def _group_by_day(ipp: IPP, unique_days: Time) -> list[IPP]:
    if not unique_days.shape:
        return ipp
    ipps = []
    for day in unique_days.mjd:
        indices = np.floor(ipp.times.mjd) == day
        ipps.append(
            IPP(
                times=ipp.times[indices],
                loc=ipp.loc[indices],
                los=ipp.los[indices],
                airmass=ipp.airmass[indices],
                altaz=ipp.altaz[indices],
            )
        )
    return ipps


def _read_ionex_stuff(ipp: IPP, server: str | None = None) -> ArrayLike:
    sorted_ionex_paths = get_ionex_file(ipp.times, server=server)
    if not sorted_ionex_paths.unique_days.shape:
        ionex = read_ionex(sorted_ionex_paths.ionex_list[0])
        return interpolate_ionex(ionex, ipp.loc.lon.deg, ipp.loc.lat.deg, ipp.times)
    day_groups = _group_by_day(ipp, sorted_ionex_paths.unique_days)
    tec = []
    for u_ipp, ionex_file in zip(day_groups, sorted_ionex_paths.ionex_list):
        if ionex_file is None:
            msg = "No ionex file found!"
            raise FileNotFoundError(msg)
        ionex = read_ionex(ionex_file)
        tec.append(
            interpolate_ionex(ionex, u_ipp.loc.lon.deg, u_ipp.loc.lat.deg, u_ipp.times)
        )
    return np.concatenate(tec)  # wrong if order has changed
