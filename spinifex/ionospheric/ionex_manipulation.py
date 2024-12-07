"""Module of ionex manipulation tools"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.coordinates import EarthLocation
from astropy.time import Time

from spinifex.ionospheric.index_tools import (
    _compute_index_and_weights,
    get_indices_axis,
)
from spinifex.ionospheric.ionex_parser import IonexData, read_ionex


def interpolate_ionex(
    ionex: IonexData,
    lons: np.ndarray,
    lats: np.ndarray,
    times: Time,
    apply_earth_rotation: float = 1,
) -> np.ndarray:
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
    lons : np.ndarray
        longitudes (deg) of all points to interpolate to
    lats : np.ndarray
        lattitudes (deg) of all points to interpolate to
    times : Time
        times of all points to interpolate to
    apply_earth_rotation : float, optional
        specify (with a number between 0 and 1) how much of the earth rotation
        is taken in to account in the interpolation step., by default 1

    Returns
    -------
    np.ndarray
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


def get_ionex_file(time: Time, server: str | None = None, prefix: str = "CODG") -> str:
    """Find ionex file locally or online.

    Parameters
    ----------
    time : Time
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
    doy = time.datetime.timetuple().tm_yday
    yy = time.ymdhms[0]
    yy = yy - 2000 if yy > 2000 else yy - 1990
    if server is None:
        datapath = Path(__file__).parent.parent.parent / "tests" / "data"
        return datapath / f"{prefix}{doy:03d}0.{yy:02d}I"
    return None


def _read_ionex_stuff(loc: EarthLocation, times: Time, **kwargs) -> np.ndarray:
    server = kwargs.get("server")
    ionex_file = get_ionex_file(times[0], server=server)
    ionex = read_ionex(ionex_file)
    return interpolate_ionex(ionex, loc.lon.deg, loc.lat.deg, times)
