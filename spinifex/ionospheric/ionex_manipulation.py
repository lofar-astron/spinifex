"""Module of ionex manipulation tools"""

from __future__ import annotations

from typing import NamedTuple

import astropy.units as u
import numpy as np
from astropy.time import Time
from numpy.typing import NDArray

from spinifex.geometry.get_ipp import IPP
from spinifex.ionospheric.index_tools import (
    compute_index_and_weights,
    get_indices_axis,
)
from spinifex.ionospheric.ionex_download import (
    _download_ionex,
)
from spinifex.ionospheric.ionex_parser import (
    IonexData,
    read_ionex,
    unique_days_from_ionex_files,
)
from spinifex.ionospheric.tec_data import ElectronDensity, IonexOptions
from spinifex.times import get_indexlist_unique_days


class GroupedIPPs(NamedTuple):
    """Grouped IPPs"""

    ipps: list[IPP]
    indices: list[NDArray[np.int32]]


def interpolate_ionex(
    ionex: IonexData,
    lons: NDArray[np.float64],
    lats: NDArray[np.float64],
    times: Time,
    apply_earth_rotation: float = 1,
    get_rms: bool = False,
) -> NDArray[np.float64]:
    """Interpolate ionex data to a given lon/lat/height grid.
    lons, lats, times all should have the same length

    apply_earth_rotation:
    This is assuming that the TEC maps move according to the rotation
    Earth (following method 3 of interpolation described in the IONEX
    document). Experiments with high time resolution ROB data show that
    this is not really the case, resulting in strange wavelike structures
    when applying this smart interpolation.
    TODO: implement smoothing filters for interpolation

    Parameters
    ----------
    ionex : IonexData
        ionex object containing the information of the ionex file
    lons : NDArray
        longitudes (deg) of all points to interpolate to
    lats : NDArray
        lattitudes (deg) of all points to interpolate to
    times : Time
        times of all points to interpolate to
    apply_earth_rotation : float, optional
        specify (with a number between 0 and 1) how much of the earth rotation
        is taken in to account in the interpolation step., by default 1
    get_rms : bool, optional
        use rms values instead of tec values

    Returns
    -------
    NDArray
        array with interpolated tec values
    """
    value_array = ionex.rms if get_rms else ionex.tec

    timeindex = compute_index_and_weights(ionex.times.mjd, times.mjd)
    latindex = compute_index_and_weights(ionex.lats, lats)
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
        value_array[timeindex.idx1, lonindex1.idx1, latindex.idx1]
        * lonindex1.w1
        * timeindex.w1
        * latindex.w1
    )
    tecdata += (
        value_array[timeindex.idx1, lonindex1.idx2, latindex.idx1]
        * lonindex1.w2
        * timeindex.w1
        * latindex.w1
    )

    tecdata += (
        value_array[timeindex.idx2, lonindex2.idx1, latindex.idx1]
        * lonindex2.w1
        * timeindex.w2
        * latindex.w1
    )
    tecdata += (
        value_array[timeindex.idx2, lonindex2.idx2, latindex.idx1]
        * lonindex2.w2
        * timeindex.w2
        * latindex.w1
    )

    tecdata += (
        value_array[timeindex.idx1, lonindex1.idx1, latindex.idx2]
        * lonindex1.w1
        * timeindex.w1
        * latindex.w2
    )
    tecdata += (
        value_array[timeindex.idx1, lonindex1.idx2, latindex.idx2]
        * lonindex1.w2
        * timeindex.w1
        * latindex.w2
    )

    tecdata += (
        value_array[timeindex.idx2, lonindex2.idx1, latindex.idx2]
        * lonindex2.w1
        * timeindex.w2
        * latindex.w2
    )
    tecdata += (
        value_array[timeindex.idx2, lonindex2.idx2, latindex.idx2]
        * lonindex2.w2
        * timeindex.w2
        * latindex.w2
    )
    return tecdata


def get_density_ionex(
    ipp: IPP, ionex_options: IonexOptions | None = None
) -> ElectronDensity:
    """read ionex files and interpolate values to ipp locations/times

    Parameters
    ----------
    ipp : IPP
        ionospheric piercepoints
    ionex_options : IonexOptions, optional
        optional arguments for the ionospheric model, by default None

    Returns
    -------
    ElectronDensity
        object with arrays of tec  and tec_rms values for every entry in ipp

    Raises
    ------
    FileNotFoundError
        if ionex file cannot be downloaded
    """
    if ionex_options is None:
        ionex_options = IonexOptions()
    # TODO: apply_earth_rotation as option
    sorted_ionex_paths = _download_ionex(times=ipp.times, options=ionex_options)
    # also download data for next day, to remove midnight jumps
    sorted_next_day_paths = (
        _download_ionex(times=ipp.times + 1 * u.day, options=ionex_options)
        if ionex_options.remove_midnight_jumps
        else [
            None,  # type: ignore[list-item]
        ]
        * len(sorted_ionex_paths)
    )
    unique_days = unique_days_from_ionex_files(sorted_ionex_paths)
    if not unique_days.shape:
        ionex = read_ionex(
            sorted_ionex_paths[0], sorted_next_day_paths[0], options=ionex_options
        )
        tec = interpolate_ionex(ionex, ipp.loc.lon.deg, ipp.loc.lat.deg, ipp.times)
        electron_density_error = interpolate_ionex(
            ionex, ipp.loc.lon.deg, ipp.loc.lat.deg, ipp.times, get_rms=True
        )
        return ElectronDensity(
            electron_density=tec, electron_density_error=electron_density_error
        )
    group_indices = get_indexlist_unique_days(unique_days, ipp.times)
    tec = np.zeros(ipp.loc.shape, dtype=float)
    electron_density_error = np.full(ipp.loc.shape, np.nan)
    for indices, ionex_file, next_day_file in zip(
        group_indices, sorted_ionex_paths, sorted_next_day_paths
    ):
        if not ionex_file.exists():
            msg = f"Ionex file {ionex_file} not found!"
            raise FileNotFoundError(msg)
        u_loc = ipp.loc[indices]
        u_times = ipp.times[indices]
        ionex = read_ionex(ionex_file, next_day_file, options=ionex_options)
        tec[indices] = interpolate_ionex(ionex, u_loc.lon.deg, u_loc.lat.deg, u_times)
        electron_density_error[indices] = interpolate_ionex(
            ionex, u_loc.lon.deg, u_loc.lat.deg, u_times, get_rms=True
        )

    return ElectronDensity(
        electron_density=tec, electron_density_error=electron_density_error
    )
