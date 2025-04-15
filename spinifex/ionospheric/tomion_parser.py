"""Module to parse the UPC-Ionsat tomion data format"""

from __future__ import annotations

import gzip
from datetime import datetime
from pathlib import Path
from typing import Any, NamedTuple

import astropy.units as u
import numpy as np
import requests
from astropy.table import Table
from astropy.time import Time
from numpy.typing import NDArray
from pydantic import Field

from spinifex.exceptions import IonexError, TomionError
from spinifex.geometry import IPP
from spinifex.ionospheric.index_tools import (
    compute_index_and_weights,
    get_interpol,
    get_sorted_indices,
)
from spinifex.logger import logger
from spinifex.options import Options
from spinifex.times import get_indexlist_unique_days, get_unique_days

TOMOION_FORMAT_DICT: dict[str, Any] = {
    "mjd": float,
    "index": int,
    "value": float,
    "stddev": float,
    "type": str,
    "number_of_observations": int,
    "height": float,
    "ra": float,
    "dec": float,
    "i": int,
    "j": int,
    "k": int,
    "label": str,
    "longitude": float,
    "lst": float,
    "year": int,
    "doy": int,
    "month": int,
    "dom": int,
}
MAX_INTERPOL_POINTS: int = 4  # number of points used for lon/lat interpolation
TOMION_HEIGHTS: u.Quantity = [450, 1150] * u.km


class TomionOptions(Options):
    """Options for tomion model"""

    output_directory: Path | None = Field(
        None, description="Output directory for tomion files"
    )


class TomionData(NamedTuple):
    """Object containing all necessary information from Tomion data"""

    lons: NDArray[np.float64]
    """array with longitude values (degrees)"""
    lats: NDArray[np.float64]
    """array with latitude values (degrees)"""
    available_times: Time
    """array with unique times"""
    times: Time
    """times"""
    h: NDArray[np.float64]
    """heights (km)"""
    tec: NDArray[np.float64]
    """array with voxel tecvalues(TECU)"""
    h_idx: NDArray[np.int128]
    """array with index of height"""


def _read_tomion(fname: Path) -> TomionData:
    """reads a tomion format file and returns a TomionData object

    Parameters
    ----------
    fname : Path
        filename

    Returns
    -------
    TomionData
        object with data and axes of the data
    """
    try:
        tomion_data = Table.read(
            fname,
            format="ascii",
            names=list(TOMOION_FORMAT_DICT.keys()),
        )
    except Exception as e:
        msg = f"Could not read tomion file {fname}"
        raise TomionError(msg) from e

    time_column = tomion_data["mjd"].value
    available_times = Time(np.unique(time_column), format="mjd")
    times = Time(time_column, format="mjd")
    height = tomion_data["height"].value
    layer_height = np.ptp(height)
    layer_index = tomion_data["k"].value
    # TODO: What are these numbers?
    electron_density = (10.0 / 1.05) * layer_height * tomion_data["value"].value
    available_lat = tomion_data["dec"].value
    available_lon = tomion_data["longitude"].value
    return TomionData(
        lons=available_lon,
        lats=available_lat,
        available_times=available_times,
        times=times,
        h=height,
        h_idx=layer_index,
        tec=electron_density,
    )


def get_tomion_paths(
    unique_days: Time, tomion_options: TomionOptions | None = None
) -> list[Path]:
    """download tomion data for all unique days

    Parameters
    ----------
    unique_days : Time
        days for which to download data
    tomion_options : TomionOptions | None, optional
        options for the tomion model, by default None

    Returns
    -------
    list[Any]
        list of paths to the files
    """
    tomion_paths: list[Path] = []
    for day in unique_days:
        url, nefull_name = _tomion_format(day)
        msg = f"downloading {url} to {nefull_name}"
        logger.info(msg)
        tomion_paths.append(
            _download_tomion_file(
                url=url, nefull_name=nefull_name, tomion_options=tomion_options
            )
        )
    return tomion_paths


def _tomion_format(time: Time) -> tuple[str, str]:
    """helper function to get the url of the tomion data

    Parameters
    ----------
    time : Time
        day for which to get the url

    Returns
    -------
    tuple[Any, Any]
        url and the name how the data will be stored on disc
    """
    assert time.isscalar, "Only one time is supported"
    url_stem = "http://cabrera.upc.es/upc_ionex_GPSonly-RINEXv3/"
    dtime: datetime = time.to_datetime()
    doy = time.datetime.timetuple().tm_yday
    # YYYY/DDD_YYMMDD.15min/.bias_dens
    nefull_name = f"NeFull.{dtime.year:04d}{doy:03d}"
    yy = f"{dtime.year:02d}"[-2:]
    file_name = f"bias_dens.0001.{dtime.year:04d}{doy:03d}.gz"
    directory_name = f"{dtime.year:04d}/{doy:03d}_{yy}{dtime.month:02d}{dtime.day:02d}.15min/.bias_dens"
    return f"{url_stem}/{directory_name}/{file_name}", nefull_name


def _download_tomion_file(
    url: str,
    nefull_name: str,
    tomion_options: TomionOptions | None = None,
    timeout_seconds: int = 30,
    chunk_size: int = 1000,
) -> Path:
    """download and convert tomion file to readable nefull text file or if the file already exists
     just return a pointer to the file

    Parameters
    ----------
    url : str
        url of the file to download
    nefull_name : str
        name of the extracted nefull file
    tomion_options : TomionOptions | None, optional
        options for the ionospheric model, by default None
    timeout_seconds : int, optional
        time out for downloading, by default 30
    chunk_size : int, optional
        chunksize for downloading, by default 1000

    Returns
    -------
    Path
        pointer to the nefull file

    Raises
    ------
    IonexError
        error if the download times out
    """
    # TODO: Consider merging with download tools
    if tomion_options is None or tomion_options.output_directory is None:
        output_directory = Path.cwd() / "tomion_files"

    output_directory.mkdir(exist_ok=True)
    output_file = output_directory / nefull_name

    if output_file.exists():
        msg = f"File {output_file} already exists. Skipping download."
        logger.info(msg)
        return output_file
    msg = f"Downloading from {url}"
    logger.info(msg)
    try:
        response = requests.get(url, timeout=timeout_seconds)
    except requests.exceptions.Timeout as e:
        msg = "Timed out connecting to server"
        logger.error(msg)
        raise IonexError(msg) from e
    response.raise_for_status()

    tomion_file = output_directory / url.split("/")[-1]
    with tomion_file.open("wb") as file_desc:
        for chunk in response.iter_content(chunk_size=chunk_size):
            file_desc.write(chunk)
    _extract_nefull(tomion_file, output_file)
    return output_file


def _extract_nefull(
    tomion_file: Path, output_file: Path, search_term: str = "NeFull"
) -> Path:
    """helper function to get the relevant information from the large tomion file.
    Deletes the tomion file if the extraction was successful

    Parameters
    ----------
    tomion_file : Path
        pointer to the tomion file
    output_file : Path
        name of the nefull file
    search_term : str, optional
        the indicator for useful data in the tomion file, by default "NeFull"

    Returns
    -------
    Path
        pointer to the nefull file
    """
    with gzip.open(tomion_file, "rt") as f_in, Path.open(output_file, "w") as f_out:
        for line in f_in:
            if search_term in line:
                f_out.write(line)
    # Verify data was written correctly
    if Path.exists(output_file) and Path(output_file).stat().st_size > 0:
        msg = f"Extraction successful! Output saved in {output_file}"
        logger.info(msg)
        tomion_file.unlink()
        return output_file

    msg = f"Could not convert {tomion_file} to {output_file}"
    raise TomionError(msg)


def interpolate_tomion(
    tomion: TomionData,
    lons: NDArray[np.float64],
    lats: NDArray[np.float64],
    times: Time,
) -> NDArray[np.float64]:
    """Interpolate tomion data to the requested lons/lats/times

    Parameters
    ----------
    tomion : TomionData
        data object
    lons : NDArray[np.float64]
        array of longitudes at the two TOMION_HEIGHTS, shape (2,)
    lats : NDArray[np.float64]
        array of latitudes at the two TOMION_HEIGHTS, shape (2,)
    times : Time
        time

    Returns
    -------
    NDArray[np.float64]
        electron density values at two TOMION_HEIGHTS, shape (2,)
    """
    # TODO: implement apply_earth_rotation
    # TODO: implement this function directly for an array of times
    timeindex = compute_index_and_weights(tomion.available_times.mjd, times.mjd)
    time1 = tomion.available_times.mjd[timeindex.idx1]
    time2 = tomion.available_times.mjd[timeindex.idx2]
    timeselect1 = tomion.times.mjd == time1
    timeselect2 = tomion.times.mjd == time2
    timeselect = [timeselect1, timeselect2]
    # get data for two layers for two times
    layers_lo = [tomion.h_idx[timeselect1] == 1, tomion.h_idx[timeselect2] == 1]
    layers_hi = [tomion.h_idx[timeselect1] == 2, tomion.h_idx[timeselect2] == 2]
    # get lon,lat idx for these
    tec_lo = []
    tec_hi = []
    tec = np.zeros((2,), dtype=float)
    for lo, tms in zip(layers_lo, timeselect):
        isorted_low = get_sorted_indices(
            lon=lons[0],
            lat=lats[0],
            avail_lon=tomion.lons[tms][lo],
            avail_lat=tomion.lats[tms][lo],
        )
        tec_lo.append(
            get_interpol(
                tomion.tec[tms][lo][isorted_low.indices[:MAX_INTERPOL_POINTS]],
                isorted_low.distance[:MAX_INTERPOL_POINTS],
            )
        )
    for hi, tms in zip(layers_hi, timeselect):
        isorted_hi = get_sorted_indices(
            lon=lons[1],
            lat=lats[1],
            avail_lon=tomion.lons[tms][hi],
            avail_lat=tomion.lats[tms][hi],
        )
        tec_hi.append(
            get_interpol(
                tomion.tec[tms][hi][isorted_hi.indices[:MAX_INTERPOL_POINTS]],
                isorted_hi.distance[:MAX_INTERPOL_POINTS],
            )
        )

    tec[0] = tec_lo[0] * timeindex.w1[0] + tec_lo[1] * timeindex.w2[0]
    tec[1] = tec_hi[0] * timeindex.w1[0] + tec_hi[1] * timeindex.w2[0]

    return tec


def get_density_dual_layer(
    ipp: IPP, tomion_options: TomionOptions | None = None
) -> NDArray[np.float64]:
    """extracts electron densities for the two TOMION_HEIGHTS for all times in ipp. The returned array
    will have zeros every where apart from the two altitudes closest to TOMION_HEIGHTS

    Parameters
    ----------
    ipp : IPP
        input piercepoint locations
    tomion_options : TomionOptions | None, optional
        optional ionospheric model options, by default None

    Returns
    -------
    NDArray[np.float64]
        array with electron densities shape (ipp.loc.shape). The array will only have values at the
        altitudes closest to the TOMION_HEIGHTS and zeros everywhere else.

    Raises
    ------
    FileNotFoundError
        error if the tomion files are not available locally nor online
    """
    # TODO: no need to go through all the burden, just make sure that the ipps are correct for this model?
    h_index = [
        np.argmin(np.abs(ipp.loc[0].height.to(u.km).value - h.to(u.km).value))
        for h in TOMION_HEIGHTS
    ]
    selected_ipp = ipp.loc[:, h_index]
    tec = np.zeros(ipp.loc.shape, dtype=float)
    unique_days = get_unique_days(times=ipp.times)
    if not unique_days.shape:
        url, nefull_name = _tomion_format(time=unique_days)
        tomion_path = _download_tomion_file(
            url=url, nefull_name=nefull_name, tomion_options=tomion_options
        )
        if not tomion_path.exists():
            msg = f"Tomion file {tomion_path} not found!"
            raise FileNotFoundError(msg)
        tomion = _read_tomion(tomion_path)
        for ippi in range(tec.shape[0]):
            tec[ippi, h_index] = interpolate_tomion(
                tomion,
                selected_ipp[ippi].lon.deg,
                selected_ipp[ippi].lat.deg,
                ipp.times[ippi],
            )

        return tec
    group_indices = get_indexlist_unique_days(unique_days, ipp.times)
    sorted_tomion_paths = get_tomion_paths(unique_days)
    for indices, tomion_file in zip(group_indices, sorted_tomion_paths):
        if not tomion_file.exists():
            msg = f"Tomion file {tomion_file} not found!"
            raise FileNotFoundError(msg)
        u_loc = selected_ipp[indices]
        u_times = ipp.times[indices]
        tomion = _read_tomion(tomion_file)
        for idxi, ippi in enumerate(np.arange(tec.shape[0])[indices]):
            tec[ippi, h_index] = interpolate_tomion(
                tomion, u_loc[idxi].lon.deg, u_loc[idxi].lat.deg, u_times[idxi]
            )

    return tec
