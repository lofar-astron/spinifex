#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module to parse the IONosphere map EXchange (IONEX) data format,
as described in Schaer and Gurtner (1998)"""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple, TextIO

import astropy.units as u
import numpy as np
from astropy.time import Time


class IonexData(NamedTuple):
    """Object containing all necessary information from Ionex data"""

    lons: np.ndarray[float]
    """array with available longitude values (degrees)"""
    lats: np.ndarray[float]
    """array with available latitude values (degrees)"""
    times: Time
    """available times"""
    dims: int
    """dimension of the heights (usually 1)"""
    h: np.ndarray[float]
    """available heights (km)"""
    tec: np.ndarray[float]
    """array with tecvalues times x lons x lats (TECU)"""
    rms: np.ndarray[float]
    """array with rms of tecvalues times x lons x lats (TECU, if available, zeros otherwise)"""


class IonexHeader(NamedTuple):
    """Object containing header information from ionex file"""

    lons: np.ndarray[float]
    """array with available longitude values (degrees)"""
    lats: np.ndarray[float]
    """array with available latitude values (degrees)"""
    times: Time
    """available times"""
    dims: int
    """dimension of the heights (usually 1)"""
    h: np.ndarray[float]
    """available heights (km)"""
    mfactor: float
    """multiplication factor for tec values"""


def read_ionex(ionex_filename: Path) -> IonexData:
    """Read and parse a ionex file. Returns a ionex object.

    Parameters
    ----------
    ionex_filename : str
        _description_

    Returns
    -------
    IonexData
        ionex object with data and grid

    """
    with Path.open(ionex_filename, encoding="utf-8") as myf:
        return _read_ionex_data(myf)


def _read_ionex_header(filep: TextIO) -> IonexHeader:
    """Read header from ionex file. Put filepointer to the end of the header."""
    filep.seek(0)
    h1, h2, hstep = (None,) * 3
    start_lon, end_lon, step_lon, start_lat, end_lat, step_lat = (None,) * 6
    start_time, ntimes, step_time = (None,) * 3
    mfactor, dimension = (None,) * 2
    for line in filep:
        if "END OF HEADER" in line:
            break
        label = line[60:-1]
        record = line[:60]
        if "EPOCH OF FIRST MAP" in label:
            yy, mm, day, hr, minute, second = (int(i) for i in record.strip().split())
            epoch = Time(f"{yy}-{mm}-{day}T{hr%24}:{minute}:{second}")
            start_time = epoch
        if "INTERVAL" in label:
            step_time = float(record) * u.s
        if "EXPONENT" in label:
            mfactor = 10.0 ** float(record)
        if "MAP DIMENSION" in label:
            dimension = int(record)
        if "HGT1 / HGT2 / DHGT" in label:
            h1, h2, hstep = (float(i) for i in record.split())
        if "LON1 / LON2 / DLON" in label:
            start_lon, end_lon, step_lon = (float(i) for i in record.split())
        if "LAT1 / LAT2 / DLAT" in label:
            start_lat, end_lat, step_lat = (float(i) for i in record.split())
        if "# OF MAPS IN FILE" in label:
            ntimes = int(record)
    if h1 is None or start_lon is None or start_lat is None or start_time is None:
        msg = f"Not a valid IONex file: {filep.name}"
        raise OSError(msg)

    harray = np.arange(h1, h2 + 0.5 * hstep, hstep) if hstep > 0 else np.array([h1])

    lonarray = np.arange(start_lon, end_lon + 0.5 * step_lon, step_lon)
    latarray = np.arange(start_lat, end_lat + 0.5 * step_lat, step_lat)
    timearray = start_time + np.arange(0, ntimes) * step_time

    return IonexHeader(
        mfactor=mfactor,
        lons=lonarray,
        lats=latarray,
        times=timearray,
        dims=dimension,
        h=harray,
    )


def _fill_data_record(
    data: np.ndarray,
    filep: TextIO,
    stop_label: str,
    timeidx: int,
    ionex_header: IonexHeader,
):
    """Helper function to parse a data block of a single map in ionex.
    Puts filepointer to the end of the map

    Parameters
    ----------
    data : np.ndarray
        pre allocated array to store the datablock
    filep : TextIO
        _description_
    stop_label : str
        end of the data block indicator
    timeidx : int
        index of time of the data block
    ionex_header : namedtuple
        header information
    """
    line = filep.readline()  # read EPOCH (not needed since we have the index)
    tec = []
    lonidx = 0
    latidx = 0
    for line in filep:
        label = line[60:-1]
        if stop_label in label:
            if tec:
                tec = np.array(tec) * ionex_header.mfactor
                data[timeidx, lonidx:, latidx] = tec
            return
        if "LAT/LON1/LON2/DLON/H" in label:
            if tec:
                tec = np.array(tec) * ionex_header.mfactor
                data[timeidx, lonidx:, latidx] = tec
            tec = []
            record = line[:60]
            lat, lon1, _, _, _ = (float(record[i : i + 6]) for i in range(2, 32, 6))
            latidx = np.argmin(np.abs(ionex_header.lats - lat))
            lonidx = np.argmin(np.abs(ionex_header.lons - lon1))
        else:
            record = line[:-1]
            tec += [float(record[i : i + 5]) for i in range(0, len(record), 5)]


def _read_ionex_data(filep: TextIO) -> IonexData:
    """This function parses the IONEX file.
    Some fixed structure (like data records being strings of exactly 80 characters) of the file
    is assumed. This structure is described in Schaer and Gurtner (1998).

    Parameters
    ----------
    filep : TextIO
        pointer to an ionex file

    Returns
    -------
    IonexData
        ionex object
    """
    ionex_header = _read_ionex_header(filep)
    tecarray = np.zeros(
        ionex_header.times.shape + ionex_header.lons.shape + ionex_header.lats.shape,
        dtype=float,
    )
    rmsarray = np.zeros_like(tecarray)
    for line in filep:
        # _read_ionex_header should have put the filep at the end of the header
        label = line[60:-1]
        record = line[:60]
        if "START OF TEC MAP" in label:
            timeidx = int(record) - 1
            _fill_data_record(tecarray, filep, "END OF TEC MAP", timeidx, ionex_header)
        if "START OF RMS MAP" in label:
            timeidx = int(record) - 1
            _fill_data_record(rmsarray, filep, "END OF RMS MAP", timeidx, ionex_header)

    return IonexData(
        lons=ionex_header.lons,
        lats=ionex_header.lats,
        times=ionex_header.times,
        dims=ionex_header.dims,
        h=ionex_header.dims,
        tec=tecarray,
        rms=rmsarray,
    )
