"""Module to parse the IONosphere map EXchange (IONEX) data format,
as described in Schaer and Gurtner (1998)"""

from __future__ import annotations

import gzip
import io
from pathlib import Path
from typing import Any, NamedTuple, TextIO

import astropy.units as u
import numpy as np
from astropy.time import Time
from numpy.typing import NDArray
from unlzw3 import unlzw

from spinifex.exceptions import IonexError
from spinifex.ionospheric.tec_data import IonexOptions
from spinifex.times import get_unique_days


class IonexData(NamedTuple):
    """Object containing all necessary information from Ionex data"""

    lons: NDArray[np.float64]
    """array with available longitude values (degrees)"""
    lats: NDArray[np.float64]
    """array with available latitude values (degrees)"""
    times: Time
    """available times"""
    dims: int
    """dimension of the heights (usually 1)"""
    h: NDArray[np.float64]
    """available heights (km)"""
    tec: NDArray[np.float64]
    """array with tecvalues times x lons x lats (TECU)"""
    rms: NDArray[np.float64]
    """array with rms of tecvalues times x lons x lats (TECU, if available, nan otherwise)"""


class IonexHeader(NamedTuple):
    """Object containing header information from ionex file"""

    lons: NDArray[np.float64]
    """array with available longitude values (degrees)"""
    lats: NDArray[np.float64]
    """array with available latitude values (degrees)"""
    times: Time
    """available times"""
    dims: int
    """dimension of the heights (usually 1)"""
    h: NDArray[np.float64]
    """available heights (km)"""
    mfactor: float
    """multiplication factor for tec values"""


def read_ionex(
    ionex_filename: Path,
    next_day_ionex_filename: Path | None = None,
    options: IonexOptions | None = None,
) -> IonexData:
    """Read and parse a ionex file. Returns a ionex object.

    Parameters
    ----------
    ionex_filename : str
        name of the ionex data file
    next_day_ionex_filename : str| None
        name of the ionex data file for the next day (to remove midnight jumps)
    options : IonexOptions
        Options for the ionex model.
    Returns
    -------
    IonexData
        ionex object with data and grid

    """
    if ionex_filename.suffix == ".gz":
        with gzip.open(ionex_filename, "rt", encoding="utf-8") as file_buffer:
            ionex = _read_ionex_data(file_buffer, options=options)

    elif ionex_filename.suffix == ".Z":
        data = unlzw(ionex_filename.read_bytes()).decode("utf-8")
        with io.StringIO(data) as file_buffer:
            ionex = _read_ionex_data(file_buffer, options=options)
    else:
        with ionex_filename.open(encoding="utf-8") as file_buffer:
            ionex = _read_ionex_data(file_buffer, options=options)
    if next_day_ionex_filename is not None:
        if next_day_ionex_filename.suffix == ".gz":
            with gzip.open(
                next_day_ionex_filename, "rt", encoding="utf-8"
            ) as file_buffer:
                next_day_ionex = _read_ionex_data(file_buffer, options=options)

        elif next_day_ionex_filename.suffix == ".Z":
            data = unlzw(next_day_ionex_filename.read_bytes()).decode("utf-8")
            with io.StringIO(data) as file_buffer:
                next_day_ionex = _read_ionex_data(file_buffer, options=options)
        else:
            with next_day_ionex_filename.open(encoding="utf-8") as file_buffer:
                next_day_ionex = _read_ionex_data(file_buffer, options=options)
        ionex = _replace_midnight_data(ionex, next_day_ionex)

    return ionex


def _replace_midnight_data(ionex: IonexData, next_day_ionex: IonexData) -> IonexData:
    """mtigate jumps in tec value at midnight by inserting the tec value of the next day

    Parameters
    ----------
    ionex : IonexData
        ionex data object
    next_day_ionex : IonexData
        ionex data of the next day

    Returns
    -------
    IonexData
        ionex data object with the data of last timeslot replaced by the data of the next day
    """
    tec = ionex.tec
    tec_error = ionex.rms
    tmidx = np.where(np.isclose(ionex.times.mjd - next_day_ionex.times.mjd[0], 0, 1e-6))
    tec[tmidx] = next_day_ionex.tec[0]
    tec_error[tmidx] = next_day_ionex.rms[0]
    return IonexData(
        tec=tec,
        rms=tec_error,
        lons=ionex.lons,
        lats=ionex.lats,
        dims=ionex.dims,
        h=ionex.h,
        times=ionex.times,
    )


def _read_ionex_header(filep: TextIO) -> IonexHeader:
    """Read header from ionex file. Put filepointer to the end of the header."""
    # Declare variables with types
    h1: float | None = None
    h2: float | None = None
    hstep: float | None = None
    start_lon: float | None = None
    end_lon: float | None = None
    step_lon: float | None = None
    start_lat: float | None = None
    end_lat: float | None = None
    step_lat: float | None = None
    start_time: Time | None = None
    ntimes: int | None = None
    step_time: u.Quantity | None = None
    mfactor: float | None = None
    dimension: int | None = None

    filep.seek(0)
    for line in filep:
        if "END OF HEADER" in line:
            break
        label = line[60:-1]
        record = line[:60]
        if "EPOCH OF FIRST MAP" in label:
            yy, mm, day, hr, minute, second = (
                int(float(i)) for i in record.strip().split()
            )
            epoch = Time(f"{yy}-{mm}-{day}T{hr % 24}:{minute}:{second}")
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

    # Check that all optional values are not None
    # Need to do this one by one to get MyPy to understand that they are not None
    # Should probably replace with higher level checks of the header values
    bad_msg = "Could not parse ionex header"
    if h1 is None:
        raise IonexError(bad_msg)
    if h2 is None:
        raise IonexError(bad_msg)
    if hstep is None:
        raise IonexError(bad_msg)
    if start_lon is None:
        raise IonexError(bad_msg)
    if end_lon is None:
        raise IonexError(bad_msg)
    if step_lon is None:
        raise IonexError(bad_msg)
    if start_lat is None:
        raise IonexError(bad_msg)
    if end_lat is None:
        raise IonexError(bad_msg)
    if step_lat is None:
        raise IonexError(bad_msg)
    if start_time is None:
        raise IonexError(bad_msg)
    if ntimes is None:
        raise IonexError(bad_msg)
    if step_time is None:
        raise IonexError(bad_msg)
    if mfactor is None:
        raise IonexError(bad_msg)
    if dimension is None:
        raise IonexError(bad_msg)

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
    data: NDArray[np.float64],
    filep: TextIO,
    stop_label: str,
    timeidx: int,
    ionex_header: IonexHeader,
) -> None:
    """Helper function to parse a data block of a single map in ionex.
    Puts filepointer to the end of the map

    Parameters
    ----------
    data : NDArray
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
    # TODO: Properly type tec
    tec: list[Any] = []
    lonidx = 0
    latidx = 0
    # TODO: Replace magic numbers with named constants
    for line in filep:
        label = line[60:-1]
        if stop_label in label:
            if tec:
                tec_array = np.array(tec) * ionex_header.mfactor
                data[timeidx, lonidx:, latidx] = tec_array
            return
        if "LAT/LON1/LON2/DLON/H" in label:
            if tec:
                tec_array = np.array(tec) * ionex_header.mfactor
                data[timeidx, lonidx:, latidx] = tec_array
            tec = []
            record = line[:60]
            lat, lon1, _, _, _ = (float(record[i : i + 6]) for i in range(2, 32, 6))
            latidx = int(np.argmin(np.abs(ionex_header.lats - lat)))
            lonidx = int(np.argmin(np.abs(ionex_header.lons - lon1)))
        else:
            record = line[:-1]
            tec += [
                float(record[i : i + 5])
                for i in range(0, len(record), 5)
                if record[i : i + 5].strip()
            ]


def _read_ionex_data(filep: TextIO, options: IonexOptions | None = None) -> IonexData:
    """This function parses the IONEX file.
    Some fixed structure (like data records being strings of exactly 80 characters) of the file
    is assumed. This structure is described in Schaer and Gurtner (1998).

    Parameters
    ----------
    filep : TextIO
        pointer to an ionex file
    options : IonexOptions | None, optional
        options for ionex model, by default None

    Returns
    -------
    IonexData
        ionex object
    """
    if options is None:
        options = IonexOptions()
    ionex_header = _read_ionex_header(filep)
    tecarray = np.zeros(
        ionex_header.times.shape + ionex_header.lons.shape + ionex_header.lats.shape,
        dtype=float,
    )
    rmsarray = np.full(tecarray.shape, np.nan)
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
    if options.prefix == "uqr" and options.correct_uqrg_rms:
        # apply uqr correction zhao et al. (2021)
        rmsarray = rmsarray - 6
        rmsarray[rmsarray < 1] = 1

    return IonexData(
        lons=ionex_header.lons,
        lats=ionex_header.lats,
        times=ionex_header.times,
        dims=ionex_header.dims,
        h=ionex_header.h,
        tec=tecarray,
        rms=rmsarray,
    )


def unique_days_from_ionex(ionex_data: IonexData | list[IonexData]) -> Time:
    """Get unique days from a ionex object or list of ionex objects.

    Parameters
    ----------
    ionex_data : IonexData | list[IonexData]
        ionex object or list of ionex objects

    Returns
    -------
    Time
        unique days
    """
    # Get first MJD of each ionex object
    # This avoids issues with midnight crossing
    if isinstance(ionex_data, IonexData):
        time_jd_array = ionex_data.times.sort().mjd[0]
    else:
        time_list: list[NDArray[np.int64]] = []
        for ionex in ionex_data:
            ionex_time = ionex.times.sort().mjd[0]
            time_list.append(ionex_time)
        time_jd_array = np.array(time_list)

    times = Time(
        time_jd_array,
        format="mjd",
    )
    return get_unique_days(times)


def unique_days_from_ionex_files(ionex_files: list[Path] | Path) -> Time:
    """Get unique days from a list of ionex files.

    Parameters
    ----------
    ionex_files : list[Path]
        list of ionex files

    Returns
    -------
    Time
        unique days
    """

    if isinstance(ionex_files, Path):
        ionex_data = read_ionex(ionex_files)
        return unique_days_from_ionex(ionex_data)

    ionex_data_list = [read_ionex(ionex_file) for ionex_file in ionex_files]
    return unique_days_from_ionex(ionex_data_list)
