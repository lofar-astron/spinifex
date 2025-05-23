"""Utilities for downloading IONEX files from various servers."""

from __future__ import annotations

import asyncio
from datetime import datetime
from pathlib import Path

import astropy.units as u
from astropy.time import Time

from spinifex.asyncio_wrapper import sync_wrapper
from spinifex.download import download_or_copy_url
from spinifex.exceptions import IonexError, TimeResolutionError
from spinifex.ionospheric.tec_data import (
    CENTER_NAMES,
    DEFAULT_TIME_RESOLUTIONS,
    NAME_SWITCH_WEEK,
    SOLUTION,
    IonexOptions,
    Servers,
)
from spinifex.logger import logger
from spinifex.times import get_gps_week, get_unique_days

# We need to support downloading from the following sources:
# "cddis.nasa.gov": cddis_nasa_gov,
# "chapman.upc.es": chapman_upc_es
# "igsiono.uwm.edu.pl": igsiono_uwm_edu_pl


def new_cddis_format(
    time: Time,
    prefix: str = "cod",
    time_resolution: u.Quantity | None = None,
    url_stem: str | None = None,
    solution: SOLUTION = "final",
) -> str:
    """Get the URL for a new IONEX file from CDDIS.

    Parameters
    ----------
    time : Time
        Time to get the URL for.
    prefix : str, optional
        Analysis centre prefix, by default "cod". Must be a supported analysis centre.
    time_resolution : u.Quantity | None, optional
        Time resolution, by default None, will default to the server time resolution.
    url_stem : str | None, optional
        URL steam, by default None, will default to CDDIS.
    solution : SOLUTION, optional
        Solution type, by default "final", must be "final" or "rapid".

    Returns
    -------
    str
        File URL
    """
    # Name Format Since GPS Week 2238
    # WWWW/IGS0OPSTYP_YYYYDDDHHMM_01D_SMP_CNT.INX.gz
    # Code	Meaning
    # WWWW	GPS week
    # TYP	solution type identifier -FIN (Final solution combination)
    # RAP (Rapid solution combination
    # YYYY	4-digit year
    # DDD	3-digit day of year
    # HH	2-digit hour
    # MM	2-digit minute
    # SMP	temporal product sampling resolution
    # CNT	content type -GIM (global ionosphere (TEC) maps)
    # ROT (rate of TEC index maps)
    # .gz	gzip compressed file
    assert time.isscalar, "Only one time is supported"
    if prefix not in CENTER_NAMES:
        msg = f"prefix {prefix} is not supported. Supported prefixes are {CENTER_NAMES}"
        raise IonexError(msg)

    # Parse time resolution
    if time_resolution is None:
        time_resolution = DEFAULT_TIME_RESOLUTIONS.get(prefix)
        if time_resolution is None:
            msg = f"Time resolution not defined for {prefix}"
            raise TimeResolutionError(msg)
        msg = f"Using default time resolution {time_resolution} for {prefix}"
        logger.info(msg)
    if not time_resolution.value.is_integer():
        error = f"Time resolution must be an integer. Got {time_resolution}"
        raise TimeResolutionError(error)
    if time_resolution.to(u.min).value % 15 != 0:
        msg = f"Time resolution on CDDIS is multiples 15 minutes. Please check the time resolution ({time_resolution})."
        logger.warning(msg)
    if time_resolution < 1 * u.hour:
        time_resolution = time_resolution.to(u.min)
    else:
        time_resolution = time_resolution.to(u.hour)
    time_res_str = (
        f"{int(time_resolution.value):02d}{str(time_resolution.unit)[0].upper()}"
    )

    if solution == "final":
        solution_str = "FIN"
    elif solution == "rapid":
        solution_str = "RAP"
    else:
        msg = f"Solution {solution} is not supported. Supported solutions are ['final', 'rapid']"  # type: ignore[unreachable]
        raise IonexError(msg)

    dtime: datetime = time.to_datetime()
    doy = time.datetime.timetuple().tm_yday
    prefix_formatted = prefix.upper()

    # WWWW/IGS0OPSTYP_YYYYDDDHHMM_01D_SMP_CNT.INX.gz
    file_name = f"{prefix_formatted}0OPS{solution_str}_{dtime.year:03d}{doy:03d}0000_01D_{time_res_str}_GIM.INX.gz"
    directory = f"{dtime.year:04d}/{doy:03d}"

    if url_stem is None:
        url_stem = Servers.get_url(Servers.CDDIS)

    return f"{url_stem}/{directory}/{file_name}"


def old_cddis_format(
    time: Time,
    prefix: str = "cod",
    time_resolution: u.Quantity | None = None,
    url_stem: str | None = None,
    solution: SOLUTION = "final",
) -> str:
    """Get the URL for an old IONEX file from CDDIS.

    Parameters
    ----------
    time : Time
        Time to get the URL for.
    prefix : str, optional
        Analysis centre prefix, by default "cod"
    time_resolution : u.Quantity | None, optional
        Time resolution, by default None, will default to the server time resolution.
    url_stem : str | None, optional
        URL steam, by default None, will default to CDDIS.
    solution : SOLUTION, optional
        Solution type, by default "final", must be "final" or "rapid".

    Returns
    -------
    str
        File URL

    Raises
    ------
    IonexError
        If the prefix is not a supported analysis centre.
    """
    # Name Format Before GPS Week 2237
    # YYYY/DDD/AAAgDDD#.YYi.Z - Vertical total electron content (TEC) maps
    # Code	Meaning
    # YYYY	4-digit year
    # DDD	3-digit day of year
    # AAA	Analysis center name
    # #	file number for the day, typically 0
    # YY	2-digit year
    # .Z	Unix compressed file
    if time_resolution is not None:
        msg = f"Time resolution is not used by the old CDDIS format (got `{time_resolution}`). Ignoring."
        logger.warning(msg)

    assert time.isscalar, "Only one time is supported"
    if prefix not in CENTER_NAMES:
        msg = f"prefix {prefix} is not supported. Supported prefixes are {CENTER_NAMES}"
        raise IonexError(msg)

    if solution == "final":
        prefix_str = prefix
    elif solution == "rapid":
        prefix_str = prefix[:-1] + "r"
    else:
        msg = f"Solution {solution} is not supported. Supported solutions are ['final', 'rapid']"  # type: ignore[unreachable]
        raise IonexError(msg)

    dtime: datetime = time.to_datetime()
    doy = time.datetime.timetuple().tm_yday
    # YYYY/DDD/AAAgDDD#.YYi.Z
    yy = f"{dtime.year:02d}"[-2:]
    file_name = f"{prefix_str}g{doy:03d}0.{yy}i.Z"
    directory_name = f"{dtime.year:04d}/{doy:03d}"
    if url_stem is None:
        url_stem = Servers.get_url(Servers.CDDIS)

    return f"{url_stem}/{directory_name}/{file_name}"


async def download_from_cddis(
    times: Time,
    prefix: str = "cod",
    time_resolution: u.Quantity | None = None,
    url_stem: str | None = None,
    solution: SOLUTION = "final",
    output_directory: Path | None = None,
) -> list[Path]:
    """Download IONEX files from CDDIS.

    Parameters
    ----------
    times : Time
        Times to download for.
    prefix : str, optional
        Analysis centre prefix, by default "cod"
    time_resolution : u.Quantity | None, optional
        Time resolution, by default None, will default to the server time resolution.
    url_stem : str | None, optional
        URL steam, by default None, will default to CDDIS.
    solution : SOLUTION, optional
        Solution type, by default "final", must be "final" or "rapid".
    output_directory : Path | None, optional
        Output directory path, by default None, will default to `ionex_files` in the current working directory.

    Returns
    -------
    list[Path]
        List of downloaded files
    """
    unique_days: Time = get_unique_days(times)

    coros = []
    for day in unique_days:
        if get_gps_week(day) >= NAME_SWITCH_WEEK:
            formatter = new_cddis_format
        else:
            formatter = old_cddis_format
        url = formatter(
            day,
            prefix=prefix,
            url_stem=url_stem,
            time_resolution=time_resolution,
            solution=solution,
        )
        coros.append(download_or_copy_url(url, output_directory=output_directory))

    return await asyncio.gather(*coros)


def chapman_format(
    time: Time,
    prefix: str = "uqr",
    url_stem: str | None = None,
) -> str:
    """Format a URL for an IRTG file from Chapman.

    Parameters
    ----------
    time : Time
        Time to get the URL for.
    prefix : str, optional
        Ionex prefix, by default "irt"
    url_stem : str | None, optional
        URL stem, by default None

    Returns
    -------
    str
        File URL
    """

    assert time.isscalar, "Only one time is supported"

    dtime: datetime = time.to_datetime()
    doy = int(time.datetime.timetuple().tm_yday)
    yy = f"{dtime.year:02d}"[-2:]
    # e.g. uqrg3260.96i.Z
    file_name = f"{prefix}g{doy:03d}0.{yy}i.Z"
    # e.g. 1996/326_961121.15min/
    # Hard coded to 15 minute resolution - Chapman only provides 15 minute resolution right now
    # Could be updated to be more flexible in the future
    directory_name = (
        f"{dtime.year:04d}/{doy:03d}_{yy}{dtime.month:02d}{dtime.day:02d}.15min"
    )
    if url_stem is None:
        url_stem = Servers.get_url(Servers.CHAPMAN)

    return f"{url_stem}/{directory_name}/{file_name}"


async def download_from_chapman(
    times: Time,
    prefix: str = "uqr",
    url_stem: str | None = None,
    output_directory: Path | None = None,
) -> list[Path]:
    """Download IONEX files from Chapman.

    Parameters
    ----------
    times : Time
        Times to download for.
    prefix : str, optional
        Ionex prefix, by default "irt"
    url_stem : str | None, optional
        URL stem, by default None
    output_directory : Path | None, optional
        Output directory path, by default None

    Returns
    -------
    list[Path]
        List of downloaded file paths

    Raises
    ------
    IonexError
        If the date is before 2021/024.
    """
    unique_days: Time = get_unique_days(times)
    for day in unique_days:
        # Chapman goes back to 326th day of 1996 == 21/11/1996
        _chapman_start = datetime(1996, 11, 21)
        if day.datetime < _chapman_start:
            msg = f"Date {day.datetime} is before {_chapman_start}. Chapman only goes back to this date."
            raise IonexError(msg)

    coros = []
    for day in unique_days:
        url = chapman_format(day, prefix=prefix, url_stem=url_stem)
        coros.append(download_or_copy_url(url, output_directory=output_directory))

    return await asyncio.gather(*coros)


def igsiono_format(
    time: Time,
    prefix: str = "igs",
    time_resolution: u.Quantity | None = None,
    url_stem: str | None = None,
    solution: SOLUTION = "final",
) -> str:
    """Format a URL for an IGS IONO file from IGSIONO.

    Parameters
    ----------
    time : Time
        Time to get the URL for.
    prefix : str, optional
        Ionex prefix, by default "igs"
    time_resolution : u.Quantity | None, optional
        Ionex resolution, by default None
    url_stem : str | None, optional
        URL stem, by default None
    solution : SOLUTION, optional
        Type of solution, by default "final"

    Returns
    -------
    str
        File URL
    """
    assert time.isscalar, "Only one time is supported"

    dtime: datetime = time.to_datetime()
    doy = time.datetime.timetuple().tm_yday
    yy = f"{dtime.year:02d}"[-2:]
    directory_name = f"pub/gps_data/GPS_IONO/cmpcmb/{yy}{doy:03d}"
    if url_stem is None:
        url_stem = Servers.get_url(Servers.IGSIONO)

    # File name matches CDDIS format
    if get_gps_week(time) >= NAME_SWITCH_WEEK:
        formatter = new_cddis_format
    else:
        formatter = old_cddis_format
    file_name = Path(
        formatter(
            time,
            prefix=prefix,
            url_stem=url_stem,
            time_resolution=time_resolution,
            solution=solution,
        )
    ).name
    return f"{url_stem}/{directory_name}/{file_name}"


async def download_from_igsiono(
    times: Time,
    prefix: str = "igs",
    time_resolution: u.Quantity | None = None,
    url_stem: str | None = None,
    solution: SOLUTION = "final",
    output_directory: Path | None = None,
) -> list[Path]:
    """Download IONEX files from IGSIONO.

    Parameters
    ----------
    times : Time
        Times to download for.
    prefix : str, optional
        URL prefix, by default "igs"
    time_resolution : u.Quantity | None, optional
        Ionex resolution, by default None
    url_stem : str | None, optional
        URL stem, by default None
    solution : SOLUTION, optional
        Type of solution, by default "final"
    output_directory : Path | None, optional
        Output directory path, by default None

    Returns
    -------
    list[Path]
        List of downloaded file paths
    """
    unique_days: Time = get_unique_days(times)

    coros = []
    for day in unique_days:
        url = igsiono_format(
            day,
            prefix=prefix,
            time_resolution=time_resolution,
            url_stem=url_stem,
            solution=solution,
        )
        coros.append(download_or_copy_url(url, output_directory=output_directory))

    return await asyncio.gather(*coros)


async def _download_ionex_coro(
    times: Time,
    options: IonexOptions | None = None,
) -> list[Path]:
    """Download IONEX files from a server.

    Parameters
    ----------
    times : Time
        Times to download for.
    options : IonexOptions
        Options for the ionex model.

    Returns
    -------
    list[Path]
        List of downloaded files

    Raises
    ------
    IonexError
        If the server is not supported.
    NotImplementedError
        If the server is not implemented yet.
    """
    if options is None:
        options = IonexOptions()

    # TODO: Consider refactor to pass around the options object
    if options.server == Servers.CDDIS:
        return await download_from_cddis(
            times,
            prefix=options.prefix,
            url_stem=options.url_stem,
            time_resolution=options.time_resolution,
            solution=options.solution,
            output_directory=options.output_directory,
        )
    if options.server == Servers.CHAPMAN:
        if options.prefix != "uqr":
            msg = f"Chapman server primarily supports uqr prefix. You provided {options.prefix}, this may fail"
            logger.warning(msg)

        return await download_from_chapman(
            times,
            prefix=options.prefix,
            url_stem=options.url_stem,
            output_directory=options.output_directory,
        )

    # Must be IGSIONO
    assert options.server == Servers.IGSIONO, "Server not supported"
    if options.prefix != "igs":
        msg = f"IGS IONO server primarily supports igs prefix. You provided {options.prefix}, this may fail"
        logger.warning(msg)
    return await download_from_igsiono(
        times,
        prefix=options.prefix,
        time_resolution=options.time_resolution,
        url_stem=options.url_stem,
        solution=options.solution,
        output_directory=options.output_directory,
    )


async def download_ionex_coro(
    times: Time,
    server: str = "chapman",
    prefix: str = "uqr",
    url_stem: str | None = None,
    time_resolution: u.Quantity | None = None,
    solution: SOLUTION = "final",
    output_directory: Path | None = None,
) -> list[Path]:
    """Download IONEX files from a server.

    Parameters
    ----------
    times : Time
        Times to download for.
    server : str
        Server to download from, by default "cddis". Must be a supported server.
    prefix : str, optional
        Analysis centre prefix, by default "cod". Must be a supported analysis centre.
    url_stem : str | None, optional
        URL stem, by default None, will default to the server URL.
    time_resolution : u.Quantity | None, optional
        Time resolution, by default None, will default to the server time resolution.
    solution : SOLUTION, optional
        Solution type, by default "final", must be "final" or "rapid".
    output_directory : Path | None, optional
        Output directory path, by default None, will default to `ionex_files` in the current working directory.

    Returns
    -------
    list[Path]
        List of downloaded files

    Raises
    ------
    IonexError
        If the server is not supported.
    NotImplementedError
        If the server is not implemented yet.
    """

    options = IonexOptions(
        server=server,
        prefix=prefix,
        url_stem=url_stem,
        time_resolution=time_resolution,
        solution=solution,
        output_directory=output_directory,
        correct_uqrg_rms=False,
        height=350 * u.km,
        remove_midnight_jumps=False,
    )

    return await _download_ionex_coro(times, options)


# Create synchronous wrapper of download_ionex_coro
download_ionex = sync_wrapper(download_ionex_coro)
_download_ionex = sync_wrapper(_download_ionex_coro)
