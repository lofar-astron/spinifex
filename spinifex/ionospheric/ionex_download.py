from __future__ import annotations

import asyncio
import shutil
from datetime import datetime
from pathlib import Path
from typing import Literal
from urllib.parse import urlparse

import astropy.units as u
import numpy as np
import requests
from astropy.time import Time
from numpy.typing import ArrayLike

from spinifex.exceptions import IonexError, TimeResolutionError
from spinifex.logger import logger

# We need to support downloading from the following sources:
# "cddis.nasa.gov": cddis_nasa_gov,
# "http://chapman.upc.es": chapman_upc_es
# "igsiono.uwm.edu.pl": igsiono_uwm_edu_pl

# TODO: Add support for chapman and igsiono

CENTER_NAMES = {
    "cod",
    "esa",
    "igs",
    "jpl",
    "upc",
}
SERVERS = {
    "cddis": "https://cddis.nasa.gov/archive/gnss/products/ionex",
    "chapman": "",
    "igsiono": "",
}
NAME_SWITCH_WEEK = 2238  # GPS Week where the naming convention changed

SOLUTION = Literal["final", "rapid"]


async def download_file(
    url: str,
    output_file: Path,
    timeout_seconds: int = 30,
    chunk_size: int = 1000,
) -> None:
    """Download a file from a given URL using asyncio.

    Parameters
    ----------
    url : str
        URL to download.
    output_file : Path
        Output file path.
    timeout_seconds : int, optional
        Seconds to wait for request timeout, by default 30
    chunk_size : int, optional
        Chunks of data to download, by default 1000

    Raises
    ------
    IonexError
        If the download times out.
    """
    msg = f"Downloading from {url}"
    logger.info(msg)
    try:
        response = await asyncio.to_thread(requests.get, url, timeout=timeout_seconds)
    except requests.exceptions.Timeout as e:
        msg = "Timed out connecting to server"
        logger.error(msg)
        raise IonexError(msg) from e

    response.raise_for_status()
    with output_file.open("wb") as file_desc:
        for chunk in response.iter_content(chunk_size=chunk_size):
            await asyncio.to_thread(file_desc.write, chunk)


async def download_or_copy_url(
    url: str,
    output_directory: Path | None = None,
    chunk_size: int = 1000,
    timeout_seconds: int = 30,
) -> Path:
    """Download a file from a given URL.

    If the URL is a file URL (i.e. starting with `file://`), it will be copied to the output directory.

    Parameters
    ----------
    url : str
        URL to download.
    output_directory : Path | None, optional
        Output directory, by default None. If None, will default to `ionex_files` in the current working directory.
    chunk_size : int, optional
        Download chunks, by default 1000
    timeout_seconds : int, optional
        Request timeout in seconds, by default 30

    Returns
    -------
    Path
        Output file path

    Raises
    ------
    FileNotFoundError
        If the .netrc file is not found when downloading from CDDIS.
    """
    if output_directory is None:
        output_directory = Path.cwd() / "ionex_files"

    output_directory.mkdir(exist_ok=True)

    url_parsed = urlparse(url)
    url_path = Path(url_parsed.path)
    file_name = url_path.name
    output_file = output_directory / file_name

    if output_file.exists():
        msg = f"File {output_file} already exists. Skipping download."
        logger.info(msg)
        return output_file

    if url_parsed.scheme == "file":
        msg = f"URL scheme {url_parsed.scheme} is not supported"
        logger.info(msg)
        result = await asyncio.to_thread(shutil.copy, url_path, output_file)
        return Path(result)

    if url.startswith("https://cddis.nasa.gov"):
        # CDDIS requires a .netrc file to download
        netrc = Path("~/.netrc").expanduser()
        if not netrc.exists():
            msg = "See: https://cddis.nasa.gov/Data_and_Derived_Products/CreateNetrcFile.html"
            logger.error(msg)
            msg = "Please add your NASA Earthdata login credentials to ~/.netrc"
            logger.error(msg)
            raise FileNotFoundError(msg)

    await download_file(
        url, output_file, timeout_seconds=timeout_seconds, chunk_size=chunk_size
    )

    return output_file


def get_gps_week(time: Time) -> ArrayLike:
    """Get the GPS week from a time.

    Parameters
    ----------
    time : Time
        Time(s) to get the GPS week from

    Returns
    -------
    ArrayLike
        GPS week(s)
    """
    return np.floor((time.gps * u.s).to(u.week).value).astype(int)


def get_unique_days(times: Time) -> Time:
    """Get the unique days from a list of times.

    Parameters
    ----------
    times : Time
        Times to get the unique days from.

    Returns
    -------
    Time
        Unique days
    """
    return Time(np.unique(np.floor(times.mjd)), format="mjd")


def new_cddis_format(
    time: Time,
    prefix: str = "cod",
    url_stem: str | None = None,
    time_resolution: u.Quantity = 2 * u.hour,
    solution: SOLUTION = "final",
) -> str:
    """Get the URL for a new IONEX file from CDDIS.

    Parameters
    ----------
    time : Time
        Time to get the URL for.
    url_stem : str | None, optional
        URL steam, by default None, will default to CDDIS.

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
    dtime: datetime = time.to_datetime()
    prefix = prefix.upper()

    # Parse time resolution
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

    # WWWW/IGS0OPSTYP_YYYYDDDHHMM_01D_SMP_CNT.INX.gz
    file_name = f"{prefix}0OPS{solution_str}_{dtime.year:03d}{dtime.day:03d}0000_01D_{time_res_str}_GIM.INX.gz"
    directory = f"{dtime.year:04d}/{dtime.day:03d}"

    if url_stem is None:
        url_stem = SERVERS.get("cddis")

    return f"{url_stem}/{directory}/{file_name}"


def old_cddis_format(
    time: Time,
    prefix: str = "cod",
    url_stem: str | None = None,
    time_resolution: u.Quantity = 2 * u.hour,
    solution: SOLUTION = "final",
) -> str:
    """Get the URL for an old IONEX file from CDDIS.

    Parameters
    ----------
    time : Time
        Time to get the URL for.
    prefix : str, optional
        Analysis centre prefix, by default "cod"
    url_stem : str | None, optional
        URL steam, by default None, will default to CDDIS.

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
    _ = time_resolution  # Unused
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
    # YYYY/DDD/AAAgDDD#.YYi.Z
    yy = f"{dtime.year:02d}"[-2:]
    file_name = f"{prefix_str}g{dtime.day:03d}0.{yy}i.Z"
    directory_name = f"{dtime.year:04d}/{dtime.day:03d}"
    if url_stem is None:
        url_stem = SERVERS.get("cddis")
    return f"{url_stem}/{directory_name}/{file_name}"


async def download_from_cddis(
    times: Time,
    prefix: str = "cod",
    url_stem: str | None = None,
    time_resolution: u.Quantity = 2 * u.hour,
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
    url_stem : str | None, optional
        URL steam, by default None, will default to CDDIS.
    output_directory : Path | None, optional
        Output directory path, by default None, will default to `ionex_files` in the current working directory.

    Returns
    -------
    list[Path]
        List of downloaded files
    """
    unique_days = get_unique_days(times)

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


def download_ionex(
    server: str,
    times: Time,
    prefix: str = "cod",
    url_stem: str | None = None,
    time_resolution: u.Quantity = 2 * u.hour,
    solution: SOLUTION = "final",
    output_directory: Path | None = None,
) -> list[Path]:
    """Download IONEX files from a server.

    Parameters
    ----------
    server : str
        Server to download from. Must be a supported server.
    times : Time
        Times to download for.
    prefix : str, optional
        Analysis centre prefix, by default "cod". Must be a supported analysis centre.
    url_stem : str | None, optional
        URL stem, by default None, will default to the server URL.
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

    if server not in SERVERS:
        msg = f"Server {server} is not supported. Supported servers are {list(SERVERS.keys())}"
        raise IonexError(msg)

    if server == "cddis":
        return asyncio.run(
            download_from_cddis(
                times,
                prefix=prefix,
                url_stem=url_stem,
                time_resolution=time_resolution,
                solution=solution,
                output_directory=output_directory,
            )
        )

    msg = f"Server {server} is not implemented yet"
    raise NotImplementedError(msg)
