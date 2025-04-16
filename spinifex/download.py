"""Generic download utils"""

from __future__ import annotations

import asyncio
import shutil
from ftplib import FTP
from pathlib import Path
from urllib.parse import urlparse

import requests

from spinifex.exceptions import IonexError
from spinifex.logger import logger


def _ftp_download_and_quit(ftp: FTP, file_name: str, output_file: Path) -> None:
    """Download a file from an FTP server and quit the connection.

    Parameters
    ----------
    ftp : FTP
        FTP connection.
    file_name : str
        File name to download.
    output_file : Path
        Output file path.

    Raises
    ------
    e
        If the download fails.
    """
    try:
        with output_file.open("wb") as file_desc:
            ftp.retrbinary(f"RETR {file_name}", file_desc.write)
    except Exception as e:
        output_file.unlink(missing_ok=True)
        raise e
    finally:
        ftp.quit()


async def download_file_ftp(
    url: str,
    output_file: Path,
) -> None:
    """Download a file from a given URL using asyncio.

    Parameters
    ----------
    url : str
        URL to download.
    output_file : Path
        Output file path.
    """
    url_parsed = urlparse(url)
    url_path = Path(url_parsed.path)
    file_name = url_path.name
    directory_name = url_path.parent.as_posix()[1:]  # Remove leading slash

    ftp = FTP(url_parsed.netloc)
    # Anonymous login
    ftp.login()
    ftp.cwd(directory_name)

    await asyncio.to_thread(_ftp_download_and_quit, ftp, file_name, output_file)


async def download_file_http(
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

    if url_parsed.scheme == "ftp":
        await download_file_ftp(url, output_file)
        return output_file

    if url.startswith("https://cddis.nasa.gov"):
        # CDDIS requires a .netrc file to download
        netrc = Path("~/.netrc").expanduser()
        if not netrc.exists():
            msg = "See: https://cddis.nasa.gov/Data_and_Derived_Products/CreateNetrcFile.html"
            logger.error(msg)
            msg = "Please add your NASA Earthdata login credentials to ~/.netrc"
            logger.error(msg)
            raise FileNotFoundError(msg)

    await download_file_http(
        url, output_file, timeout_seconds=timeout_seconds, chunk_size=chunk_size
    )

    return output_file
