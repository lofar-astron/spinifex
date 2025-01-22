from __future__ import annotations

import filecmp
import shutil
from importlib import resources
from pathlib import Path

from astropy.utils import iers

iers.conf.auto_download = False

import astropy.units as u
import pytest
from astropy.time import Time
from spinifex.exceptions import IonexError, TimeResolutionError
from spinifex.ionospheric import ionex_download
from unlzw3 import unlzw

pytest_plugins = ("pytest_asyncio",)


@pytest.fixture
def times() -> Time:
    times_str = [
        "1999-02-01T00:00:00.00",
        "2010-02-01T00:00:00",
        "2024-02-01T00:00:00",
        "2024-02-01T01:00:00",
    ]
    return Time(times_str)


@pytest.fixture
def old_time() -> Time:
    times_str = [
        "1989-02-01T00:00:00.00",
    ]
    return Time(times_str)


@pytest.fixture
def new_time() -> Time:
    times_str = [
        "2025-01-01T00:00:00.00",
    ]
    return Time(times_str)


@pytest.fixture
def igsiono_time() -> Time:
    times_str = "2024-12-14T00:00:00.00"
    return Time(times_str)


def test_old_cddis_format(times):
    time = times[0]
    prefix = "cod"

    # time resolution is not used!
    for time_resolution in (None, 30 * u.min, 2 * u.hour):
        url = ionex_download.old_cddis_format(
            time, prefix=prefix, time_resolution=time_resolution
        )
        assert (
            url
            == "https://cddis.nasa.gov/archive/gnss/products/ionex/1999/032/codg0320.99i.Z"
        )
    prefix = "bad"
    with pytest.raises(IonexError):
        url = ionex_download.old_cddis_format(time, prefix=prefix)

    prefix = "esa"
    url = ionex_download.old_cddis_format(time, prefix=prefix)
    assert (
        url
        == "https://cddis.nasa.gov/archive/gnss/products/ionex/1999/032/esag0320.99i.Z"
    )

    url_stem = "my_stem"
    url = ionex_download.old_cddis_format(time, url_stem=url_stem)
    assert url == "my_stem/1999/032/codg0320.99i.Z"

    url = ionex_download.old_cddis_format(time, prefix="cod", solution="rapid")
    assert (
        url
        == "https://cddis.nasa.gov/archive/gnss/products/ionex/1999/032/corg0320.99i.Z"
    )


def test_new_cddis_format(times):
    time = times[-1]
    url = ionex_download.new_cddis_format(time)
    assert (
        url
        == "https://cddis.nasa.gov/archive/gnss/products/ionex/2024/032/COD0OPSFIN_20240320000_01D_01H_GIM.INX.gz"
    )

    url_stem = "my_stem"
    url = ionex_download.new_cddis_format(time, url_stem=url_stem)
    assert url == "my_stem/2024/032/COD0OPSFIN_20240320000_01D_01H_GIM.INX.gz"

    time_resolution = 2 * u.hour
    url = ionex_download.new_cddis_format(time, time_resolution=time_resolution)
    assert (
        url
        == "https://cddis.nasa.gov/archive/gnss/products/ionex/2024/032/COD0OPSFIN_20240320000_01D_02H_GIM.INX.gz"
    )

    time_resolution = 30 * u.min
    url = ionex_download.new_cddis_format(time, time_resolution=time_resolution)
    assert (
        url
        == "https://cddis.nasa.gov/archive/gnss/products/ionex/2024/032/COD0OPSFIN_20240320000_01D_30M_GIM.INX.gz"
    )

    url = ionex_download.new_cddis_format(time, solution="rapid")
    assert (
        url
        == "https://cddis.nasa.gov/archive/gnss/products/ionex/2024/032/COD0OPSRAP_20240320000_01D_01H_GIM.INX.gz"
    )

    url = ionex_download.new_cddis_format(time, prefix="esa")
    assert (
        url
        == "https://cddis.nasa.gov/archive/gnss/products/ionex/2024/032/ESA0OPSFIN_20240320000_01D_02H_GIM.INX.gz"
    )

    with pytest.raises(IonexError):
        url = ionex_download.new_cddis_format(time, prefix="bad")

    with pytest.raises(TimeResolutionError):
        url = ionex_download.new_cddis_format(time, time_resolution=1.5 * u.min)


def test_chapman_format(times):
    base_url = "http://chapman.upc.es/tomion/rapid"
    expected_urls = [
        "1999/032_990201.15min/uqrg0320.99i.Z",
        "2010/032_100201.15min/uqrg0320.10i.Z",
        "2024/032_240201.15min/uqrg0320.24i.Z",
        "2024/032_240201.15min/uqrg0320.24i.Z",
        "2025/001_250101.15min/uqrg0010.25i.Z",
    ]
    for time, expected_url in zip(times, expected_urls):
        url = ionex_download.chapman_format(time, url_stem=base_url)
        assert url == f"{base_url}/{expected_url}"


@pytest.mark.asyncio
async def test_chapman_download(tmpdir, old_time, new_time):
    with pytest.raises(IonexError):
        await ionex_download.download_from_chapman(
            times=old_time,
            url_stem=None,
            output_directory=Path(tmpdir),
        )

    downloaded_files = await ionex_download.download_from_chapman(
        times=new_time,
        url_stem=None,
        output_directory=Path(tmpdir),
    )
    downloaded_file = downloaded_files[0]

    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        expected_file_compressed = datapath / "uqrg0010.25i.truncated.zip"
        shutil.unpack_archive(expected_file_compressed, tmpdir)

    expected_file = Path(tmpdir) / "uqrg0010.25i.truncated"

    downloaded_data = (
        unlzw(
            downloaded_file.read_bytes(),
        )
        .decode("utf-8")
        .splitlines()
    )

    expected_data = expected_file.read_text().splitlines()
    n_lines_tructed = len(expected_data) - 1  # Last line is truncated

    # Compare first 49 lines of the downloaded file
    for i, (line1, line2) in enumerate(zip(downloaded_data, expected_data)):
        if i == n_lines_tructed:
            break
        assert line1 == line2

    # Clean up
    downloaded_file.unlink(missing_ok=True)
    expected_file.unlink(missing_ok=True)


def test_igsiono_format(igsiono_time):
    url = ionex_download.igsiono_format(igsiono_time)
    assert (
        url
        == "ftp://igs-final.man.olsztyn.pl/pub/gps_data/GPS_IONO/cmpcmb/24349/IGS0OPSFIN_20243490000_01D_02H_GIM.INX.gz"
    )


@pytest.mark.asyncio
async def test_igsiono_download(tmpdir, igsiono_time):
    downloaded_files = await ionex_download.download_from_igsiono(
        times=igsiono_time,
        output_directory=Path(tmpdir),
    )
    downloaded_file = downloaded_files[0]

    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        expected_file = datapath / "IGS0OPSFIN_20243490000_01D_02H_GIM.INX.gz"

    assert filecmp.cmp(downloaded_file, expected_file)

    # Clean up
    downloaded_file.unlink(missing_ok=True)


def test_download_ionex_igsiono(tmpdir, igsiono_time):
    downloaded_files = ionex_download.download_ionex(
        server="igsiono",
        times=igsiono_time,
        prefix="igs",
        url_stem=None,
        time_resolution=None,
        solution="final",
        output_directory=Path(tmpdir),
    )
    downloaded_file = downloaded_files[0]

    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        expected_file = datapath / "IGS0OPSFIN_20243490000_01D_02H_GIM.INX.gz"

    assert filecmp.cmp(downloaded_file, expected_file)

    # Clean up
    downloaded_file.unlink(missing_ok=True)


def test_download_ionex_chapman(tmpdir, new_time):
    downloaded_files = ionex_download.download_ionex(
        server="chapman",
        times=new_time,
        prefix="uqr",
        url_stem=None,
        time_resolution=None,
        solution="final",
        output_directory=Path(tmpdir),
    )
    downloaded_file = downloaded_files[0]

    with resources.as_file(resources.files("spinifex.data.tests")) as datapath:
        expected_file_compressed = datapath / "uqrg0010.25i.truncated.zip"
        shutil.unpack_archive(expected_file_compressed, tmpdir)

    expected_file = Path(tmpdir) / "uqrg0010.25i.truncated"

    downloaded_data = (
        unlzw(
            downloaded_file.read_bytes(),
        )
        .decode("utf-8")
        .splitlines()
    )

    expected_data = expected_file.read_text().splitlines()
    n_lines_tructed = len(expected_data) - 1  # Last line is truncated

    # Compare first 49 lines of the downloaded file
    for i, (line1, line2) in enumerate(zip(downloaded_data, expected_data)):
        if i == n_lines_tructed:
            break
        assert line1 == line2

    # Clean up
    downloaded_file.unlink(missing_ok=True)
    expected_file.unlink(missing_ok=True)
