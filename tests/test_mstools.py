"""Testing of mstools"""

from __future__ import annotations

import shutil
from importlib import resources
from pathlib import Path

from astropy.utils import iers

iers.conf.auto_download = False

import importlib.util

import astropy.units as u
import pytest

if not importlib.util.find_spec("casacore"):
    pytest.skip(
        "Skipping tests of ms_tools because casacore is not available.",
        allow_module_level=True,
    )

from spinifex.vis_tools.ms_tools import (
    get_columns_from_ms,
    get_dtec_from_ms,
    get_rm_from_ms,
)


@pytest.fixture
def unzip_ms(tmpdir) -> Path:  # type: ignore[misc]
    with resources.as_file(resources.files("spinifex.data.tests")) as test_data:
        zipped_ms = test_data / "test.ms.zip"
    shutil.unpack_archive(zipped_ms, tmpdir)

    yield Path(tmpdir / "test.MS")

    shutil.rmtree(tmpdir / "test.MS")


def test_unzip_worked(unzip_ms: Path):
    # Check that the unzipped directory exists
    assert unzip_ms.exists()


def test_mstools(unzip_ms: Path) -> None:
    cols = get_columns_from_ms(unzip_ms)
    assert "ANTENNA1" in cols
    with resources.as_file(resources.files("spinifex.data.tests")) as test_data:
        rms = get_rm_from_ms(
            unzip_ms,
            output_directory=test_data,
            prefix="esa",
            server="cddis",
            use_stations=["CS002HBA0"],
            timestep=20 * u.s,
        )
        assert "CS002HBA0" in rms
        dtec = get_dtec_from_ms(
            unzip_ms,
            output_directory=test_data,
            prefix="esa",
            server="cddis",
            use_stations=["CS002HBA0"],
            timestep=20 * u.s,
        )
        assert "CS002HBA0" in dtec


def test_station_selection(unzip_ms: Path) -> None:
    with resources.as_file(resources.files("spinifex.data.tests")) as test_data:
        rms = get_rm_from_ms(
            unzip_ms,
            use_stations=["CS002HBA0"],
            timestep=20 * u.s,
            output_directory=test_data,
            prefix="esa",
            server="cddis",
        )
        assert "CS002HBA0" in rms
        assert len(rms) == 1

        rms = get_rm_from_ms(
            unzip_ms,
            use_stations="all",
            timestep=20 * u.s,
            output_directory=test_data,
            prefix="esa",
            server="cddis",
        )
        stations = ["CS002HBA0", "RS210HBA", "RS509HBA"]
        assert all(station in rms for station in stations)
        assert len(rms) == 3

        rms = get_rm_from_ms(
            unzip_ms,
            use_stations="average",
            timestep=20 * u.s,
            output_directory=test_data,
            prefix="esa",
            server="cddis",
        )
        assert "average_station_pos" in rms
        assert len(rms) == 1
