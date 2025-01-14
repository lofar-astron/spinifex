#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0
# pylint: disable=duplicate-code
# pylint: disable=redefined-outer-name
"""Testing of mstools"""

from __future__ import annotations

import shutil
from importlib import resources
from pathlib import Path

from astropy.utils import iers  # pylint: disable=ungrouped-imports

iers.conf.auto_download = False

import astropy.units as u
import pytest
from spinifex.vis_tools.ms_tools import get_columns_from_ms, get_rm_from_ms


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
            iono_kwargs={
                "output_directory": test_data,
                "prefix": "esa",
                "server": "cddis",
            },
            use_stations=["CS002HBA0"],
            timestep=20 * u.s,
        )
    assert "CS002HBA0" in rms
