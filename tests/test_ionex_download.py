# pylint: disable=redefined-outer-name
from __future__ import annotations

import numpy as np
import pytest
from astropy.time import Time
from spinifex.exceptions import IonexError
from spinifex.ionospheric import ionex_download


@pytest.fixture
def times() -> Time:
    times_str = [
        "1999-01-01T00:00:00.00",
        "2010-01-01T00:00:00",
        "2024-02-01T00:00:00",
        "2024-02-01T01:00:00",
    ]
    return Time(times_str)


def test_unique_days(times):
    unique_days = ionex_download.get_unique_days(times)
    assert len(unique_days) == 3


def test_gps_week(times):
    gps_weeks = ionex_download.get_gps_week(times)
    test_weeks = np.array([990, 1564, 2299, 2299])
    assert np.all(gps_weeks == test_weeks)


def test_old_cddis_format(times):
    time = times[0]
    prefix = "cod"
    url = ionex_download.old_cddis_format(time, prefix=prefix)
    assert (
        url
        == "https://cddis.nasa.gov/archive/gnss/products/ionex/1999/001/codg0010.99i.Z"
    )
    prefix = "bad"
    with pytest.raises(IonexError):
        url = ionex_download.old_cddis_format(time, prefix=prefix)

    prefix = "esa"
    url = ionex_download.old_cddis_format(time, prefix=prefix)
    assert (
        url
        == "https://cddis.nasa.gov/archive/gnss/products/ionex/1999/001/esag0010.99i.Z"
    )

    url_stem = "my_stem"
    url = ionex_download.old_cddis_format(time, url_stem=url_stem)
    assert url == "my_stem/1999/001/codg0010.99i.Z"


def test_new_cddis_format(times):
    time = times[-1]
    url = ionex_download.new_cddis_format(time)
    assert (
        url
        == "https://cddis.nasa.gov/archive/gnss/products/ionex/2024/001/IGS0OPSFIN_20240010000_01D_02H_GIM.INX.gz"
    )

    url_stem = "my_stem"
    url = ionex_download.new_cddis_format(time, url_stem=url_stem)
    assert url == "my_stem/2024/001/IGS0OPSFIN_20240010000_01D_02H_GIM.INX.gz"
