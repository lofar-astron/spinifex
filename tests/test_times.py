from __future__ import annotations

import numpy as np
import pytest
from astropy.utils import iers

iers.conf.auto_download = False

from astropy.time import Time
from spinifex.times import (
    get_consecutive_days,
    get_gps_week,
    get_indexlist_unique_days,
    get_unique_days,
    get_unique_days_index,
)


@pytest.fixture
def times() -> Time:
    times_str = [
        "1999-02-01T00:00:00.00",
        "2010-02-01T00:00:00",
        "2010-02-01T14:27:01",
        "2024-02-01T00:00:00",
        "2024-02-01T01:00:00",
    ]
    return Time(times_str)


def test_unique_days(times) -> None:
    unique_days = get_unique_days(times)
    assert len(unique_days) == 3


def test_gps_week(times):
    gps_weeks = get_gps_week(times)
    test_weeks = np.array([995, 1569, 1569, 2299, 2299])
    assert np.all(gps_weeks == test_weeks)


def test_unique_days_index(times):
    unique_days_index = get_unique_days_index(times)
    test_index = np.array([0, 1, 1, 2, 2])
    assert np.all(unique_days_index == test_index)


def test_get_indexlist_unique_days(times):
    unique_days = get_unique_days(times)
    test_index = np.zeros((3, 5), dtype=bool)
    test_index[0, 0] = True
    test_index[1, (1, 2)] = True
    test_index[2, 3:] = True
    indexlist_unique_days = get_indexlist_unique_days(unique_days, times)
    assert np.all(indexlist_unique_days == test_index)


def test_consecutive_days(times):
    unique_days = get_unique_days(times)
    consecutive_days_index = get_consecutive_days(unique_days)
    test_index = np.array([0, 1, 2])
    assert np.all(consecutive_days_index == test_index)
