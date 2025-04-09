"""Testing of ionospheric indexing tools"""

from __future__ import annotations

import numpy as np
import pytest
from spinifex.exceptions import ArrayShapeError
from spinifex.ionospheric.index_tools import (
    compute_index_and_weights,
    get_indices,
    get_indices_axis,
    get_interpol,
    get_sorted_indices,
)


def test_get_indices():
    """test that indexing tools behave as expected"""
    idx = get_indices(5, np.arange(-10, 10))
    assert idx.idx1 == 15
    assert idx.idx2 == 16
    assert np.isclose(idx.w1, 1.0, 0.001)
    idx = get_indices(355, np.arange(0, 360, 10), wrap_unit=360.0)
    assert idx.idx1 == 0
    assert idx.idx2 == 35
    assert np.isclose(idx.w1, 0.5, 0.001)

    idx = get_indices_axis(np.array(355), np.arange(0, 360, 10), wrap_unit=360.0)
    assert idx.idx1 == 35
    assert idx.idx2 == 0
    assert np.isclose(idx.w1, 0.5, 0.001)
    idx = get_indices_axis(
        np.array([355, 700, 2]), np.arange(0, 360, 10), wrap_unit=360
    )
    assert np.array_equal(idx.idx1, np.array([35, 33, 0]))
    idx = get_indices_axis(
        np.arange(360).reshape(10, 36), np.arange(0, 360, 10), wrap_unit=360
    )
    assert idx.idx1.shape == (10, 36)


def test_get_sorted_indices():
    """test that sorted_indices behaves as expected"""
    with pytest.raises(ArrayShapeError):
        idx = get_sorted_indices(
            lon=350,
            lat=20,
            avail_lat=np.linspace(0, 180, 10),
            avail_lon=np.linspace(0, 360, 36),
            wrap_unit=360,
        )

    idx = get_sorted_indices(
        lon=350,
        lat=20,
        avail_lat=np.linspace(0, 180, 10),
        avail_lon=np.linspace(0, 360, 10),
        wrap_unit=360,
    )
    assert np.array_equal(idx.indices, np.array([0, 1, 2, 3, 7, 8, 6, 9, 5, 4]))


def test_compute_index_and_weights():
    """test that indices and weights are correctly computed"""
    idx = compute_index_and_weights(
        maparray=np.linspace(0, 180, 10), mapvalues=np.linspace(5, 175, 4)
    )
    assert np.array_equal(idx.idx1, np.array([0, 3, 5, 8]))
    assert idx.w2[0] == 0.25
    idx = compute_index_and_weights(maparray=np.linspace(0, 180, 10), mapvalues=55)
    assert idx.idx1.shape == (1,)


def test_get_interpol():
    """tests if interpolation behaves as expected"""
    value = get_interpol(
        data=np.linspace(0, 180, 10), dist=np.abs(np.arange(-50, 50, 10))
    )
    assert value == 100
    value = get_interpol(
        data=np.linspace(0, 180, 10), dist=np.array([0, 1, 2, 3, 4, 0, 1, 2, 3, 4])
    )
    assert value == 50
