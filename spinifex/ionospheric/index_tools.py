"""Some useful multipurpose functions to interpolate on grids"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from numpy.typing import NDArray

from spinifex.exceptions import ArrayShapeError


class Indices(NamedTuple):
    """Indices of the closest two points in a possibly wrapping selection and the inverse distance weights"""

    idx1: NDArray[np.int64]
    """Index of the first closest point"""
    idx2: NDArray[np.int64]
    """Index of the second closest point"""
    w1: NDArray[np.float64]
    """Weight of the first closest point"""
    w2: NDArray[np.float64]
    """Weight of the second closest point"""


class SortedIndices(NamedTuple):
    """Indices of the closest two points in a possibly wrapping selection"""

    indices: NDArray[np.int64]
    """Indices sorted by distance"""
    distance: NDArray[np.float64]
    """Distances"""


class Weights(NamedTuple):
    """Weights of the closest two points in a possibly wrapping selection"""

    w1: NDArray[np.float64]
    """Weight of the first closest point"""
    w2: NDArray[np.float64]
    """Weight of the second closest point"""


def get_indices_axis(
    goal: NDArray[np.float64], selection: NDArray[np.float64], wrap_unit: float = 0
) -> Indices:
    """get indices of the closest two points in a possibly wrapping selection for an array
    of goals

    Parameters
    ----------
    goal : NDArray[np.float64]
        array of points for which to get the indices
    selection : NDArray[np.float64]
        array from which to get the indices
    wrap_unit : float, optional
        set if selection is a wrapping entity (e.g. angles), by default 0

    Returns
    -------
    Indices
        object with indices and weights
    """

    if wrap_unit > 0:
        goal = np.remainder(goal + 0.5 * wrap_unit, wrap_unit) - 0.5 * wrap_unit
        idx1 = np.argmin(
            np.absolute(
                np.remainder(
                    goal[..., np.newaxis] - selection + 0.5 * wrap_unit, wrap_unit
                )
                - 0.5 * wrap_unit
            ),
            axis=-1,
        ).astype(np.int64)
    else:
        idx1 = np.argmin(np.absolute(goal[..., np.newaxis] - selection), axis=-1)
    if np.isscalar(idx1):
        if goal < selection[idx1]:
            idx2 = idx1
            idx1 -= 1
            if idx1 < 0:
                idx1 = selection.shape[0] - 1 if wrap_unit > 0 else 0
        else:
            idx2 = idx1 + 1
            if idx2 >= selection.shape[0]:
                if wrap_unit > 0:
                    idx2 = 0
                else:
                    idx2 -= 1
    else:
        idx2 = np.copy(idx1)
        idx1[goal < selection[idx1]] = idx1[goal < selection[idx1]] - 1
        idx2[goal >= selection[idx2]] = idx2[goal >= selection[idx2]] + 1
        if wrap_unit > 0:
            idx1[idx1 < 0] = selection.shape[0] - 1
            idx2[idx2 >= selection.shape[0]] = 0
        else:
            idx1[idx1 < 0] = 0
            idx2[idx2 >= selection.shape[0]] -= 1
    weights = _get_weights(goal, idx1, idx2, selection, wrap_unit)
    return Indices(idx1=idx1, idx2=idx2, w1=weights.w1, w2=weights.w2)


def _get_weights(
    goal: NDArray[float],
    index1: NDArray[np.int64],
    index2: NDArray[np.int64],
    selection: NDArray[np.float64],
    wrap_unit: float = 0,
) -> Weights:
    """Calculate weights based on distance of goal to selection

    Parameters
    ----------
    goal : NDArray
        array of points to get weights for
    index1 : NDArray
        indices in selection for goals (index1, index2) per goal
    index2 : NDArray
        indices in selection for goals (index1, index2) per goal
    selection : NDArray
        array to select from
    wrap_unit : float, optional
        if goal/selection is a wrapable (e.g. angle) set this unit (e.g. 360), by default 0

    Returns
    -------
    namedtuple
        tuple with weights
    """
    if wrap_unit > 0:
        distance1 = np.abs(wrap_around_zero(selection[index1] - goal, wrap_unit))
        distance2 = np.abs(wrap_around_zero(selection[index2] - goal, wrap_unit))
    else:
        distance1 = np.abs(selection[index1] - goal)
        distance2 = np.abs(selection[index2] - goal)

    sumdist = distance1 + distance2
    if np.any(sumdist == 0):
        # indices are equal and distance = 0
        distance1[sumdist == 0] = 0.5
        distance2[sumdist == 0] = 0.5
        sumdist[sumdist == 0] = 1
    return Weights(w1=1 - distance1 / sumdist, w2=1 - distance2 / sumdist)


def get_indices(
    goal: float, selection: NDArray[np.float64], wrap_unit: float = 0
) -> Indices:
    """find the indices of the closest two points in a possibly wrapping array selection

    Parameters
    ----------
    goal : float
        location of point
    selection : NDArray
        array of points
    wrap_unit : float, optional
        if goal/selection is a wrapping entity (e.g. angles) set this to the wrap value (e.g. 360),
        by default 0

    Returns
    -------
    Indices:
        sorted list of index1 and index2
    """

    if wrap_unit > 0:
        goal = np.remainder(goal + 0.5 * wrap_unit, wrap_unit) - 0.5 * wrap_unit
        idx1 = np.argmin(
            np.absolute(
                np.remainder(goal - selection + 0.5 * wrap_unit, wrap_unit)
                - 0.5 * wrap_unit
            )
        )
        idx2 = (
            idx1 - 1
            if np.remainder(goal - selection[idx1] + 0.5 * wrap_unit, wrap_unit)
            - 0.5 * wrap_unit
            < 0
            else idx1 + 1
        )
        if idx2 < 0:
            idx2 = selection.shape[0] - 1
        if idx2 >= selection.shape[0]:
            idx2 = 0
    else:
        idx1 = np.argmin(
            np.absolute(goal - selection),
        )
        idx2 = idx1 - 1 if goal < selection[idx1] else idx1 + 1
        if idx2 < 0 or idx2 >= selection.shape[0]:
            idx2 = idx1
    weights = _get_weights(goal, idx1, idx2, selection, wrap_unit)
    return Indices(idx1=idx1, idx2=idx2, w1=weights.w1, w2=weights.w2)


def get_sorted_indices(
    lon: float,
    lat: float,
    avail_lon: NDArray[np.float64],
    avail_lat: NDArray[np.float64],
    wrap_unit: float = 360.0,
) -> SortedIndices:
    """find distances of a lon/lat grid to a point and return sorted list of indices

    Parameters
    ----------
    lon : float
        longitude
    lat : float
        latitude
    avail_lon : NDArray[np.float64]
        array of available longitudes (must have same length as avail_lon)
    avail_lat : NDArray[np.float64]
        array of available latitudes (must have same length as avail_lat)
    wrap_unit : float, optional
        if goal/selection is a wrapping entity (e.g. angles) set this to the wrap value (e.g. 360),
        by default 360.0

    Returns
    -------
    SortedIndices
        sorted indices and distances in the array

    Raises
    ------
    ArrayShapeError
        if the shape of avail_lon is not equal to the shape of avail_lat
    """

    if avail_lon.shape != avail_lat.shape:
        msg = f"shapes of longitudes {avail_lon.shape} and lattiudes {avail_lat.shape} need to be equal"
        raise ArrayShapeError(msg)
    # some messed up thinking here
    distance = (
        wrap_around_zero(avail_lon - lon, wrap_unit) ** 2
        + wrap_around_zero(avail_lat - lat, wrap_unit) ** 2
    )
    sorted_idx = np.argsort(distance)
    return SortedIndices(indices=sorted_idx, distance=distance[sorted_idx])


def get_interpol(data: NDArray[np.float64], dist: NDArray[np.float64]) -> float:
    """get distance weighted sum of data

    Parameters
    ----------
    data : NDArray[np.float64]
        input data
    dist : NDArray[np.float64]
        distances (inverse weights)

    Returns
    -------
    float
        weighted sum of data
    """

    if np.any(dist == 0):
        w = np.zeros(dist.shape, dtype=float)
        w[dist == 0] = 1.0 / np.sum(dist == 0)
    else:
        w = 1.0 / dist
    w /= np.sum(w)
    return float(np.sum(data * w))


def wrap_around_zero(
    data: NDArray[np.float64], wrap_unit: float = 2 * np.pi
) -> NDArray[np.float64]:
    """Function to calculate the remainder of data such that this is centered around zero

    Parameters
    ----------
    data : NDArray[np.float64]
        input data
    wrap_unit : float, optional
        unit for wrapping, by default 2*np.pi

    Returns
    -------
    NDArray[np.float64]
        wrapped data
    """

    return np.remainder(data + 0.5 * wrap_unit, wrap_unit) - 0.5 * wrap_unit


def compute_index_and_weights(
    maparray: NDArray[np.float64], mapvalues: float | NDArray[np.float64]
) -> Indices:
    """helper function  to get indices and weights for interpolating tecmaps
    Only works for the non wrapping axes. So time and latitude


    Args:
        maparray (NDArray) : array to get indices in
        mapvalues (float | np.array ):  values to get indices for
    Returns:
        Tuple[np.array, np.array, np.array]: idx1,idx2 and weights for idx2,
                                             idx2 is always >= idx1

    """
    is_reverse = maparray[1] < maparray[0]
    idx1 = np.argmin(
        np.absolute(maparray[np.newaxis] - np.atleast_1d(mapvalues)[:, np.newaxis]),
        axis=1,
    )
    idx2 = idx1.copy()
    if not is_reverse:
        idx1[maparray[idx1] > mapvalues] -= 1
        idx2[maparray[idx2] < mapvalues] += 1
    else:
        idx1[maparray[idx1] < mapvalues] -= 1
        idx2[maparray[idx2] > mapvalues] += 1
    idx1[idx1 < 0] = 0
    idx2[idx2 < 0] = 0
    idx1[idx1 >= maparray.shape[0]] = maparray.shape[0] - 1
    idx2[idx2 >= maparray.shape[0]] = maparray.shape[0] - 1
    _steps = np.absolute(maparray[idx2] - maparray[idx1])
    weights = np.absolute(mapvalues - maparray[idx1])
    weights[_steps == 0] = 1.0
    weights[_steps != 0] = weights[_steps != 0] / _steps[_steps != 0]
    return Indices(idx1=idx1, idx2=idx2, w1=1 - weights, w2=weights)
