"""Some useful multipurpose functions to interpolate on grids"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np


class Indices(NamedTuple):
    """Indices of the closest two points in a possibly wrapping selection and the inverse distance weights"""

    idx1: np.ndarray[int]
    """Index of the first closest point"""
    idx2: np.ndarray[int]
    """Index of the second closest point"""
    w1: np.ndarray[float]
    """Weight of the first closest point"""
    w2: np.ndarray[float]
    """Weight of the second closest point"""


class SortedIndices(NamedTuple):
    """Indices of the closest two points in a possibly wrapping selection"""

    indices: np.ndarray
    """Index of the first closest point"""
    distance: np.ndarray
    """Index of the second closest point"""


class Weights(NamedTuple):
    """Weights of the closest two points in a possibly wrapping selection"""

    w1: np.ndarray
    """Weight of the first closest point"""
    w2: np.ndarray
    """Weight of the second closest point"""


def get_indices_axis(
    goal: np.ndarray, selection: np.ndarray, wrap_unit: float = 0
) -> Indices:
    """get indices of the closest two points in a possibly wrapping selection for an array
    of goals"""
    if wrap_unit > 0:
        idx1 = np.argmin(
            np.absolute(
                np.remainder(
                    goal[..., np.newaxis] - selection + 0.5 * wrap_unit, wrap_unit
                )
                - 0.5 * wrap_unit
            ),
            axis=-1,
        )
    else:
        idx1 = np.argmin(np.absolute(goal[..., np.newaxis] - selection), axis=-1)
    if np.isscalar(idx1):
        if goal < selection[idx1]:
            idx2 = idx1
            idx1 -= 1
            if idx1 < 0:
                idx1 = selection.shape[0] if wrap_unit > 0 else 0
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
    goal: np.ndarray,
    index1: np.ndarray,
    index2: np.ndarray,
    selection: np.ndarray,
    wrap_unit: float = 0,
) -> Weights:
    """Calculate weights based on distance of goal to selection

    Parameters
    ----------
    goal : np.ndarray
        array of points to get weights for
    index1 : np.ndarray
        indices in selection for goals (index1, index2) per goal
    index2 : np.ndarray
        indices in selection for goals (index1, index2) per goal
    selection : np.ndarray
        array to select from
    wrap_unit : float, optional
        if goal/selection is a wrapable (e.g. angle) set this unit (e.g. 360), by default 0

    Returns
    -------
    namedtuple
        tuple with weights
    """
    distance1 = np.abs(wrap_around_zero(selection[index1] - goal, wrap_unit))
    distance2 = np.abs(wrap_around_zero(selection[index2] - goal, wrap_unit))
    sumdist = distance1 + distance2
    return Weights(w1=1 - distance1 / sumdist, w2=1 - distance2 / sumdist)


def get_indices(goal: float, selection: np.ndarray, wrap_unit: float = 0) -> Indices:
    """find the indices of the closest two points in a possibly wrapping array selection

    Parameters
    ----------
    goal : float
        location of point
    selection : np.ndarray
        array of points
    wrap_unit : float, optional
        if goal/selection is a wrapping entity (e.g. angles) set this to the wrap value (e.g. 360), by default 0

    Returns
    -------
    Indices:
        sorted list of index1 and index2
    """

    if wrap_unit > 0:
        idx1 = np.argmin(
            np.absolute(
                np.remainder(goal - selection + 0.5 * wrap_unit, wrap_unit)
                - 0.5 * wrap_unit
            )
        )
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
    avail_lon: np.ndarray,
    avail_lat: np.ndarray,
    wrap_unit: float = 360.0,
) -> SortedIndices:
    """find distances of a lon/lat grid to a point and return sorted list of indices"""
    distance = (
        wrap_around_zero(avail_lon - lon, wrap_unit) ** 2
        + wrap_around_zero(avail_lat - lat, wrap_unit) ** 2
    )
    sorted_idx = np.argsort(distance)
    return SortedIndices(indices=sorted_idx, distance=distance[sorted_idx])


def get_interpol(data: np.ndarray, dist: np.ndarray):
    """get distance weighted sum of data"""
    w = 1.0 / dist
    w /= np.sum(w)
    return np.sum(data * w)


def wrap_around_zero(data: np.ndarray, wrap_unit: float = 2 * np.pi):
    """Function to calculate the remainder of data such that this is centered around zero"""
    return np.remainder(data + 0.5 * wrap_unit, wrap_unit) - 0.5 * wrap_unit


def _compute_index_and_weights(maparray: np.ndarray, mapvalues: np.ndarray) -> Indices:
    """helper function  to get indices and weights for interpolating tecmaps


    Args:
        maparray (np.ndarray) : array to get indices in
        mapvalues (Union[float,np.array]) :  values to get indices for
    Returns:
        Tuple[np.array, np.array, np.array]: idx1,idx2 and weights for idx2,
                                             idx2 is always >= idx1

    """
    is_reverse = maparray[1] < maparray[0]
    idx1 = np.argmin(
        np.absolute(maparray[np.newaxis] - mapvalues[:, np.newaxis]), axis=1
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
