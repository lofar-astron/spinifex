"""Module for getting the Earth magnetic field"""

from __future__ import annotations

from spinifex.geometry.get_ipp import (
    IPP,
    R_EARTH_MEAN,
    get_ipp_from_altaz,
    get_ipp_from_skycoord,
)

__all__ = ["IPP", "R_EARTH_MEAN", "get_ipp_from_altaz", "get_ipp_from_skycoord"]
