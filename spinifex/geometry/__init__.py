#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module for getting the Earth magnetic field"""

from __future__ import annotations

from spinifex.geometry.get_ipp import IPP, get_ipp_from_altaz, get_ipp_from_skycoord

__all__ = ["IPP", "get_ipp_from_altaz", "get_ipp_from_skycoord"]
