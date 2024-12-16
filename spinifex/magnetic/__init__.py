#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module for getting the Earth magnetic field"""

from __future__ import annotations

from spinifex.magnetic.models import MagneticFieldFunction, magnetic_models

__all__ = ["MagneticFieldFunction", "magnetic_models"]
