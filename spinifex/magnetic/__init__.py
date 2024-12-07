#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module for getting the Earth magnetic field"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from spinifex.magnetic.magnetic_models import get_ppigrf_magnetic_field


@dataclass
class MagneticModels:
    ppigrf: Callable


models = MagneticModels(ppigrf=get_ppigrf_magnetic_field)
