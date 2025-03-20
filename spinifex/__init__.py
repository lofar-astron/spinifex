#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module for getting ionospheric Faraday rotation from external geomagnetic and ionospheric ionospheric_models"""

from __future__ import annotations

from importlib import metadata

__version__ = metadata.version("spinifex")
from spinifex import options
from spinifex.ionospheric import ModelDensityFunction, ionospheric_models

__all__ = ["ModelDensityFunction", "ionospheric_models", "options"]
