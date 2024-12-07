#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module for getting ionospheric Faraday rotation from external geomagnetic and ionospheric models"""

from __future__ import annotations

from importlib import metadata
from typing import Callable

__version__ = metadata.version("spinifex")
from spinifex.ionospheric import models


def get_rm(iono_model: Callable = models.IONEX):
    """Prints a nice message"""
    # IPP = geometry.getIPP()
    # iono = get_density_profile(IPP, time, iono_model)

    return iono_model


get_rm(models.IONEX)
