#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module for getting ionospheric models"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import astropy.units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time

from spinifex.ionospheric.iono_models import (
    get_density_ionex_iri,
    get_density_ionex_single_layer,
)

__all__ = ["get_ionosphere"]


@dataclass
class IonosphericModels:
    """Names space for different ionospheric models. An ionospheric model should be
    a callable get_density"""

    ionex: Callable
    ionex_iri: Callable


models = IonosphericModels(
    ionex=get_density_ionex_single_layer,
    ionex_iri=get_density_ionex_iri,
)


def get_ionosphere(
    loc: EarthLocation,
    time: Time,
    iono_model: Callable = models.ionex,
) -> u.Quantity[u.m**-3]:
    """Prints a nice message"""

    return iono_model(loc=loc, times=time)
