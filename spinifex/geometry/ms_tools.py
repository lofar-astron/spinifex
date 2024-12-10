#  Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Module for getting the Ionospheric Piercepoints"""

from __future__ import annotations

from typing import NamedTuple

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from numpy.typing import ArrayLike

from spinifex.geometry.get_ipp import IPP, get_ipp_from_skycoord


class MsMetaData(NamedTuple):
    """Metadata from a Measurement Set"""
    times: Time
    location: EarthLocation
    name: str
    source: SkyCoord


def get_ipp_from_ms(
    ms: str, height_array: ArrayLike, timestep: int = 1
) -> IPP:  # depends on casacore,
    # an ms has different stations, you optionally want to return an ipp object for every station
    msmetadata = get_metadata_from_ms(ms, timestep)
    return get_ipp_from_skycoord(
        loc=msmetadata.location,
        times=msmetadata.times,
        source=msmetadata.source,
        height_array=height_array,
    )


def get_metadata_from_ms(ms: str, timestep: int) -> MsMetaData:
    """Dummy: needs to be implemented"""
    times = Time.now() + np.arange(10) * u.min
    return MsMetaData(
        times=times[::timestep],
        location=EarthLocation(),
        name=ms,
        source=SkyCoord.from_name("Cas A"),
    )
