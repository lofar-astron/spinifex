from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from numpy.typing import NDArray

from spinifex.geometry import IPP, get_ipp_from_skycoord

try:
    from casacore.tables import table
except ImportError as e:
    MSG = "casacore is not installed! To operate on MeasurementSets, install spinifex[casa]."
    raise ImportError(MSG) from e


class MsMetaData(NamedTuple):
    """Metadata from a Measurement Set"""

    times: Time
    location: EarthLocation
    name: str
    source: SkyCoord


def get_ipp_from_ms(
    ms: str, height_array: NDArray[np.float64], timestep: int = 1
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


def get_columns_from_ms(ms_path: Path) -> list[str]:
    """Get the columns from a MeasurementSet"""
    with table(ms_path.as_posix()) as tab:
        return list(tab.colnames())
