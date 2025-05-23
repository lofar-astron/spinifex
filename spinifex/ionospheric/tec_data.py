"""data objects and options for ionospheric models"""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Literal, NamedTuple

import astropy.units as u
import numpy as np
from numpy.typing import NDArray
from pydantic import Field

from spinifex.options import Options

SOLUTION = Literal["final", "rapid"]

CENTER_NAMES = {
    "cod",
    "esa",
    "igs",
    "jpl",
    "upc",
    "irt",
    "uqr",
}
DEFAULT_TIME_RESOLUTIONS: dict[str, u.Quantity] = {
    "cod": 1 * u.hour,
    "esa": 2 * u.hour,
    "igs": 2 * u.hour,
    "jpl": 2 * u.hour,
    "upc": 2 * u.hour,
    "irt": 2 * u.hour,
    "uqr": 15 * u.min,  # Chapman only provides 15 minute resolution right now
}

assert set(DEFAULT_TIME_RESOLUTIONS.keys()) == CENTER_NAMES, (
    "Time resolutions must be defined for all analysis centres"
)

NAME_SWITCH_WEEK = 2238  # GPS Week where the naming convention changed


class Servers(str, Enum):
    CDDIS = ("cddis", "https://cddis.nasa.gov/archive/gnss/products/ionex")
    CHAPMAN = ("chapman", "http://chapman.upc.es/tomion/rapid")
    IGSIONO = ("igsiono", "ftp://igs-final.man.olsztyn.pl")

    def __new__(cls, value: str, _url: str) -> Servers:
        obj = str.__new__(cls, value)
        obj._value_ = value
        obj._url = _url  # type: ignore[attr-defined]
        return obj

    @property
    def url(self) -> str:
        return str(self._url)  # type: ignore[attr-defined]

    @classmethod
    def get_url(cls, server: Servers) -> str:
        return str(server.url)  # Retrieve URL given an Enum value


class IonexOptions(Options):
    """Options for ionex model"""

    server: Servers = Field(
        Servers.CHAPMAN, description="Server to download ionex files from"
    )
    prefix: str = Field("uqr", description="Analysis centre prefix")
    url_stem: str | None = Field(None, description="URL stem")
    time_resolution: u.Quantity | None = Field(
        None, description="Time resolution for ionex files"
    )
    solution: SOLUTION = Field("final", description="Solution type")
    output_directory: Path | None = Field(
        None, description="Output directory for ionex files"
    )
    correct_uqrg_rms: bool = Field(
        True, description="Correct overestimated rms of uqr maps"
    )
    height: u.Quantity = Field(
        350 * u.km, description="altitude of single layer ionosphere"
    )
    remove_midnight_jumps: bool = Field(
        True,
        description="mitigate midnight jumps in the ionex files by inserting the data of the next day",
    )


class TomionOptions(Options):
    """Options for tomion model"""

    output_directory: Path | None = Field(
        None, description="Output directory for tomion files"
    )


class ElectronDensity(NamedTuple):
    """object containing interpolated electron density values and their estimated uncertainty"""

    electron_density: NDArray[np.float64]
    """electron density in TECU"""
    electron_density_error: NDArray[np.float64]
    """uncertainty in TECU"""
