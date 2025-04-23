"""data objects and options for ionospheric models"""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import numpy as np
from numpy.typing import NDArray
from pydantic import Field

from spinifex.options import Options


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
