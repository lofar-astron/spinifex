"""data objects for ionospheric data"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from numpy.typing import NDArray


class ElectronDensity(NamedTuple):
    """object containing interpolated electron density values and their estimated uncertainty"""

    electron_density: NDArray[np.float64]
    """electron density in TECU"""
    electron_density_error: NDArray[np.float64]
    """uncertainty in TECU"""
