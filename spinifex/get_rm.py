from __future__ import annotations

from typing import NamedTuple

import numpy as np
from astropy.time import Time


class RM(NamedTuple):
    """object with all rotation measures"""

    rm: np.ndarray[float]
    """rotation measures"""
    times: Time
    """time axis"""


def get_rm() -> RM:
    raise NotImplementedError
