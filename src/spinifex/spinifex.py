"""All the stuff that is not in submodules"""

from __future__ import annotations

import numpy as np


def get_bla(my_input: float) -> float:
    """Compute three times input"""
    return my_input * 3


def get_blo(my_input: np.ndarray) -> np.ndarray:
    """Compute three times input"""
    return my_input * 3
