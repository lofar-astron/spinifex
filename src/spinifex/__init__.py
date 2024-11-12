"""
Copyright (c) 2024 Maaijke Mevius, Alec Thomson. All rights reserved.

spinifex: Tool for calculating ionospheric Faraday rotation
"""

from __future__ import annotations

from ._version import version as __version__
from .spinifex import get_bla, get_blo

__all__ = ["__version__", "get_bla", "get_blo"]
