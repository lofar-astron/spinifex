from __future__ import annotations

try:
    from casacore.tables import table
except ImportError as e:
    MSG = "casacore is not installed! To operate on MeasurementSets, install spinifex[casa]."
    raise ImportError(MSG) from e

from pathlib import Path


def get_columns_from_ms(ms_path: Path) -> list[str]:
    """Get the columns from a MeasurementSet"""
    with table(ms_path.as_posix()) as tab:
        return list(tab.colnames())
