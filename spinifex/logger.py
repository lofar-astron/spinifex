"""Logging utilities for spinifex."""

from __future__ import annotations

import io
import logging

logging.captureWarnings(True)


formatter = logging.Formatter(
    fmt="%(levelname)s %(module)s - %(funcName)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


class TqdmToLogger(io.StringIO):
    """
    Output stream for TQDM which will output to logger module instead of
    the StdOut.
    """

    logger = None
    level = None
    buf = ""

    def __init__(self, logger: logging.Logger | None, level: int | None = None) -> None:
        super().__init__()
        self.logger = logger
        self.level = level or logging.INFO

    def write(self, buf: str) -> int:
        self.buf = buf.strip("\r\n\t ")
        return len(buf)

    def flush(self) -> None:
        if self.logger is not None and isinstance(self.level, int):
            self.logger.log(self.level, self.buf)


def set_verbosity(verbosity: int) -> None:
    """Set the logger verbosity.

    Parameters
    ----------
    verbosity : int
        Verbosity level
    """
    if verbosity == 0:
        level = logging.WARNING
    elif verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG
    else:
        level = logging.CRITICAL

    logging.getLogger().setLevel(level)
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    ch.setLevel(level)
    logging.getLogger().addHandler(ch)


logger = logging.getLogger("spinifex")

ch = logging.StreamHandler()
ch.setFormatter(formatter)
ch.setLevel(logging.INFO)
logger.addHandler(ch)
logger.setLevel(logging.INFO)
