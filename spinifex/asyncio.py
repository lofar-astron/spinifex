"""Async utility functions for Spinifex."""

from __future__ import annotations

import asyncio
import sys
from collections.abc import Coroutine
from functools import wraps
from typing import Callable, TypeVar

import nest_asyncio

if sys.version_info < (3, 10):
    from typing_extensions import ParamSpec
else:
    from typing import ParamSpec

P = ParamSpec("P")
T = TypeVar("T")


def sync_wrapper(coro: Callable[P, Coroutine[None, None, T]]) -> Callable[P, T]:
    nest_asyncio.apply()

    @wraps(coro)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
        try:
            loop = asyncio.get_running_loop()
        except RuntimeError:
            # No running loop: safe to use asyncio.run
            return asyncio.run(coro(*args, **kwargs))
        else:
            return loop.run_until_complete(coro(*args, **kwargs))

    return wrapper
