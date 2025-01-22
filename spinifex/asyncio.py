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
    # Apply the `nest_asyncio` magic to allow nested event loops
    # Need for running async functions in a Jupyter notebook
    nest_asyncio.apply()

    @wraps(coro)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
        try:
            # Check if we are in an event loop
            loop = asyncio.get_running_loop()
        except RuntimeError:
            # Nope, we are not in an event loop
            # Create a new event loop with `asyncio.run`
            return asyncio.run(coro(*args, **kwargs))
        # We are in an event loop
        # Use the current event loop to run the coroutine
        # This is useful when we are in a Jupyter notebook
        return loop.run_until_complete(coro(*args, **kwargs))

    return wrapper
