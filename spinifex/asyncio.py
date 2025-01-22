"""Async utility functions for Spinifex."""

from __future__ import annotations

import asyncio
import sys
from collections.abc import Coroutine
from functools import wraps
from typing import Callable, TypeVar

if sys.version_info < (3, 10):
    from typing_extensions import ParamSpec
else:
    from typing import ParamSpec

P = ParamSpec("P")
T = TypeVar("T")


def sync_wrapper(coro: Callable[P, Coroutine[None, None, T]]) -> Callable[P, T]:
    @wraps(coro)
    def wrapper(
        *args: P.args,
        **kwargs: P.kwargs,
    ) -> T:
        # First see if we are already in an event loop
        try:
            loop = asyncio.get_event_loop()
        except RuntimeError:
            # Nope - we can go ahead and create one
            return asyncio.run(coro(*args, **kwargs))
        else:
            # We are in an event loop, so we can't use run
            # Submit the coroutine to the event loop
            return asyncio.run_coroutine_threadsafe(
                coro(*args, **kwargs), loop
            ).result()

    # Keep the function docs correct
    wrapper.__doc__ = coro.__doc__
    return wrapper
