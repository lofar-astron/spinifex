"""Utility functions for Spinifex."""

from __future__ import annotations

import asyncio
import inspect
from functools import wraps
from typing import Callable


def sync_wrapper(coro: Callable) -> Callable:
    """Creates a synchronous wrapper for an asynchronous coroutine."""

    @wraps(coro)
    def wrapper(*args, **kwargs):
        # Call the coroutine in an event loop
        return asyncio.run(coro(*args, **kwargs))

    # Dynamically copy the signature
    wrapper.__signature__ = inspect.signature(coro)
    return wrapper
