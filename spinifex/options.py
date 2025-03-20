from __future__ import annotations

from pydantic import BaseModel, ConfigDict


class Options(BaseModel):
    """Base model options class (immutable by default)."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True, extra="forbid")
