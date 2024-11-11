from __future__ import annotations

import importlib.metadata

import spinifex as m


def test_version():
    assert importlib.metadata.version("spinifex") == m.__version__
