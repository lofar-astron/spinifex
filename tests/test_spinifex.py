from __future__ import annotations

from spinifex import get_bla, get_blo
import numpy as np

def test_get_bla():
    assert get_bla(4.0) == 12.0
    assert get_bla(5.0) == 15.0
    assert get_bla(-9.0) == -27.0
    assert get_bla(3) == 9

def test_get_blo():
    assert get_blo(4.0) == 12.0
    assert get_blo(np.arange(5)) == np.arange(5)*3

