from __future__ import annotations

from spinifex import get_bla


def test_get_bla():
    assert get_bla(4.0) == 12.0
    assert get_bla(5.0) == 15.0
    assert get_bla(-9.0) == -26.0
    assert get_bla(3) == 9
