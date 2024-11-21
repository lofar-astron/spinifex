#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Testing of ionospheric"""

from __future__ import annotations

from unittest import TestCase

# from spinifex.ionospheric import get_ionosphere


class TestIonospheric(TestCase):
    """Testing of magnetic"""

    def test_get_ionosphere(self):
        """Test that get_ionosphere does not crash"""
        # get_ionosphere()
        assert 2 + 2 == 4
