#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Testing of ionospheric"""
from unittest import TestCase

from spinifex.ionospheric import get_ionosphere


class TestIonospheric(TestCase):
    """Testing of magnetic"""

    def test_get_ionosphere(self):
        """Test that get_ionosphere does not crash"""
        get_ionosphere()
        self.assertEqual(2 + 2, 4)
