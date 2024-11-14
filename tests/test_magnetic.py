#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Testing of magnetic"""
from unittest import TestCase

from spinifex.magnetic import get_magnetic_field


class TestMagnetic(TestCase):
    """Testing of magnetic"""

    def test_get_magnetic_field(self):
        """Test that get_magnetic does not crash"""
        get_magnetic_field()
        self.assertEqual(2 + 2, 4)
