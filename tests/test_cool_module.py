#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Testing of the Cool Module"""
from unittest import TestCase

from spinifex.cool_module import greeter


class TestCoolModule(TestCase):
    """Test Case of the Cool Module"""

    def test_greeter(self):
        """Testing that the greeter does not crash"""
        greeter()
        self.assertEqual(2 + 2, 4)
