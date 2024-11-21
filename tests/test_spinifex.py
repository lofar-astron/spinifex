#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0

"""Testing of the main module"""

from __future__ import annotations

from unittest import TestCase

# from spinifex import get_rm


class TestSpinifex(TestCase):
    """Test Case of the Cool Module"""

    def test_get_rm(self):
        """Test that get_rm does not crash"""
        # get_rm()
        assert 2 + 2 == 4
