#!/usr/bin/env python

"""Tests for `mt2` package."""

import pytest

from mt2 import get_mt2


def test_basic():
    with pytest.raises(NotImplementedError):
        get_mt2(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
