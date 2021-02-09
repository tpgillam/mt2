#!/usr/bin/env python

"""Tests for `mt2` package."""

import pytest

from mt2 import get_mT2


def test_basic():
    with pytest.raises(NotImplementedError):
        get_mT2(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
