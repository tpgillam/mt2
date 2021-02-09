#!/usr/bin/env python

"""Tests for `mt2` package."""

import pytest

from mt2 import get_mT2

def test_simple_lester():
    mVisA = 100
    pxA = 410
    pyA = 20

    mVisB = 150
    pxB = -210
    pyB = -300

    pxMiss = -200
    pyMiss = 280

    chiA = 100
    chiB = 100

    computed_val = get_mT2(mVisA, pxA, pyA, mVisB, pxB, pyB, pxMiss, pyMiss, chiA, chiB)

    assert computed_val == pytest.approx(412.628)
