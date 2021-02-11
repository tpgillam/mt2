#!/usr/bin/env python

"""Tests for `mt2` package."""

import pytest

from mt2 import mt2

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

    computed_val = mt2(mVisA, pxA, pyA, mVisB, pxB, pyB, pxMiss, pyMiss, chiA, chiB)

    assert computed_val == pytest.approx(412.628)

def test_near_massless_lester():
    
    # This test is based on Fig 5 of https://arxiv.org/pdf/1411.4312.pdf
    
    mVisA = 0
    pxA = -42.017340486
    pyA = -146.365340528

    mVisB = 0.087252259
    pxB = -9.625614206
    pyB = 145.757295514

    pxMiss = -16.692279406
    pyMiss = -14.730240471

    chiA = 0
    chiB = 0

    computed_val = mt2(mVisA, pxA, pyA, mVisB, pxB, pyB, pxMiss, pyMiss, chiA, chiB)

    assert computed_val == pytest.approx(0.09719971)
