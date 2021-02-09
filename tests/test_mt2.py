#!/usr/bin/env python

"""Tests for `mt2` package."""

import pytest

from mt2 import get_mt2

def test_basic():
    num_times = 100

    rep = lambda x: [x] * num_times

    pxA = 410
    pyA = rep(20)
    mVisA = rep(100)
    pxB = rep(-210)
    pyB = rep(-300)
    mVisB = rep(150)
    pxMiss = rep(-200)
    pyMiss = rep(280)
    chiA = rep(100)
    chiB = rep(100)
    computed_val = get_mt2(mVisA, pxA, pyA, mVisB, pxB, pyB, pxMiss, pyMiss, chiA, chiB, 0, False)
    expected_val = 412.628
    print("Testing Tom's MT2 PiPI thing.")
    print("Expected mT2 = " + str(expected_val))
    print("Computed mT2 = " + str(computed_val))

test_basic()
