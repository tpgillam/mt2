"""Tests for the variant of MT2 by Colin Lally."""

import pytest

from mt2 import mt2_lally


def test_simple_example():
    computed_val = mt2_lally(100, 410, 20, 150, -210, -300, -200, 280, 100, 100)
    assert computed_val == pytest.approx(412.628)


def test_near_massless():
    # This test is based on Fig 5 of https://arxiv.org/pdf/1411.4312.pdf
    m_vis_a = 0
    px_a = -42.017340486
    py_a = -146.365340528

    m_vis_b = 0.087252259
    px_b = -9.625614206
    py_b = 145.757295514

    px_miss = -16.692279406
    py_miss = -14.730240471

    chi_a = 0
    chi_b = 0

    computed_val = mt2_lally(
        m_vis_a, px_a, py_a, m_vis_b, px_b, py_b, px_miss, py_miss, chi_a, chi_b
    )
    assert computed_val == pytest.approx(0.09719971)
