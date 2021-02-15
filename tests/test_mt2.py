"""Tests for `mt2` package."""
import math
import random

import pytest

from mt2 import mt2


def test_simple_example():
    computed_val = mt2(100, 410, 20, 150, -210, -300, -200, 280, 100, 100)
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

    computed_val = mt2(
        m_vis_a, px_a, py_a, m_vis_b, px_b, py_b, px_miss, py_miss, chi_a, chi_b
    )
    assert computed_val == pytest.approx(0.09719971)


def test_collinear_endpoint_cases():
    random.seed(0)  # If the test fails, we want to be able to repeat the test!
    for i in range(10000):
        m_vis_a = random.uniform(0, 10)
        m_vis_b = random.uniform(0, 10)
        m_invis_a = random.uniform(0, 10)
        m_invis_b = random.uniform(0, 10)

        # Addition of random positive number prevents possibility of division
        # by zero at the locations below marked "MOO"
        m_parent = max(m_vis_a + m_invis_a, m_vis_b + m_invis_b) + random.uniform(
            0.1, 10
        )
        p_parent_a = random.uniform(0, 10)
        p_parent_b = random.uniform(0, 10)
        e_parent_a = math.sqrt(p_parent_a ** 2 + m_parent ** 2)
        e_parent_b = math.sqrt(p_parent_b ** 2 + m_parent ** 2)
        beta_a = p_parent_a / e_parent_a  # MOO
        beta_b = p_parent_b / e_parent_b  # MOO
        gamma_a = 1.0 / math.sqrt(1 - beta_a ** 2)
        gamma_b = 1.0 / math.sqrt(1 - beta_b ** 2)
        pA = math.sqrt(
            (m_parent - m_vis_a - m_invis_a)
            * (m_parent + m_vis_a - m_invis_a)
            * (m_parent - m_vis_a + m_invis_a)
            * (m_parent + m_vis_a + m_invis_a)
        ) / (2 * m_parent)
        pB = math.sqrt(
            (m_parent - m_vis_b - m_invis_b)
            * (m_parent + m_vis_b - m_invis_b)
            * (m_parent - m_vis_b + m_invis_b)
            * (m_parent + m_vis_b + m_invis_b)
        ) / (2 * m_parent)
        p_vis_a_boosted = gamma_a * (beta_a * math.sqrt(m_vis_a ** 2 + pA ** 2) + pA)
        p_vis_b_boosted = gamma_b * (beta_b * math.sqrt(m_vis_b ** 2 + pB ** 2) + pB)
        p_invis_a_boosted = gamma_a * (
            beta_a * math.sqrt(m_invis_a ** 2 + pA ** 2) - pA
        )
        p_invis_b_boosted = gamma_b * (
            beta_b * math.sqrt(m_invis_b ** 2 + pB ** 2) - pB
        )
        p_miss = p_invis_a_boosted + p_invis_b_boosted
        theta = random.uniform(0, math.tau)
        c = math.cos(theta)
        s = math.sin(theta)
        px_miss, py_miss = p_miss * c, p_miss * s
        ax, ay = p_vis_a_boosted * c, p_vis_a_boosted * s
        bx, by = p_vis_b_boosted * c, p_vis_b_boosted * s
        val = mt2(
            m_vis_a, ax, ay, m_vis_b, bx, by, px_miss, py_miss, m_invis_a, m_invis_b
        )

        # passes with rel=1e-12 but sporadically fails with rel=1e-13
        assert val == pytest.approx(m_parent, rel=1e-12), (
            f"WARNING! Expected {m_parent} from collinear event but instead got "
            f"{val} for mt2({m_vis_a},{ax},{ay}, {m_vis_b},{bx},{by}, "
            f"{px_miss},{py_miss}, {m_invis_a},{m_invis_b}) in test case {i}."
        )
