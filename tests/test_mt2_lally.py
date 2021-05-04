"""Tests for the variant of MT2 by Colin Lally."""

import numpy
import pytest

from .common import mt2_lally, mt2_lester


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


@pytest.mark.skip(reason="Currently failing due to inconsistencies")
def test_fuzz():
    batch_size = 100
    num_tests = 1000

    numpy.random.seed(42)

    def _random_batch(min_, max_):
        return numpy.random.uniform(min_, max_, (batch_size,))

    for _ in range(num_tests):
        m_vis_1 = _random_batch(0, 100)
        px_vis_1 = _random_batch(-100, 100)
        py_vis_1 = _random_batch(-100, 100)
        m_vis_2 = _random_batch(0, 100)
        px_vis_2 = _random_batch(-100, 100)
        py_vis_2 = _random_batch(-100, 100)
        px_miss = _random_batch(-100, 100)
        py_miss = _random_batch(-100, 100)
        m_invis_1 = _random_batch(0, 100)
        m_invis_2 = _random_batch(0, 100)

        args = (
            m_vis_1,
            px_vis_1,
            py_vis_1,
            m_vis_2,
            px_vis_2,
            py_vis_2,
            px_miss,
            py_miss,
            m_invis_1,
            m_invis_2,
        )

        result_lester = mt2_lester(*args)
        result_lally = mt2_lally(*args)

        numpy.testing.assert_allclose(result_lester, result_lally, rtol=1e-12)
