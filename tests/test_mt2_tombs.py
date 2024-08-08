"""Tests for the variant of MT2 by Rupert Tombs."""

import unittest

import numpy

from tests.common import mt2_lester, mt2_tombs


class TestTombs(unittest.TestCase):
    def test_simple_example(self):
        computed_val = mt2_tombs(100, 410, 20, 150, -210, -300, -200, 280, 100, 100)
        self.assertAlmostEqual(computed_val, 412.627668458219)

    def test_near_massless(self):
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

        computed_val = mt2_tombs(
            m_vis_a, px_a, py_a, m_vis_b, px_b, py_b, px_miss, py_miss, chi_a, chi_b
        )
        self.assertAlmostEqual(computed_val, 0.09719971)

    def test_fuzz(self):
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
            result_tombs = mt2_tombs(*args)

            numpy.testing.assert_allclose(result_lester, result_tombs, rtol=1e-12)

    def test_scale_invariance(self):
        example_args = numpy.array((100, 410, 20, 150, -210, -300, -200, 280, 100, 100))
        example_val = mt2_tombs(*example_args)

        # mt2 scales with its arguments; check over some orders of magnitude.
        for i in range(-100, 100, 10):
            scale = 10.0**i
            with numpy.errstate(over="ignore"):
                # Suppress overflow warnings when performing the evaluation; we're happy
                # so long as we match approximately in the test below.
                computed_val = mt2_tombs(*(example_args * scale))
            numpy.testing.assert_allclose(computed_val, example_val * scale)

    def test_negative_masses(self):
        # Any negative mass is unphysical.
        # These arguments use negative masses to make both initial bounds negative.
        # Check that the result is neither positive nor an infinite loop.
        computed_val = mt2_tombs(1, 2, 3, 4, 5, 6, 7, 8, -90, -100)
        self.assertLessEqual(computed_val, 0)
