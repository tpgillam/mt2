"""Tests for the variant of MT2 by Rupert Tombs."""

import math
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
        # Any negative mass is un-physical. We clip negative values to zero,
        # such that all non-positive mass arguments are equivalent.
        # And by equivalent, I mean equivalent: strict floating point equality.
        negative = mt2_tombs(1, 2, 3, 4, 5, 6, 7, 8, -90, -100)
        zero = mt2_tombs(1, 2, 3, 4, 5, 6, 7, 8, 0, 0)
        self.assertEqual(negative, zero)
        negative_zero = mt2_tombs(1, 2, 3, 4, 5, 6, 7, 8, -0.0, -0.0)
        self.assertEqual(negative_zero, zero)

    def test_zero_mass(self):
        # Zero masses are perfectly fine.
        # Thanks to Sebastian Rutherford Colmenares for this example in:
        #     https://github.com/tpgillam/mt2/issues/75#issuecomment-2656551786

        # MT2 must be positive and finite for all-zero masses.
        zero = mt2_tombs(
            0.0,  # Visible 1 mass
            -30500.0,
            34500.0,
            0.0,  # Visible 2 mass
            -29100.0,
            -55400.0,
            58900.0,
            20300.0,
            0.0,  # Invisible 1 mass
            0.0,  # Invisible 2 mass
        )
        self.assertGreater(zero, 0)
        self.assertTrue(math.isfinite(zero))

        # MT2 must change continuously between from zero mass to small mass.
        small = mt2_tombs(
            0.5,  # Visible 1 mass
            -30500.0,
            34500.0,
            0.5,  # Visible 2 mass
            -29100.0,
            -55400.0,
            58900.0,
            20300.0,
            0.5,  # Invisible 1 mass
            0.5,  # Invisible 2 mass
        )
        self.assertAlmostEqual(zero, small, delta=1e-3)
        self.assertGreater(small, zero)
