import unittest
import numpy as np
import scipy as sc

from goph420_lab01.integration import integrate_gauss

class TestGaussLegendreInvalidInitializers(unittest.TestCase):

    def test_callable_object(self):
        with self.assertRaises(TypeError):
            integrate_gauss(f := 13.0, lims = [1, 2])

    def test_length_lims(self):
        with self.assertRaises(ValueError):
            integrate_gauss(f = lambda f: np.sin(f), lims = [1])

    def test_dtype_lims(self):
        with self.assertRaises(ValueError):
            integrate_gauss(f = lambda f: np.sin(f), lims = ["one", "two"])

    def test_npts(self):
        with self.assertRaises(ValueError):
            integrate_gauss(f = lambda f: np.sin(f), lims = [1, 2], npts = 7)

class TestGaussLegendreFromPolynnomials(unittest.TestCase):

    def setUp(self):
        pass

class TestGaussLegendreAgainstScipy(unittest.TestCase):

    def setUp(self):
        pass
