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

class TestGaussLegendreFirstOrder(unittest.TestCase):

    def setUp(self):
        self.f = lambda x: x + 6
    
    def test_npts(self):
        expected = 80
        self.assertAlmostEqual(integrate_gauss(self.f, lims = [0, 8], npts = 1), expected)

class TestGaussLegendreSecondOrder(unittest.TestCase):

    def setUp(self):
        self.f = lambda x: x ** 2 + 2
    
    def test_npts(self):
        expected = 155 / 3
        self.assertAlmostEqual(integrate_gauss(self.f, lims = [0, 5], npts = 2), expected)

class TestGaussLegendreThirdOrder(unittest.TestCase):

    def setUp(self):
        self.f = lambda x: x ** 3 + 1
    
    def test_npts(self):
        expected = 161.25
        self.assertAlmostEqual(integrate_gauss(self.f, lims = [0, 5], npts = 2), expected)

class TestGaussLegendreFourthOrder(unittest.TestCase):

    def setUp(self):
        self.f = lambda x: 3 * x ** 4 + 2
    
    def test_npts(self):
        expected = 23.2
        self.assertAlmostEqual(integrate_gauss(self.f, lims = [0, 2], npts = 3), expected)

class TestGaussLegendreFifthOrder(unittest.TestCase):

    def setUp(self):
        self.f = lambda x: 2 * x ** 5
    
    def test_npts(self):
        expected = 243
        self.assertAlmostEqual(integrate_gauss(self.f, lims = [0, 3], npts = 3), expected)

class TestGaussLegendreSixthOrder(unittest.TestCase):

    def setUp(self):
        self.f = lambda x: x ** 6 + 3
    
    def test_npts(self):
        expected = 2250 / 7
        self.assertAlmostEqual(integrate_gauss(self.f, lims = [0, 3], npts = 4), expected)

class TestGaussLegendreSeventhOrder(unittest.TestCase):

    def setUp(self):
        self.f = lambda x: x ** 7 + 2
    
    def test_npts(self):
        expected = 36
        self.assertAlmostEqual(integrate_gauss(self.f, lims = [0, 2], npts = 4), expected)

class TestGaussLegendreEighthOrder(unittest.TestCase):

    def setUp(self):
        self.f = lambda x: x ** 8 + 2
    
    def test_npts(self):
        expected = 19 / 9
        self.assertAlmostEqual(integrate_gauss(self.f, lims = [0, 1], npts = 5), expected)

class TestGaussLegendreNinthOrder(unittest.TestCase):

    def setUp(self):
        self.f = lambda x: x ** 9 + 8
    
    def test_npts(self):
        expected = 8.1
        self.assertAlmostEqual(integrate_gauss(self.f, lims = [0, 1], npts = 5), expected)

