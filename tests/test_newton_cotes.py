import unittest
import numpy as np
import scipy as sc

from goph420_lab01.integration import integrate_newton

class TestNewtonCotesInvalidInitializers(unittest.TestCase):

    def test_string_flag(self):
        with self.assertRaises(ValueError):
            integrate_newton(x := np.linspace(1, 10), f = np.exp(x), alg = 'quad')
    
    def test_length_inputs(self):
        with self.assertRaises(ValueError):
            integrate_newton(x = np.linspace(1, 5, 3), f = np.linspace(1, 10, 10))

    def test_dimensions_inputs(self):
        with self.assertRaises(ValueError):
            integrate_newton(x = np.linspace(0, 10), f = np.array([[1, 2, 3], [1, 2, 3]]))
            integrate_newton(x = np.array([[1, 2, 3], [1, 2, 3]]), f = np.array([[1, 2, 3], [1, 2, 3]]))
class TestTrapRuleLinear(unittest.TestCase):
    # f(x) = x
    def setUp(self):
        self.x = [0, 1, 2]
        self.f = [0, 1, 2]

    def test_value(self):
        expected = sc.integrate.trapezoid(self.f, self.x)
        self.assertAlmostEqual(integrate_newton(self.x, self.f, alg = 'Trap '), expected)

class TestSimpRuleQuadOdd(unittest.TestCase):

    def setUp(self):
        self.x = np.linspace(0, 11, 9)
        self.f = 2 * (self.x ** 2) + 5
    
    def test_value(self):
        expected = sc.integrate.simpson(self.f, self.x)
        self.assertAlmostEqual(integrate_newton(self.x, self.f, alg = 'Simp'), expected)
    

class TestSimpRuleQuadEven(unittest.TestCase):

    def setUp(self):
        self.x = np.linspace(0, 4, 4)
        self.f = 4 * (self.x ** 2) + 3
    
    def test_value(self):
        expected = sc.integrate.simpson(self.f, self.x)
        self.assertAlmostEqual(integrate_newton(self.x, self.f, alg = 'Simp'), expected)

if __name__ == "__main__":
    unittest.main()