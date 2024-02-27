import unittest
import numpy as np

from goph420_lab01.integration import integrate_newton

class TestNewtonCotesInvalidInitializers(unittest.TestCase):

    def test_string_flag(self):
        with self.assertRaises(ValueError):
            integrate_newton(x := np.arange(1, 10), f = np.exp(x), alg = 'quad')
    
    def test_length_inputs(self):
        with self.assertRaises(ValueError):
            integrate_newton(x = np.arange(1, 5), f = np.arange(1, 10))

    def test_dimensions_inputs(self):
        with self.assertRaises(ValueError):
            integrate_newton(x = np.arange(0, 10), f = np.array([[1, 2, 3], [1, 2, 3]]))
            integrate_newton(x = np.array([[1, 2, 3], [1, 2, 3]]), f = np.array([[1, 2, 3], [1, 2, 3]]))
class TestTrapRuleLinear(unittest.TestCase):
    # f(x) = x
    def setUp(self):
        self.x = [0, 0.5, 1, 1.5, 2]
        self.f = [0, 0.5, 1, 1.5, 2]

    def test_value(self):
        expected = 2.0
        self.assertAlmostEqual(integrate_newton(self.x, self.f, alg = 'Trap '), expected, delta=1e-15)

class TestSimpRuleQuadOdd(unittest.TestCase):

    def setUp(self):
        self.x = np.arange(1, 10)
        self.f = self.x ** 4 + 2
    
    def test_value(self):
        expected = 1064 / 5
        self.assertAlmostEqual(integrate_newton(self.x, self.f, alg = 'Simp'), expected, delta=1e-15)

if __name__ == "__main__":
    unittest.main()