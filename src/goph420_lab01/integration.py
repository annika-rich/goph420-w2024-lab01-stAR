# Utility Functions for Numerical Integration

import numpy as np

def integrate_newton(x, f, alg = "trap"):
    """Performs numerical integration of discrete data using Newton-Cotes rules.

    Parameters
    ----------
    x:  array-like
        contains the coordinates and values of the sample points

    f:  array_like
        contains the coordinates and values of the sample points
        must be the same shape as x

    alg:    optional string flag
            indicates integration rule used in implementation
            acceptable inputs are trap or simp


    Returns
    -------
    float, float-like
        Provides the integral estimate

    Raises
    ------
    ValueError:
        If the string flag contains a str other than trap or simp
        If the dimensions of x and f are incompatible
    """
    # make sure inputs are arrays
    x = np.array(x, dtype = float)
    f = np.array(f, dtype = float)

    # check that x and f have the same length
    if (m := len(x)) != (n := len(f)):
        raise ValueError(f"The length of x ({m}) is not equal to the length of f ({n}, x and f must be the same shape")

    # check dimensions of inputs
    mdim = len(x.shape)
    ndim = len(f.shape)
    if (mdim := len(x.shape)) != (ndim := len(f.shape)) or (mdim not in [1]):
        raise ValueError(f"x has {mdim} dimensions and f has {ndim} dimensions, dimensions must be equal and 1D")

    # check string flag
    alg = alg.strip().lower()
    if alg not in ["trap", "simp"]:
        raise ValueError(f"The algorithm string flag, contains a string ({alg}) other than 'trap' or 'simp'.")

    if alg == "trap":
        def trap():
            pass

    # may need to differentiate based on sample points...
    if alg == "simp":
        def simp():
            pass


    def integrate_gauss(f, lims, npts = 3):
        """Performs numerical integration of a function using Gauss-Legendre quadrature.

        Parameters
        ----------
        f:  reference to a callable object
            e.g., funciton of class that implements the __call__() method

        lims:   object, len == 2
                contains the lower and upper bounds of integration (x = a, x = b)

        npts:   integer, optional
                gives number of integration points to use
                Possible inputs are 1, 2, 3, 4, 5

        Returns
        -------
        float, float-like
            Provides the integral estimate

        Raises
        ------
        TypeError:
            If f is not callable

        ValueError:
            If lims does not have len == 2
            If lims[0] or lims[1] is not convertible to float
            If npts is not in [1, 2, 3, 4, 5]
        """