# Utility Functions for Numerical Integration

import numpy as np

def integrate_newton(x, f, alg = "trap"):
    """Performs numerical integration of discrete data using Newton-Cotes rules.

    Parameters
    ----------
    x:  array-like
        contains the coordinates and values of the sample points
        Sample points must be evenly spaced.

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
        If the string flag contains a string other than trap or simp
        If the dimensions of x and f are incompatible
    """
    # make sure inputs are arrays
    x = np.array(x, dtype = float)
    f = np.array(f, dtype = float)

    # check that x and f have the same length
    if (m := len(x)) != (n := len(f)):
        raise ValueError(f"The length of x ({m}) is not equal to the length of f ({n}, x and f must be the same shape")

    # check dimensions of inputs
    if (mdim := len(x.shape)) != (ndim := len(f.shape)) or (mdim not in [1]):
        raise ValueError(f"x has {mdim} dimensions and f has {ndim} dimensions, dimensions must be equal and 1D")

    # check string flag
    alg = alg.strip().lower()
    if alg not in ["trap", "simp"]:
        raise ValueError(f"The algorithm string flag, contains a string ({alg}) other than 'trap' or 'simp'.")
    
    # set spacing between data points
    delta_x = x[1] - x[0]

    if alg == "trap":
        integral = (delta_x / 2) * (f[0] + 2 * np.sum(f[1:-1]) + f[-1])

    elif alg == "simp":
        # Simpson's 1/3 rule (odd number of segments)
        if m % 2:
            integral = (delta_x / 3) * (f[0] + 2 * np.sum(f[2:-1:2]) + 4 * np.sum(f[1:-1:2]) + f[-1])
       # Simpson's 1/3 and 3/8th rule (even number of segments)
        else:
            integral = 0.0
            if m > 4:
                integral = (delta_x / 3) * (f[0] + 2 * np.sum(f[2:-4:2]) + 4 * np.sum(f[1:-4:2]) + f[-4])
            integral += (0.375 * delta_x) * (f[-4] + 3 * f[-3] + 3 * f[-2] + f[-1])
    return integral

def integrate_gauss(f, lims, npts = 3):
    """Performs numerical integration of a function using Gauss-Legendre quadrature.

    Parameters
    ----------
    f:  reference to a callable object
        e.g., function of class that implements the __call__() method

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
    if callable(f) == False:
        raise TypeError(f"{f} is not callable, f must be a callable object such as a function of a class that implements the __call__() method.")
    
    if (m := len(lims)) not in [2]:
        raise ValueError(f"Lims has length {m}, lims must be of length 2.")
    # set upper and lower bounds of integration
    a = float(lims[0])
    b = float(lims[1])

    if npts not in [1, 2, 3, 4, 5]:
        raise ValueError(f"The number of integration points to use (npts) must be in [1, 2, 3, 4, 5], {npts} is not a valid argument.")

def normal_density(x, mu, sigma):
    f = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mu)**2 / (sigma ** 2)))
