# Utility Functions for Numerical Integration

def integrate_newton(x, f, alg = "trap"):
    """Performs numerical integration of discrete data using Newton-Cotes rules.

    Parameters
    ----------
    x:  array-like, shape(m,n)
        contains the coordinates and values of the sample points

    f:  array_like, shape(m,n)
        contains the coordinates and values of the sample points

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