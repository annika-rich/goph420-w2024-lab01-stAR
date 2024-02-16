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

    def inegrate_gauss(f, lims, npts = 3):
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