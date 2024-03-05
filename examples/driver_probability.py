# Driver Script that uses a known function to perform numerical integration
import numpy as np
import matplotlib.pyplot as plt
from goph420_lab01.integration import integrate_gauss

def main():
    # intialize annual mean, standard deviation
    sigma = 0.5
    mu = 1.5
    f = lambda x: 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mu) ** 2 / (sigma ** 2)))
    npts = [1, 2, 3, 4, 5]
    conv_ints = []
    error = []
    I_old = 0.0
    

    # compute exceedance probability of a M > 4.0 relative event
    for i in npts:
        I_z = integrate_gauss(f, lims = [1.5, 4], npts=i)
        eps_a = np.abs((I_z - I_old) / I_z)
        conv_ints.append(I_z)
        error.append(eps_a)
        I_old = I_z
    
    print(f"Exceedance Probability of a seismic event with M> 4.0 for npts {npts}: \n{0.5 - np.array(conv_ints)}\n")

    plt.loglog(npts, error, 'r')
    plt.title('P(M >= 4) for a Seismic Event')
    plt.ylabel('Approximate Relative Error')
    plt.xlabel('Number of Points (npts)')
    plt.savefig('figures/exceedance_probability_convergence.png')
    plt.close('all')

    # distance calculations
    sigma = 0.05
    mu = 10.28
    f = lambda x: 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mu) ** 2 / (sigma ** 2)))
    conv_ints = []
    error = []
    I_old = 0.0

    # compute confidence interval
    for i in npts:
        I_z = integrate_gauss(f, lims = [10.25, 10.35], npts=i)
        eps_a = np.abs((I_z - I_old) / I_z)
        conv_ints.append(I_z)
        error.append(eps_a)
        I_old = I_z

    print(f"Confidence Interval Probability that length (L) falls between 10. 25 - 10.23 m using npts {npts}: \n{2 * np.array(conv_ints)}")

    plt.loglog(npts, error, 'k')
    plt.title('P(10.25 <= L <= 10.35) for Length (L)')
    plt.ylabel('Approximate Relative Error')
    plt.xlabel('Number of Points (npts)')
    plt.savefig('figures/confidence_interval_convergence.png')
    plt.close('all')

if __name__ == '__main__':
    main()