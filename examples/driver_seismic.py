# Driver Script that uses discrete data to perform numerical integration
import numpy as np
import matplotlib.pyplot as plt
from goph420_lab01.integration import integrate_newton

def main():
    
    vel_data = np.loadtxt('data/s_wave_data.txt')
    # set up x (time) and y (velocity) values
    t = vel_data[:,0] # seconds (units)
    v = vel_data[:,1] # mm/s (units)

    # Determine length of event (period, T (s)), i.e., the integration interval
    v_max = np.abs(np.max(v))
    v_end = v_max * 0.005

    for i, _ in enumerate(v):
        if np.abs(v[i]) > v_end:
            int_limit = i

    print(f"Integration Index Limit = {int_limit}")
    print(f"0.5% of max event velocity (mm/s) = {v_end}\nEvent Period (T) = {t[int_limit]}\n")

    # plot raw data
    plt.plot(t, v, 'k-', label = 'S-Wave Arrivals', linewidth = 0.7)
    plt.vlines(x = t[int_limit], ymin = -0.1, ymax = 0.1, color = 'c', ls = '--', label = 'event endtime')
    plt.vlines(x = t[0], ymin = -0.1, ymax = 0.1, color = 'r', ls = '--', label = 'event starttime')
    plt.ylabel('velocity (mm/s)')
    plt.xlabel('time (s)')
    plt.title('Raw Data - S Wave Arrivals')
    plt.legend()
    plt.grid()
    plt.savefig('figures/raw_swave_data.png')
    plt.close('all')

    # integration using the trapezoid rule, and simpson's 1/3 and 3/8 rules
    T = t[int_limit] # determine length of event (period)
    t = t[:int_limit] # slice list to selected event period
    v = (v[:int_limit]) ** 2

    steps = [1, 2, 4, 8, 16, 32, 64, 128, 256]
    stepsize = np.array(steps) * 0.01
    
    # initialize integral and approx error variables for trapezoid rule
    int_trap = []
    int_simp = []

    for step in steps:
        Itrap = integrate_newton(t[::step], v[::step], alg= "trap")
        Isimp = integrate_newton(t[::step], v[::step], alg= "simp")
        int_trap.append(Itrap)
        int_simp.append(Isimp)

    eps_trap = np.abs(np.diff(int_trap)) / int_trap[:-1]
    eps_simp = np.abs(np.diff(int_simp) / int_simp[:-1])

    #plt.figure(figsize=(16, 16))
    plt.loglog(stepsize[:-1], eps_trap, label = 'trapezoid rule')
    plt.loglog(stepsize[:-1], eps_simp, label = 'simpson\'s rule')
    plt.legend()
    plt.ylabel('Approximate Relative Error (' + r'$\epsilon$' + 'a)')
    plt.title('Error Convergence')
    plt.xlabel('Sampling Interval (' + r'$\Delta$' + 't)')
    plt.savefig('figures/error_convergence.png')
    plt.close('all')

    print(f"Trapezoid Integral Estimates for stepsizes {stepsize}: \n{1 / T * np.array(int_trap)}\n")
    print(f"Simpson's Integral Estimate for stepsizes {stepsize}: \n{1/T * np.array(int_simp)}\n")


if __name__ == '__main__':
    main()