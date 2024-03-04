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
    print(f"0.5% of max event velocity (mm/s) = {v_end}\nTime at end of interval = {t[int_limit]}")

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
    v = v[:int_limit]

    # initialize integral and approx error variables for trapezoid rule
    Itrap = 0.0
    segment = 0.0
    error_trap = []
    int_trap = []
    delta_t = 0.0

    for i in range(len(t) - 1):
        segment = integrate_newton(t[i:i+2], v[i:i+2], alg= "trap")
        Itrap += segment
        eps_a = np.abs(segment / Itrap)
        error_trap.append(eps_a)
        delta_t += np.round(t[i+1] - t[i], decimals = 2)
        int_trap.append(delta_t)
    print(f"Trapezoid Integral Estimate: {Itrap}")
    # integral and approx error variables for simp rule
    Isimp = 0.0
    error_simp = []
    int_simp = []
    delta_t = 0.0

    for i in range(0,len(t) - 2, 4):
        segment = integrate_newton(t[i:i+5], v[i:i+5], alg= "simp")
        Isimp += segment
        eps_a = np.abs(segment / Isimp)
        error_simp.append(eps_a)
        delta_t += np.round(t[i+1] - t[i], decimals = 2)
        int_simp.append(delta_t)

    print(f"Simpson's Integral Estimate: {Isimp}")
    print(f"Check simp: {integrate_newton(t, v, alg = "simp")}")

    plt.figure(figsize=(16, 16))
    plt.loglog(int_trap, error_trap, label = 'trapezoid rule')
    plt.loglog(int_simp, error_simp, label = 'simpson\'s rule')
    plt.legend()
    plt.ylabel('Approximate Relative Error (' + r'$\epsilon$' + 'a)')
    plt.title('Error Convergence')
    plt.xlabel('Sampling Interval (' + r'$\Delta$' + 't)')
    plt.savefig('figures/error_convergence.png')


if __name__ == '__main__':
    main()