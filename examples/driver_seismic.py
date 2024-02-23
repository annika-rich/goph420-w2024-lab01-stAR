# Driver Script that uses discrete data to perform numerical integration
import numpy as np
import matplotlib.pyplot as plt

def main():
    
    vel_data = np.loadtxt('data/s_wave_data.txt')
    # set up x (time) and y (velocity) values
    t = vel_data[:,0] # seconds (units)
    v = vel_data[:,1] # mm/s (units)

    # Determine length of event (period, T (s))
    v_max = np.abs(np.max(v))
    v_end = v_max * 0.005
    print(f"0.5% of max event velocity (mm/s) = {v_end}\nTime at end of interval = {t[436]}")

    # plot raw data
    plt.plot(t, v, 'k-', label = 'S-Wave Arrivals', linewidth = 0.7)
    plt.vlines(x = t[436], ymin = -0.1, ymax = 0.1, color = 'c', ls = '--', label = 'event endtime')
    plt.vlines(x = t[0], ymin = -0.1, ymax = 0.1, color = 'r', ls = '--', label = 'event starttime')
    plt.ylabel('velocity (mm/s)')
    plt.xlabel('time (s)')
    plt.title('Raw Data - S Wave Arrivals')
    plt.legend()
    plt.grid()
    plt.savefig('figures/raw_swave_data.png')

if __name__ == '__main__':
    main()