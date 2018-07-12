import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt
from hodhux import HodHux
import stimulus as stim

import numpy as np

def sim():
    DC = 35
    t_i_start = 25
    t_i_stop = 200
    t_end = 250
    Itau = 1
    dt = 0.025
    V0 = 0
 
    hh = HodHux()

    I = stim.step(t_i_start, t_i_stop, t_end, DC, dt)

    if Itau:
        I = stim.filter_I_tau(I, Itau, dt)

    V = hh.simulate(I, dt, V0)
    t = np.arange(0.0, t_end, dt)

    fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(5,6))
    ax1.plot(t, I)
    ax1.set_ylabel('stimulus (nA/cm^2)')
    
    ax2.plot(t, V)
    ax2.set_xlabel('time (ms)')
    ax2.set_ylabel('response (mV)')

    plt.suptitle('Hodgkin-Huxley Step Response')
    plt.show()
    plt.savefig('example.png')


if __name__ == "__main__": sim()
