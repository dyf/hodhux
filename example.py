import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt
from hodhux import HodHux
import numpy as np

def sim():
    DC = 35
    t_i_start = 25
    t_i_stop = 200
    t_end = 250
    Itau = 1
 
    hh = HodHux(C=1,
                GKMax=36, GNaMax=120, Gm=0.3,
                EK=-12, ENa=115,
                VRest=10.613, V0=0)

    t = np.arange(0.0, t_end, hh.dt)
    V = np.zeros(t.shape)
    Iraw = np.zeros(t.shape)
    Iraw[(t>=t_i_start)&(t<t_i_stop)] = DC
    Iinj = hh.filter_I(Iraw, Itau)

    for i,tt in enumerate(t):
        iraw = Iraw[i]

        V[i] = hh.step(Iinj[i])

    fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(5,6))
    ax1.plot(t, Iinj)
    ax1.set_ylabel('stimulus (nA/cm^2)')
    
    ax2.plot(t, V)
    ax2.set_xlabel('time (ms)')
    ax2.set_ylabel('response (mV)')

    plt.suptitle('Hodgkin-Huxley Step Response')
    plt.show()
    plt.savefig('example.png')


if __name__ == "__main__": sim()
