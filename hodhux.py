import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt

import numpy as np

class HodHux:
    def __init__(self, C, GKMax, GNaMax,
                 EK, ENa, Gm, VRest, V0):

        self.C = C #in muF/cm^2
        self.GKMax = GKMax
        self.GNaMax = GNaMax
        self.EK = EK # mV
        self.ENa = ENa
        self.Gm = Gm
        self.VRest = VRest
        self.V = V0

        self.dt = 0.025 # ms
        self.n = 0.32
        self.m = 0.05
        self.h = 0.60
        self.INa = 0
        self.IK = 0
        self.Im = 0

    def filter_I(self, I, Itau):
        Iout = np.zeros_like(I)
        
        Icurr = I[0]
        for idx, i in enumerate(I):
            Icurr += self.dt/Itau * (i - Icurr)
            Iout[idx] = Icurr

        return Iout
         
    def alphaN(self, V):
        if V==10: return self.alphaN(V+0.001) # 0/0 -> NaN
        return (10-V) / (100*(np.exp((10-V)/10)-1))

    def betaN(self, V):
        return 0.125 * np.exp(-V/80)

    def alphaM(self, V):
        if V==25: return self.alphaM(V+0.001) # 0/0 -> NaN
        return (25-V) / (10 * (np.exp((25-V)/10)-1))

    def betaM(self, V):
        return 4 * np.exp(-V/18)

    def alphaH(self, V):
        return 0.07 * np.exp(-V/20)
                    
    def betaH(self, V):
        return 1 / (np.exp((30-V)/10)+1)

    def step(self, Iinj):
        aN = self.alphaN(self.V)
        bN = self.betaN(self.V)
        aM = self.alphaM(self.V)
        bM = self.betaM(self.V)
        aH = self.alphaH(self.V)
        bH = self.betaH(self.V)
        
        tauN = 1 / (aN + bN)
        tauM = 1 / (aM + bM)
        tauH = 1 / (aH + bH)
        nInf = aN * tauN
        mInf = aM * tauM
        hInf = aH * tauH
        
        self.n += self.dt / tauN * (nInf - self.n)
        self.m += self.dt / tauM * (mInf - self.m)
        self.h += self.dt / tauH * (hInf - self.h)
        
        self.INa = self.GNaMax * self.m * self.m * self.m * self.h * (self.ENa - self.V)
        self.IK = self.GKMax * self.n * self.n * self.n * self.n * (self.EK - self.V)
        self.Im = self.Gm * (self.VRest - self.V)
        
        self.V += self.dt / self.C * (self.INa + self.IK + self.Im + Iinj)

        return self.V


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
             

