import numpy as np

class HHDefaults:
    C = 1
    GKMax = 36
    GNaMax = 120
    Gm = 0.3
    EK = -12
    ENa = 115
    VRest = 10.613
    n = 0.32
    m = 0.05
    h = 0.6

    @staticmethod
    def random(s):
        return dict(C = HHDefaults.C * np.random.normal(1,s),
                    GKMax = HHDefaults.GKMax * np.random.normal(1,s),
                    GNaMax = HHDefaults.GNaMax * np.random.normal(1,s),
                    Gm = HHDefaults.Gm * np.random.normal(1,s),
                    EK = HHDefaults.EK * np.random.normal(1,s),
                    ENa = HHDefaults.ENa * np.random.normal(1,s),
                    VRest = HHDefaults.VRest * np.random.normal(1,s),
                    n = HHDefaults.n * np.random.normal(1,s),
                    m = HHDefaults.m * np.random.normal(1,s),
                    h = HHDefaults.h * np.random.normal(1,s))
                    
                    
    
class HodHux:
    def __init__(self, C = HHDefaults.C,
                 GKMax = HHDefaults.GKMax,
                 GNaMax = HHDefaults.GNaMax,
                 Gm = HHDefaults.Gm,
                 EK = HHDefaults.EK,
                 ENa = HHDefaults.ENa,
                 VRest = HHDefaults.VRest,
                 n = HHDefaults.n,
                 m = HHDefaults.m,
                 h = HHDefaults.h):

        self.C = C # in muF/cm^2
        self.GKMax = GKMax
        self.GNaMax = GNaMax
        self.EK = EK # mV
        self.ENa = ENa
        self.Gm = Gm
        self.VRest = VRest
        
        self.n = n
        self.m = m
        self.h = h

        self.reset(dt=None, V0=None)

    def reset(self, dt, V0):
        self.dt = dt
        self.V = V0
        
        self.INa = 0
        self.IK = 0
        self.Im = 0
         
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
        
    def simulate(self, I, dt, V0):
        self.reset(dt, V0)

        V = np.zeros_like(I)
        for i in range(len(I)):
            V[i] = self.step(I[i])

        return V
        
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

