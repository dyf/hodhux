import numpy as np

def step(t_start, t_stop, t_end, DC, dt):
    t = np.arange(0.0, t_end, dt)

    I = np.zeros(t.shape)
    I[(t>=t_start)&(t<t_stop)] = DC
    
    return I

def ramp(t_start, t_stop, t_end, imax, dt):
    t = np.arange(0.0, t_end, dt)

    I = np.zeros(t.shape)
    I[(t>=t_start)&(t<t_stop)] = np.linspace(0, imax, (t_stop-t_start)/dt)
    return I
    
    
def filter_I_tau(I, Itau, dt):
    Iout = np.zeros_like(I)
        
    Icurr = I[0]
    for idx, i in enumerate(I):
        Icurr += dt/Itau * (i - Icurr)
        Iout[idx] = Icurr

    return Iout
