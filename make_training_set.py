
import hodhux as hh
import stimulus as stim
import numpy as np
import h5py

def sim():
    Itau = 1
    dt = 0.025
    V0 = 0
    t_end = 500
    
    amp = np.random.uniform(-20, 50)
    t_i_start = int(np.random.uniform(0,100))
    t_i_stop = float("inf")
    while t_i_stop > t_end:
        t_i_stop = int(t_i_start + np.random.uniform(100,400))

    params = hh.HHDefaults.random(.1)
    h = hh.HodHux(**params)

    I = stim.step(t_i_start, t_i_stop, t_end, amp, dt)

    if Itau:
        I = stim.filter_I_tau(I, Itau, dt)

    V = h.simulate(I, dt, V0)

    return I, V, dict(amp=amp, t_i_start=t_i_start, t_i_stop=t_i_stop)

def main():
    N = 1000
    Is = []
    Vs = []
    
    for i in range(N):
        print("%d/%d" % (i+1, N))
        I,V,S = sim()
        Is.append(I)
        Vs.append(V)

    with h5py.File("hh_training.h5","w") as f:
        f.create_dataset("I", data=np.vstack(Is))
        f.create_dataset("V", data=np.vstack(Vs))

                


if __name__ == "__main__": main()
        
        

                      
