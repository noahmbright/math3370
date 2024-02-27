import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from math import floor

def dVdt(t, V, I, E):
    # Define your ordinary differential equation here
    dVdt = E - V + I(t)
    return dVdt

def threshold_crossed(t, V, I, E):
    # Event function to trigger when y crosses a threshold
    return V[0] - V_threshold

threshold_crossed.terminal = True

V_reset = -25
V_threshold = 0 
tf = 1000
E = -10


def integrate_with_reset(t0, tf, V0, I, E):
    V_start = [V0]
    t_start = t0
    ts = np.array([])
    Vs = np.array([])
    spike_ts = np.array([])

    while len(ts) == 0 or ts[-1] < tf:
        sol = solve_ivp(dVdt, (t_start, tf), V_start, events=(threshold_crossed), args=(I, E), method = 'LSODA')
        
        ts = np.concatenate((ts, sol.t))
        Vs = np.concatenate((Vs, sol.y[0]))
        spike_ts = np.concatenate((spike_ts, sol.t_events[0]))
        
        t_start = ts[-1]
        V_start = [V_reset]


    return ts, Vs, spike_ts

def I_step(t, t_end):
    return 1000/t_end*(t<t_end)

def square_wave(t, A, duty_cycle, T, t_end):
    if t > t_end:
        return 0

    quotient = floor(t/T)
    t_mod = t - quotient*T
    return A * (t_mod < (duty_cycle*T))

duty_cycle = .715
t_end = 75
square = lambda t : square_wave(t, 1/duty_cycle*1000/t_end, duty_cycle, 1, t_end)

omega = .00386
A = 1000/(500-np.sin(2000*omega)/4/omega)
sin2 = lambda t : A*np.sin(omega*t)**2

tau = .1
B = 1000/tau/(1 - np.exp(-1000/tau))
exp = lambda t: B*np.exp(-t/tau)

ts, Vs, spike_ts = integrate_with_reset(0, 1000, V_reset, exp, E)
plt.title(f"Spikes: {len(spike_ts)}")
plt.plot(ts, Vs)


plt.show()