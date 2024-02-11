import hh
import matplotlib.pyplot as plt
import numpy as np

hh = hh.HodgkinHuxley()
Vs = np.linspace(-80, 150, 500)
Is = [hh.I_Na(V, hh.m_inf(V), hh.h_inf(V)) for V in Vs]
Is_hyperpolarized = [hh.I_Na(V, hh.m_inf(V), hh.h_inf(V, 20)) for V in Vs]

hs = [hh.h_inf(V) for V in Vs]
hs_hyperpolarized = [hh.h_inf(V, 20) for V in Vs]

plt.figure()
plt.xlabel("V (mV)")
plt.ylabel("h")
plt.plot(Vs, hs, label = 'normal h')
plt.plot(Vs, hs_hyperpolarized, label = 'hyperpolarized')
plt.legend()
plt.show()



plt.figure()
plt.plot(Vs, Is, label = 'normal h')
plt.plot(Vs, Is_hyperpolarized, label = 'hyperpolarized h')
plt.xlabel("V (mV)")
plt.ylabel("I")
plt.legend()
plt.show()