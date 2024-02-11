import hh
import numpy as np
import matplotlib.pyplot as plt

hh = hh.HodgkinHuxley()
X = hh.integrate(I_0 = 0)
V = X[:,0]
V_eq = V[-1]

ic = lambda V : [V, hh.m_inf(V_eq), hh.h_inf(V_eq), hh.n_inf(V_eq)]

V0s = np.linspace(-100, 100, 6)
plt.figure()
for V0 in V0s:
	X = hh.integrate(I_0 = 0, ics = ic(V0), clamp = True)
	I_Na = hh.I_Na(X[:,0], X[:,1], X[:,2])
	plt.plot(hh.t, I_Na, label = f'{V0}')

plt.xlabel("time")
plt.ylabel("I")
plt.legend()
plt.show()

n = 80
V0s = np.linspace(-80, 80, n)
peak_Is = np.zeros(n)
for i in range(n):
	X = hh.integrate(I_0 = 0, ics = ic(V0s[i]), clamp = True)
	I_Na = hh.I_Na(X[:,0], X[:,1], X[:,2])
	peak_Is[i] = np.max(np.abs(I_Na))

plt.figure()
plt.plot(V0s, peak_Is)
plt.xlabel("V0")
plt.ylabel("I")
plt.show()
