import hh
import numpy as np
import matplotlib.pyplot as plt

hh = hh.HodgkinHuxley()

plt.figure()

X = hh.integrate(0, [-39, hh.m_inf(-75), hh.h_inf(-75), hh.n_inf(-75)])
g_K = hh.g_K * np.power(X[:, 3], 4.0)
g_Na = hh.g_Na * np.power(X[:, 1], 3.0) * X[:, 2]
plt.plot(hh.t, g_K, label = 'g_K')
plt.plot(hh.t, g_Na, label = "g_Na")

X_0 = hh.integrate(0, [0, hh.m_inf(-75), hh.h_inf(-75), hh.n_inf(-75)])
g_K_0 = hh.g_K * np.power(X_0[:, 3], 4.0)
g_Na_0 = hh.g_Na * np.power(X_0[:, 1], 3.0) * X_0[:, 2]
plt.plot(hh.t, g_K_0, label = 'g_K 0')
plt.plot(hh.t, g_Na_0, label = "g_Na 0")

plt.legend()
plt.xlabel("t")
plt.ylabel('g')
plt.show()


plt.figure()
X = hh.integrate(10, [-65, hh.m_inf(-75), hh.h_inf(-75), hh.n_inf(-75)])
m_infs = hh.h_inf(X[:, 0])
plt.plot(hh.t, m_infs, label = "h_{inf}(V(t))")
plt.plot(hh.t, X[:, 2], label = "h_(t)")
plt.xlabel("t")
plt.ylabel("h")
plt.legend()
plt.show()