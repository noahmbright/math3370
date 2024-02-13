import hh
import matplotlib.pyplot as plt
import numpy as np

hh = hh.HodgkinHuxley()

plt.figure()
Is = np.linspace(6.25, 6.3, 2)
for I in Is:
	X = hh.integrate(I, [-65, hh.m_inf(-65), hh.h_inf(-65), hh.n_inf(-65)])
	plt.plot(hh.t, X[:, 0], label = f"{I}")

plt.xlabel("t")
plt.ylabel("V")
plt.legend()
plt.show()