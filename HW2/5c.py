import hh
import matplotlib.pyplot as plt
import numpy as np

#model = hh.HodgkinHuxley()
model = hh.Rinzel()
#model = hh.Kepler()

plt.figure()
Is = np.linspace(5.7, 5.75, 2)
for I in Is:
	X = model.integrate(I, [-65, model.n_inf(-65)], False)
	plt.plot(X.t, X.y[0], label = f"{I}")

plt.xlabel("t")
plt.ylabel("V")
plt.legend()
plt.show()