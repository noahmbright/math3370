import hh
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np

dp = hh.DestexhePare()

dp.set_tf(300)

Is = np.linspace(5, 50, 50)
dp.set_gkm(0)
Fs = []
for I in Is:
	X = dp.integrate(I, [-83.87, 0, .8, .002, .0075], False)
	peaks = find_peaks(X.y[0])
	T = X.t[peaks[0][-1]] - X.t[peaks[0][-2]]
	Fs.append(1/T)
plt.plot(Is, Fs, label = 'gkm = 0')

dp.set_gkm(5)
Fs = []
for I in Is:
	X = dp.integrate(I, [-83.87, 0, .8, .002, .0075], False)
	peaks = find_peaks(X.y[0])
	T = X.t[peaks[0][-1]] - X.t[peaks[0][-2]]
	Fs.append(1/T)
plt.plot(Is, Fs, label = 'gkm = 5')

plt.legend()

plt.xlabel("I")
plt.ylabel("F")
plt.show()