import numpy as np
import matplotlib.pyplot as plt

def analytical_n(alpha, beta, n0, t):
	C = np.sqrt(n0) - float(alpha)/(alpha+beta)
	return (float(alpha)/(alpha+beta) + C*np.exp(-(float(alpha)+beta)*t))**(2.0)

def gillespie(alpha, beta, t_max, n0, N):
	t = 0
	n = n0
	ns = [n]
	ts = [0]
	while t < t_max:
		r = (1-n)*alpha + beta*n
		rho = np.random.uniform()
		t = -1/r * np.log(rho)
		x = np.random.uniform(0, r) 
		# need a N?
		if x < alpha(N-n):
			n = min(N, n+1)
		else:
			n = max(0, n-1)

		ts.append(t)
		ns.append(n)

	return ts, ns


# part b
ts = np.linspace(0, 200, 500)
alpha = 0.04
beta = 0.02

analytical_ns = [analytical_n(alpha, beta, .1, t) for t in ts]
plt.plot(ts, analytical_ns)
plt.xlabel("t (ms)")
plt.ylabel("n")
plt.title("Analytical n")
plt.show()


# part c
Ts = []
Ns = []
for i in range(100):
	sim_ts, sim_ns = gillespie(alpha, beta, 200, .1, 2)
	Ts.append(sim_ts)
	Ns.append(sim_ns)

rng = np.random.default_rng()
rand_ints = rng.integers(101, size=10)

