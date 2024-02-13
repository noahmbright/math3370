import numpy as np
import matplotlib.pyplot as plt

def analytical_n(alpha, beta, n0, t):
	C = np.sqrt(n0) - float(alpha)/(alpha+beta)
	return (float(alpha)/(alpha+beta) + C*np.exp(-(float(alpha)+beta)*t))**(2.0)

def gillespie(alpha, beta, t_max, n0, N):
	t = 0
	if np.random.uniform() < n0:
		n = N
	else:
		n = np.random.randint(0, N-1)

	ns = [int(n == N)]
	ts = [0]
	while t < t_max:
		r = (N-n)*alpha + beta*n
		rho = np.random.uniform()
		t += -1/r * np.log(rho)
		x = np.random.uniform(0, r)

		if x < alpha*(N-n):
			n = min(N, n+1)
		else:
			n = max(0, n-1)

		ts.append(t)
		ns.append(int(n == N))

	return np.array(ts), np.array(ns)


# part b
ts = np.linspace(0, 200, 30)
alpha = 0.04
beta = 0.02

analytical_ns = [analytical_n(alpha, beta, .1, t) for t in ts]
plt.figure()
plt.plot(ts, analytical_ns, label = 'analytical')
plt.xlabel("t (ms)")
plt.ylabel("n")
#plt.show()


# part c
n = 100
Ts = []
Ns = []
for i in range(n):
	sim_ts, sim_ns = gillespie(alpha, beta, 200, .1, 2)
	Ts.append(sim_ts)
	Ns.append(sim_ns)


rng = np.random.default_rng()
rand_ints = rng.integers(n, size=10)

ns10 = np.zeros(len(ts))
ns100 = np.zeros(len(ts))
for rand in rand_ints:
	ns10 += np.interp(ts, Ts[rand], Ns[rand])

for i in range(n):
	ns100 += np.interp(ts, Ts[i], Ns[i])

plt.plot(ts, ns10/10, label = '10 samples')
plt.plot(ts, ns100/n, label = f'{n} samples')
plt.legend()
plt.show()

