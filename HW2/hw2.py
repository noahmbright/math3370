import numpy as np
import matplotlib.pyplot as plt

def analytical_n(alpha, beta, n0, t):
	gamma = alpha +  beta
	C = alpha/beta - np.sqrt(n0)
	D = 3*alpha**2.0/gamma**2.0 - 4*alpha/gamma*np.sqrt(n0)

	return D + 2*C* ( alpha/gamma*np.exp(-gamma*t) + C/2.0 * np.exp(-2*gamma*t) )

def gillespie(alpha, beta, t_max, n0):
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
			n += 1
		else:
			n -= 1

		ts.append(t)
		ns.append(n)

	return ts, ns

