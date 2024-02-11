import numpy as np
import matplotlib.pyplot as plt

g_cl_dark = 57500000*1E-8
g_na_dark = .005
g_k_dark = 1

g_cl_light = 57500000*1E-8
g_na_light = 20
g_k_light = 1

E_cl = -50
E_na = 50
E_k = -90

V_dark = ( g_cl_dark * E_cl + g_na_dark * E_na + g_k_dark * E_k) / ( g_cl_dark + g_k_dark + g_na_dark)
V_light = ( g_cl_light * E_cl + g_na_light * E_na + g_k_light * E_k) / ( g_cl_light + g_k_light + g_na_light)

a = (g_k_dark*E_k + g_na_dark*E_na)*1E8
b = (g_k_dark + g_na_dark)*1E8
v2 = -75
new_g = ( a - b*v2 ) / ( v2 - E_cl )

print(new_g)
print(V_dark)
print(V_light)

a1 = 10E-4 #cm
h1 = 100E-4 
a2 = 25E-4
h2 = 50E-4
a3 = 20E-4
h3 = 70E-4
rL = 100

g12 = ( a1 * a2**2 ) / ( rL * h1 * ( a2**2 * h1 + a1**2 * h2 ) )
g21 = ( a1**2 * a2 ) / ( rL * h2 * ( a2**2 * h1 + a1**2 * h2 ) )
g23 = ( a2 * a3**2 ) / ( rL * h2 * ( a3**2 * h2 + a2**2 * h3 ) )
g32 = ( a2**2 * a3 ) / ( rL * h3 * ( a3**2 * h2 + a2**2 * h3 ) )
g = 1/(2E4)
print("Gs")
print(g12)
print(g21)
print(g23)
print(g32)

A = g + g21 + g23 - (g23 * g32) / (g32 + g)
B = g + g12 - (g12 * g21)/A
I1 = 1E-6
#I3 = 1E-6
I3 = 0

V1 = ( g12/A * g23*I3 / (g32 + g)  + I1) / B
V2 = g21/A * V1 + 1/A * g23*I3/(g32 + g)
V3 = (g32*V2 + I3) / (g32 + g)

print("Vs")
print(V1)
print(V2)
print(V3)

test1 = -g*V1 + g12 * (V2 - V1) + I1
test2 = -g*V2 + g21*(V1-V2) + g23*(V3 - V2)
test3 = -g*V3 + g32*(V2 - V3) + I3

print("tests")
print(test1)
print(test2)
print(test3)


# def v(x,t):
# 	return np.exp(-(x**2/t + t)) / t**(1.0/2.0)

# xs = [1, 1.25, 1.5, 2, 2.5, 3]
# ts = np.linspace(.1, 5, 100)
# for x in xs:
# 	vs = [v(x, t) for t in ts]
# 	plt.plot(ts, vs, label = f'{x}')

# plt.legend()
# plt.show()
