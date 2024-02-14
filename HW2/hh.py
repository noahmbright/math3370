import scipy as sp
import numpy as np
import pylab as plt
from scipy.integrate import solve_ivp

class HodgkinHuxley():

    #mS/cm^2
    g_Na = 120.0
    g_K  =  36.0
    g_L  =   0.3
    
    # mV
    E_Na =  50.0
    E_K  = -77.0
    E_L  = -54.387

    # TODO make phi do what phi is actually supposed to do
    phi = 1
    

    def __init__(self):
        self.t_on = 0
        self.I_p = 0
        self.pulse_width = 0    
        self.cm = 1.0 # uF/cm^2
        self.t0 = 0
        self.tf = 400

    def set_t0(self, t0):
        self.t0 = t0

    def set_tf(self, tf):
        self.tf = tf

    def set_cm(self, cm):
        self.cm = cm
    
    def alpha_m(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 0.1*(V+40.0)/(1.0 - np.exp(-(V+40.0) / 10.0))

    def beta_m(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 4.0*np.exp(-(V+65.0) / 18.0)

    def alpha_h(self, V, hyperpolarization = 0):
        """Channel gating kinetics. Functions of membrane voltage"""
        V_half = -65.0 - hyperpolarization
        return 0.07*np.exp(-(V-V_half) / 20.0)

    def beta_h(self, V, hyperpolarization = 0):
        """Channel gating kinetics. Functions of membrane voltage"""
        V_half = -35.0 - hyperpolarization
        return 1.0/(1.0 + np.exp(-(V-V_half) / 10.0))

    def alpha_n(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 0.01*(V+55.0)/(1.0 - np.exp(-(V+55.0) / 10.0))

    def beta_n(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 0.125*np.exp(-(V+65) / 80.0)

    def tau_h(self, V, hyperpolarization):
        return self.alpha_h(V, hyperpolarization) + self.beta_h(V, hyperpolarization)

    def tau_m(self, V):
        return self.alpha_m(V) + self.beta_m(V)

    def tau_n(self, V):
        #return self.alpha_n(V) + self.beta_n(V)
        return .1

    def m_inf(self, V):
        return self.alpha_m(V)/self.tau_m(V)

    def n_inf(self, V):
        return self.alpha_n(V)/self.tau_n(V)

    def h_inf(self, V, hyperpolarization=0):
        return self.alpha_h(V, hyperpolarization)/self.tau_h(V, hyperpolarization)

    def I_Na(self, V, m, h):
        """
        Membrane current (in uA/cm^2)
        Sodium (Na = element name)

        |  :param V:
        |  :param m:
        |  :param h:
        |  :return:
        """
        return self.g_Na * m**3 * h * (V - self.E_Na)

    def I_K(self, V, n):
        """
        Membrane current (in uA/cm^2)
        Potassium (K = element name)

        |  :param V:
        |  :param h:
        |  :return:
        """
        return self.g_K  * n**4 * (V - self.E_K)
    #  Leak
    def I_L(self, V):
        """
        Membrane current (in uA/cm^2)
        Leak

        |  :param V:
        |  :param h:
        |  :return:
        """
        return self.g_L * (V - self.E_L)

    #injected current
    def set_current_pulse(self, I_p, t_on, pulse_width):
        self.t_on = t_on
        self.I_p = I_p 
        self.pulse_width = pulse_width

    def I_inj(self, t):
        """
        External Current

        |  :param t: time
        |  :return: step up to 10 uA/cm^2 at t>100
        |           step down to 0 uA/cm^2 at t>200
        |   this is a pulse of current and you should modify it
        """
        return self.I_p*(t>self.t_on) - self.I_p*(t>(self.t_on + self.pulse_width)) 

    @staticmethod
    def dALLdt(t, X, self, I_0 = 0, clamp = False):
        """
        Integrate

        |  :param X:
        |  :param t:
        |  :return: calculate membrane potential & activation variables
        """
        V, m, h, n = X

        if clamp:
            dVdt = 0
        else:
            dVdt = (self.I_inj(t) + I_0 - self.I_Na(V, m, h) - self.I_K(V, n) - self.I_L(V)) / self.cm
        dmdt = self.alpha_m(V)*(1.0-m) - self.beta_m(V)*m
        dhdt = self.alpha_h(V)*(1.0-h) - self.beta_h(V)*h
        dndt = self.alpha_n(V)*(1.0-n) - self.beta_n(V)*n
        return np.array([dVdt, dmdt, dhdt, dndt])

    def integrate(self, I_0 = 0, ics = [-60, 0.09, 0.4, 0.4], clamp = False):
        #return odeint(self.dALLdt, ics, self.ts, args=(self, I_0, clamp)) # V, m , h, 
        return solve_ivp(self.dALLdt, (self.t0, self.tf), ics, args = (self, I_0, clamp))
        

    def Main(self):
        """
        Main demo for the Hodgkin Huxley neuron model
        """

        X = self.integrate(I_0 = 20)
        ts = X.t
        V = X.y[0]
        m = X.y[1]
        h = X.y[2]
        n = X.y[3]
        
        ina = self.I_Na(V, m, h)
        ik = self.I_K(V, n)
        il = self.I_L(V)

        plt.figure()

        plt.subplot(4,1,1)
        plt.title('Hodgkin-Huxley Neuron')
        plt.plot(ts, V, 'k')
        plt.ylabel('V (mV)')

        plt.subplot(4,1,2)
        plt.plot(ts, ina, 'c', label='$I_{Na}$')
        plt.plot(ts, ik, 'y', label='$I_{K}$')
        plt.plot(ts, il, 'm', label='$I_{L}$')
        plt.ylabel('Current')
        plt.legend()

        plt.subplot(4,1,3)
        plt.plot(ts, m, 'r', label='m')
        plt.plot(ts, h, 'g', label='h')
        plt.plot(ts, n, 'b', label='n')
        plt.ylabel('Gating Value')
        plt.legend()

        plt.subplot(4,1,4)
        i_inj_values = [self.I_inj(t) for t in ts]
        plt.plot(ts, i_inj_values, 'k')
        plt.xlabel('t (ms)')
        plt.ylabel('$I_{inj}$ ($\\mu{A}/cm^2$)')
        plt.ylim(-1, 40)

        plt.show()




class Rinzel(HodgkinHuxley):
    """
    inherit all the alpha/beta from HH
    use the approximation that h = h0 - n
    """
    def __init__(self):
        super().__init__()
        self.h0 = .8

    def set_h0(self, h0):
        self.h0 = h0

    def h(self, n):
        return self.h0 - n

    @staticmethod
    def dALLdt(t, X, self, I_0 = 0, clamp = False):
        V, n = X
        if clamp:
            dVdt = 0
        else:
            dVdt = (self.I_inj(t) + I_0 - self.I_Na(V, self.m_inf(V), self.h(n)) - self.I_K(V, n) - self.I_L(V)) / self.cm
            dndt = self.alpha_n(V)*(1.0-n) - self.beta_n(V)*n

        return np.array([dVdt, dndt])


class Kepler(HodgkinHuxley):
    def __init__(self):
        super().__init__()

    def dh_inf(self, V, dV = .01):
        return (self.h_inf(V + dV) - self.h_inf(V - dV))/(2*dV)

    @staticmethod
    def dALLdt(t, X, self, I_0 = 0, clamp = False):
        V, Vh = X
        if clamp:
            dVdt = 0
        else:
            dVdt = (self.I_inj(t) + I_0 - self.I_Na(V, self.m_inf(V), self.h_inf(Vh)) - self.I_K(V, self.n_inf(Vh)) - self.I_L(V)) / self.cm
            dVhdt = (self.h_inf(V) - self.h_inf(Vh))/(self.dh_inf(Vh)/(self.alpha_h(V) + self.beta_h(V)))

        return np.array([dVdt, dVhdt])



class DestexhePare(HodgkinHuxley):
    def __init__(self):
        super().__init__()
        self.E_Na = 55
        self.E_k = -85
        self.g_K = 100
        self.gkm = 2

    def set_gkm(self, gkm):
        self.gkm = gkm

    V_t = -58
    V_s = -10

    def alpha_m(self, V):
        return -.32*(V-self.V_t-13)/(np.exp(-(V-self.V_t-13)/4)-1)

    def beta_m(self, V):
        return .28*(V-self.V_t-40)/(np.exp((V-self.V_t-40)/5)-1)

    def alpha_h(self, V):
        return .128*np.exp(-(V-self.V_t-self.V_s-17)/18)

    def beta_h(self, V):
        return 4/(1+np.exp(-(V-self.V_t-self.V_s-40)/5))

    def alpha_n(self, V):
        return -.032*(V-self.V_t-15)/(np.exp(-(V-self.V_t-15)/5)-1)

    def beta_n(self, V):
        return .5*np.exp(-(V-self.V_t-10)/40)

    def alpha_km(self, V):
        return .0001*(V+30)/(1-np.exp(-(V+30)/9))

    def beta_km(self, V):
        return -.0001*(V+30)/(1-np.exp((V+30)/9))

    def I_km(self, V, m):
        return self.gkm * m * (V-self.E_k)

    @staticmethod
    def dALLdt(t, X, self, I_0 = 0, clamp = False):
        """
        Integrate

        |  :param X:
        |  :param t:
        |  :return: calculate membrane potential & activation variables
        """
        V, m, h, n, mk = X

        if clamp:
            dVdt = 0
        else:
            dVdt = (self.I_inj(t) + I_0 - self.I_Na(V, m, h) - self.I_K(V, n) - self.I_L(V) - self.I_km(V, mk)) / self.cm
        dmdt = self.alpha_m(V)*(1.0-m) - self.beta_m(V)*m
        dhdt = self.alpha_h(V)*(1.0-h) - self.beta_h(V)*h
        dndt = self.alpha_n(V)*(1.0-n) - self.beta_n(V)*n
        dmkdt = self.alpha_km(V)*(1.0-mk) - self.beta_km(V)*mk
        return np.array([dVdt, dmdt, dhdt, dndt, dmkdt])



if __name__ == '__main__':
    runner = HodgkinHuxley()
    runner.set_current_pulse(3, 100, 10)
    runner.Main()

    
