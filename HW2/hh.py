import scipy as sp
import numpy as np
import pylab as plt
from scipy.integrate import odeint

class HodgkinHuxley():
    C_m  =   1.0 # uF/cm^2

    #mS/cm^2
    g_Na = 120.0
    g_K  =  36.0
    g_L  =   0.3

    # mV
    E_Na =  50.0
    E_K  = -77.0
    E_L  = -54.387
    
    t = np.arange(0.0, 450.0, 0.01)
    

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

    def m_inf(self, V):
        return self.alpha_m(V)/(self.alpha_m(V) + self.beta_m(V))

    def n_inf(self, V):
        return self.alpha_n(V)/(self.alpha_n(V) + self.beta_n(V))

    def h_inf(self, V, hyperpolarization=0):
        return self.alpha_h(V, hyperpolarization)/(self.alpha_h(V, hyperpolarization) + self.beta_h(V, hyperpolarization))

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

    def I_inj(self, t, t_on = 100, I_p = 0.0, pulse_width = 0):
        """
        External Current

        |  :param t: time
        |  :return: step up to 10 uA/cm^2 at t>100
        |           step down to 0 uA/cm^2 at t>200
        |   this is a pulse of current and you should modify it
        """
        return I_p*(t>t_on) - I_p*(t>t_on + pulse_width) 

    @staticmethod
    def dALLdt(X, t, self, I_0, clamp = False):
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
            dVdt = (self.I_inj(t) + I_0 - self.I_Na(V, m, h) - self.I_K(V, n) - self.I_L(V)) / self.C_m
        dmdt = self.alpha_m(V)*(1.0-m) - self.beta_m(V)*m
        dhdt = self.alpha_h(V)*(1.0-h) - self.beta_h(V)*h
        dndt = self.alpha_n(V)*(1.0-n) - self.beta_n(V)*n
        return dVdt, dmdt, dhdt, dndt

    def integrate(self, I_0 = 8, ics = [-60, 0.09, 0.4, 0.4], clamp = False):
        return odeint(self.dALLdt, ics, self.t, args=(self, I_0, clamp)) # V, m , h, n
        

    def Main(self):
        """
        Main demo for the Hodgkin Huxley neuron model
        """

        X = self.integrate(I_0 = 20)
        V = X[:,0]
        m = X[:,1]
        h = X[:,2]
        n = X[:,3]
        
        ina = self.I_Na(V, m, h)
        ik = self.I_K(V, n)
        il = self.I_L(V)

        plt.figure()

        plt.subplot(4,1,1)
        plt.title('Hodgkin-Huxley Neuron')
        plt.plot(self.t, V, 'k')
        plt.ylabel('V (mV)')

        plt.subplot(4,1,2)
        plt.plot(self.t, ina, 'c', label='$I_{Na}$')
        plt.plot(self.t, ik, 'y', label='$I_{K}$')
        plt.plot(self.t, il, 'm', label='$I_{L}$')
        plt.ylabel('Current')
        plt.legend()

        plt.subplot(4,1,3)
        plt.plot(self.t, m, 'r', label='m')
        plt.plot(self.t, h, 'g', label='h')
        plt.plot(self.t, n, 'b', label='n')
        plt.ylabel('Gating Value')
        plt.legend()

        plt.subplot(4,1,4)
        i_inj_values = [self.I_inj(t) for t in self.t]
        plt.plot(self.t, i_inj_values, 'k')
        plt.xlabel('t (ms)')
        plt.ylabel('$I_{inj}$ ($\\mu{A}/cm^2$)')
        plt.ylim(-1, 40)

        plt.show()

if __name__ == '__main__':
    runner = HodgkinHuxley()
    runner.Main()

    
