#Import packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

#Physical constants
sigma = 5.67e-8
Q0 = 1368/4
eps0 = 0.62

#Numerical constants
Tmin = 220
Tmax = 310
T = np.linspace(Tmin, Tmax, 100)
steps = 1000 #Reduce this to decrease runtime

#Usefull functions
def LHS(T, Q):
    alpha = 1/2 - np.tanh((T-265)/10)/5
    return Q*(1-alpha)

def RHS(T):
    alpha = 1/2 - np.tanh((T-265)/10)/5
    return eps0*sigma*np.power(T,4)

def f(T,Q,eps):
    alpha = 1/2 - np.tanh((T-265)/10)/5
    return Q*(1-alpha) - eps*sigma*np.power(T,4)

#Central difference method to compute the derivative
def central_diff(f, a, Q, eps, h=0.01):
    return (f(a + h, Q, eps) - f(a - h, Q, eps))/(2*h)

#Function for computing the bifurcation diagram
def compute_bifurcation(T_list, P_list):
    #Temperature step size
    dT = np.abs(T_list[0] - T_list[1])

    #Lists for saving
    P_s = []
    P_us = []
    T_s = []
    T_us = []

    #Loop over all initial conditions
    for T0 in T_list:
        for P in P_list:
            #Try to find the equilibrium
            try: 
                #Find the root of the function
                if P_list[0]==0: T = T = brentq(f, T0-1, T0+1, args=(Q0, P))
                else: T = brentq(f, T0-dT, T0+dT, args=(P,eps0))

                #Determine the derivative at the root using central differences
                if P_list[0]==0: deriv = central_diff(f, T, Q0, P, 0.01)
                else: deriv = central_diff(f, T, P, eps0, 0.01)

                #Check stability, if stable:
                if deriv > 0:
                    T_us.append(T)
                    P_us.append(P)
                
                #If unstable
                else:
                    T_s.append(T)
                    P_s.append(P)
            
            #If equilibrium can not be found, continue
            except: continue
    return np.array(T_s), np.array(T_us), np.array(P_s), np.array(P_us)


#Plotting RHS and LHS for the standard values
plt.figure()
plt.title("LHS and RHS of the energy balance equilibrium equation")
plt.xlabel("T [K]")
plt.plot(T, LHS(T, Q0), color="tab:blue", label="LHS")
plt.plot(T, RHS(T), color="tab:red", label="RHS")
plt.legend()
plt.grid()
plt.savefig('BaseScenario.png')

#Values T and parameters for which we try to find the equilibrium
Q_list = np.linspace(0.5, 1.5, steps)*Q0
T_list = np.linspace(Tmin, Tmax, steps)
eps_list = np.linspace(0, 1, steps)

#Computing the Q bifurcation diagram
T_s, T_us, Q_s, Q_us = compute_bifurcation(T_list, Q_list)

#Compute some stuff for colors in the plot
m = np.argmin(np.abs(T_s - T_us[0]))
n = np.argmin(np.abs(T_s - T_us[-1]))
Q_us = np.concatenate(([Q_s[m]], Q_us, [Q_s[n]]))
T_us = np.concatenate(([T_s[m]], T_us, [T_s[n]]))

#Plot Q bifurcation diagram
plt.figure()
plt.title("Q bifurcation diagram energy balance model")
plt.plot(Q_s[:m]/Q0, T_s[:m], color="tab:purple", label="Stable")
plt.plot(Q_us/Q0, T_us, color="tab:orange", label="Unstable")
plt.plot(Q_s[n:]/Q0, T_s[n:], color="tab:purple")
plt.xlabel("Q/Q0")
plt.ylabel("T [K]")
plt.grid()
plt.legend()
plt.savefig('Q_Bifurcation.png')


#Computing the epsilon bifurcation diagram
T_s, T_us, eps_s, eps_us = compute_bifurcation(T_list, eps_list)

#Compute some stuff for colors in the plot
m = np.argmin(np.abs(T_s - T_us[0]))
n = np.argmin(np.abs(T_s - T_us[-1]))
eps_us = np.concatenate(([eps_s[m]], eps_us, [eps_s[n]]))
T_us = np.concatenate(([T_s[m]], T_us, [T_s[n]]))

#Plot epsilon bifurcation diagram
plt.figure()
plt.title("Epsilon bifurcation diagram energy balance model")
plt.plot(eps_s[:m], T_s[:m], color="tab:purple", label="Stable")
plt.plot(eps_us, T_us, color="tab:orange", label="Unstable")
plt.plot(eps_s[n:], T_s[n:], color="tab:purple")
plt.xlabel("Epsilon")
plt.ylabel("T [K]")
plt.grid()
plt.legend()
plt.savefig('eps_Bifurcation.png')






