#Import packages
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.linalg import expm

#Define xdot = f
def f(x, l):
    return l-np.abs(1-x)*x

#Euler–Maruyama scheme for dz = -dV/dx dt+ sigma dW = f(x) dt + sigma dW
def euler(z0, n, N, sigma, l, dt):
    N = int(N)
    z = np.empty((n, N))
    z[0, :] = z0
    for i in range(n-1):
        R = np.random.normal(0,1,N)
        z[i+1, :] = z[i, :] + f(z[i],l)*dt + sigma*np.sqrt(dt)*R
    return z

#Compute V
def compute_V(l, z_arr):
    V = -l*z_arr - np.sign(1-z_arr)*(np.power(z_arr,3)/3 - np.power(z_arr,2)/2 + 1/6)
    return V

#Define rho
def rho(V, sigma):
    return np.exp(-2*V/np.power(sigma,2))

#Define solution
def rho_inf(sigma, l, z_arr, dz):
    V = compute_V(l, z_arr) #This can be done analytically 3
    C = np.sum(rho(V, sigma)*dz)
    return rho(V, sigma)/C

#Function that solves the SDE
def solve(T, dt, N):
    #Compute some solutions of the sde
    z = euler(z0, n, N, sigma, l, dt)
    t = np.linspace(0, T, n)

    z_arr = np.arange(np.min(z), np.max(z), dz)
    p_inf = rho_inf(sigma, l, z_arr, dz)

    times = [0, 1, 2, 100, 200, 300, len(z)-1]
    cmap = plt.get_cmap('hsv')
    colors = cmap(np.linspace(0,1,len(times)))
    data = z[times]
    
    if N == 1: return z, z_arr, p_inf
    else: return data, times, z_arr, p_inf, colors


#Numerical Constants
sigma = 0.2 #Amplitude of the noise
dz = 1e-2 #Stepsize z
N = 1 #Ensemble size

#Initial conditions
l = 0.1 #Take between 0 and 0.4
z0 = 0.4 #Take between 0 and 1.4


if N == 1:
    T = 1e3
    dt = 1e-2
    n = int(T/dt) #Amount of steps in the simulation
    z, z_arr, p_inf = solve(T, dt, N)
    plt.plot(z_arr, p_inf, label=r'$\rho_\infty$', linewidth=1, c='black')
    counts, bins = np.histogram(z, density=True, bins=100)
    plt.stairs(counts, bins, label="Numerical")
    plt.title("Probability density plot at T = "+str(T))

else:
    z, times, z_arr, p_inf, colors = solve(10, 1e-2, N)
    plt.plot(z_arr, p_inf, label=r'$\rho_\infty$', linewidth=2, c='black')
    for i in range(len(z)):
        counts, bins = np.histogram(z[i], density=True, bins=50)
        plt.stairs(counts, bins, color=colors[i], label="timestep = "+str(times[i]))
        plt.title("Probability density ensemble N = "+str(N))
plt.xlim(np.min(z), np.max(z))
plt.xlabel("z")
plt.ylabel("Density")
plt.legend()
plt.show()




