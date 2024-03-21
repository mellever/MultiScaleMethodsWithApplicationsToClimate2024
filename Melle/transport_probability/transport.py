#Import packages
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

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
    V = []
    for z in z_arr:
        I, err = integrate.quad(lambda x: -f(x,l), 0, z)
        V.append(I)
    return np.array(V)

def rho(V, sigma):
    return np.exp(-2*V/np.power(sigma,2))

#Define solution
def rho_inf(sigma, l, z_arr, dz):
    V = compute_V(l, z_arr)
    C = np.sum(rho(V, sigma)*dz)
    return rho(V, sigma)/C

def solve(T, dt, N):
    #Numerical Constants
    n = int(T/dt) #Amount of steps in the simulation
    sigma = 0.2 #Amplitude of the noise
    dz = 1e-2 #Stepsize z

    #Initial conditions
    l = 0.1 #Take between 0 and 0.4
    z0 = 0.4 #Take between 0 and 1.4


    #Compute some solutions of the sde
    z = euler(z0, n, N, sigma, l, dt)
    t = np.linspace(0, T, n)

    z_arr = np.arange(np.min(z), np.max(z), dz)
    p_inf = rho_inf(sigma, l, z_arr, dz)

    times = [0, 1, 2, 100, 500, len(z)-1]
    data = z[times]
    
    if N == 1: return z, z_arr, p_inf
    else: return data, times, z_arr, p_inf



N = 1000

if N == 1:
    z, z_arr, p_inf = solve(5e3, 1e-2, N)
    counts, bins = np.histogram(z, density=True, bins=100)
    plt.stairs(counts, bins, label="Numerical")
else:
    z, times, z_arr, p_inf = solve(10, 1e-2, N)
    for i in range(len(z)):
        counts, bins = np.histogram(z[i], density=True, bins=50)
        plt.stairs(counts, bins, label="timestep = "+str(times[i]))
plt.plot(z_arr, p_inf, label=r'$\rho_\infty$', linewidth=5, c='r')
plt.xlim(np.min(z), np.max(z))
plt.xlabel("z")
plt.ylabel("Density")
plt.legend()
plt.show()
