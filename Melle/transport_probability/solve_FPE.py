#Import packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.sparse import csr_array

#Function of stommel 2box model
def f(x, l):
    return l-np.abs(1-x)*x

#Compute V
def compute_V(l, z_arr):
    V = -l*z_arr - np.sign(1-z_arr)*(np.power(z_arr,3)/3 - np.power(z_arr,2)/2 + 1/6)
    return V

#Rho
def rho(V, sigma):
    return np.exp(-2*V/np.power(sigma,2))

#Define solution
def rho_inf(sigma, l, z_arr, dz):
    V = compute_V(l, z_arr) #This can be done analytically
    C = np.sum(rho(V, sigma)*dz)
    return rho(V, sigma)/C


#Discretized Fokker-Plank operator
def L_FP(x, dx, l, sigma):
    #Number of x points
    nx = len(x)

    #Create arrays
    diag_f = np.diag(f(x,l))
    D1 = np.zeros((nx, nx))
    D2 = np.zeros_like(D1)

    #Discretize derivatives using central differences
    for i in range(nx):
        D2[i,i] = -2
        if i<len(D1)-1:
            D2[i+1,i] = 1
            D1[i+1,i] = -1
        if i>0:
            D2[i-1,i] = 1
            D1[i-1,i] = 1
    
    #Implement periodic boundary conditions
    D1[-1,0] = 1
    D1[0,-1] = -1
    D2[-1,0] = 1
    D2[0,-1] = 1

    #Scale by stepsize
    D1 = D1/(2*dx)
    D2 = D2/np.power(dx,2)

    #Return results
    return -np.dot(D1, diag_f) + D2*np.power(sigma,2)/2

def solve(x, dx, t, dt, rho0, l, sigma):
    #Compute discretized FP operator
    L = L_FP(x, dx, l, sigma)

    #Create array for saving the results
    rho = np.empty((len(t), len(x)))

    #Set initial distribution
    rho[0,:] = rho0

    #Solve for rho
    for n in range(len(t)-1):
        print("Progress = ", np.round(n/(len(t)-1)*100,3), "%")
        rho[n+1,:] = csr_array(expm(dt*L)).dot(rho[n,:]) #Compute rho using sparse matrix product

    #Return the result
    return rho


#Parameters
sigma = 0.2
l = 0.3
T = 1000
X = 5

#Discretization
dx = 5e-2
x = np.arange(-X, X, dx)
dt = 5
t = np.arange(0, T+dt, dt)

#Set initial distribution to delta distribution
x0 = 0.5
rho0 = np.zeros_like(x)
rho0[np.argmin(np.abs(x-x0))] = 1
C0 = np.sum(rho0*dx)
rho0 = rho0/C0

#Compute results
rho_num = solve(x, dx, t, dt, rho0, l, sigma)
rho_an = rho_inf(sigma, l, x, dx)

#Set some plot parameters
bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
stats = (f'$x_0$ = ' +str(x0) +'\n' f'$\\lambda$ = '+str(l)+'\n' f'$\\sigma$ = '+str(sigma))
t_indices = np.array([1,2,10,len(t)-1])
t_plot = t_indices*dt
fname = str(x0)+str(l)+str(sigma)+'.png'

#Plot the results
plt.figure(figsize=(19.20, 10.80), dpi=80)
plt.title("Numerical solution FPE 2-box model", fontsize=20)
plt.plot(x, rho0, label=r'$\rho_0$')
plt.plot(x, rho_an, label=r'$\rho_\infty$')
for i in range(len(t_indices)): plt.plot(x, rho_num[int(t_indices[i]),:], label="T="+str(round(t_plot[i],2)))
plt.xlabel("x")
plt.ylabel(r'$\rho$')
plt.grid()
plt.legend(fontsize=20)
plt.text(-1.9, 8.4, stats, bbox=bbox, fontsize=20)
plt.xlim(-1,2)
#plt.savefig("/home/melle/OneDrive/Master/Year1/MultiscaleMethods/transport_probability/"+fname, format='png')
plt.show()
