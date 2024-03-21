#Import packages
import numpy as np
import matplotlib.pyplot as plt

#Define xdot = f
def f(x, l):
    return l-np.abs(1-x)*x

#Euler–Maruyama scheme
def euler(x0, sigma):
    x = np.empty(n)
    x[0] = x0
    for i in range(n-1):
        R = np.random.normal(0,1)
        x[i+1] = x[i] + f(x[i],l)*dt + sigma*np.sqrt(dt)*R
    return x

#Numerical Constants
T = 10 #Time for which we want to run the simulation
dt = 1e-3 #Timestep size
n = int(T/dt) #Amount of steps in the simulation
sigma = 0.2 #Amplitude of the noise

#Initial conditions
l = 0.15 #Take between 0 and 0.4
x0 = 0.4 #Take between 0 and 1.4

#Compute some solutions of the sde
x_noiseless_1 = euler(0.4, 0)
x_noiseless_2 = euler(0.9, 0)
x_noise_1 = euler(x0, sigma)
x_noise_2 = euler(x0, sigma)
t = np.linspace(0, T, n)

#Plot the results
plt.figure()
plt.plot(t,x_noiseless_1, label="x0=0.4, sigma=0")
plt.plot(t,x_noiseless_2, label="x0=0.9, sigma=0")
plt.plot(t,x_noise_1, label="x0="+str(x0)+", sigma="+str(sigma))
plt.plot(t,x_noise_2, label="x0="+str(x0)+", sigma="+str(sigma))
plt.title("Solutions SDE 2-box model starting at lambda = "+str(l))
plt.xlabel("t")
plt.ylabel("x")
plt.grid()
plt.legend()
plt.show()


