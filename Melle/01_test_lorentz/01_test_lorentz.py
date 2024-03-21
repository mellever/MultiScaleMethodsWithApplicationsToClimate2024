#Import packages
import numpy as np
import matplotlib.pyplot as plt

#Class that implements the zero-one test for the L63 system and the logistic map.
class zero_one_test():

    #Initialization function
    def __init__(self, mu=3.97, sigma=10, rho=28, beta=8/3, steps=1e4, dt=0.01, x0=[0, 1, 1.05], tau=1, c=np.pi):
        self.steps = int(steps)
        self.dt = dt
        self.x0 = x0
        self.tau = tau
        self.sigma = sigma
        self.rho = rho
        self.beta = beta
        self.mu = mu
        self.c = c

    #Function that defines the L63 system
    def L63(self, x):
        x,y,z = x
        return np.array([self.sigma*(y-x), x*(self.rho-z)-y, x*y-self.beta*z])

    #Function that solves the L63 system using Euler forward
    def solveL63(self):
        #Create array for saving the solution
        self.sol = np.empty((self.steps+1,3))
        
        #Set initial condition
        self.sol[0] = self.x0

        #Euler forward scheme
        for i in range(self.steps):
            self.sol[i+1] = self.sol[i] + self.L63(self.sol[i])*self.dt

        self.x = self.sol[:,0]
        self.y = self.sol[:,1]
        self.z = self.sol[:,2]

    #Function that solves the logistic map
    def solveLM(self):
        x0 = np.random.rand() #Choose random initial condition
        self.sol = np.empty(self.steps+1)
        self.sol[0] = x0
        for i in range(self.steps):
            self.sol[i+1] = self.mu*self.sol[i]*(1- self.sol[i])
        self.x = self.sol

    #Function for plotting the L63 system in 3D
    def plot3D(self):
        try: 
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            ax.plot(self.x, self.y, self.z, linewidth=0.3)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            ax.set_title('L63 System');
        except: print("Try solving the system first")

    #Function for plotting the timeseries of a coordinate
    def plotTimeseries(self, timeseries='x', sample=False):
        if timeseries == 'x': phi = self.x
        elif timeries == 'y': phi = self.y
        elif timeseries == 'z': phi = self.z
        else: "Input valid timeseries, choices are 'x', 'y' and 'z'"

        #Create time for timeseries
        try:
            temp = self.y[0]
            t = np.linspace(0, self.steps*self.dt, self.steps+1)
            system = 'L63 system'
            lw = 0.5
        except:
            t = np.linspace(0, self.steps, self.steps+1)
            system = 'logistic map'
            lw = 0.1

        #Plot the results
        fig = plt.figure()
        ax = plt.axes()
        ax.plot(t, phi, label="Numerical solution", linewidth=lw)
        ax.set_xlabel('t')
        ax.set_ylabel(timeseries)
        ax.set_title('Timeseries of '+timeseries+ ' of '+system);

        #Show sampled results
        if sample:
            j = int(1/self.tau)
            phi_sample = phi[(j-1)::j]
            t_sample = t[(j-1)::j]
            ax.scatter(t_sample, phi_sample, c='r', label="sample", marker='.')
            ax.legend()


    #Function that implements the 01-test 
    def test(self, timeseries='x', print_result=True):
        #Extract timeseries phi from solution
        if timeseries == 'x': phi = self.x
        elif timeseries == 'y': phi = self.y
        elif timeseries == 'z': phi = self.z
        else: "Input valid timeseries, choices are 'x', 'y' and 'z'"

        #Sample the solution
        l = int(1/self.tau)
        phi = phi[(l-1)::l]
        N = len(phi)
        Nc = N/10 #choose Nc <= N/10 by https://arxiv.org/pdf/0906.1418.pdf  
        self.n_list = np.arange(2, Nc)
        
        #Define auxillary process
        self.p = np.zeros(N+1)
        self.q = np.zeros(N+1)
        for i in range(len(phi)):
            self.p[i+1] = self.p[i] + np.cos(self.c*i)*phi[i] 
            self.q[i+1] = self.q[i] + np.sin(self.c*i)*phi[i] 

        #Mean square average displacement
        M = np.empty_like(self.n_list)
        for j in range(len(self.n_list)):
            n = int(self.n_list[j])
            M[j] =  np.sum(np.power(self.p[n:] - self.p[:N-n+1],2) + np.power(self.q[n:] - self.q[:N-n+1],2))/(self.steps)
        
        #Compute K using covariance method
        self.K = (np.cov(M, self.n_list)/np.sqrt(np.var(M)*np.var(self.n_list)))[0,1]

        #Print the results
        if print_result:
            if round(self.K, 0) == 1: print("K = ", round(self.K,4), " Chaotic dynamics")
            elif round(self.K, 0) == 0: print("K = ", round(self.K,4), " Regular dynamics")
            else: "01-test is not conclusive. look into sampling, timesteps or number of steps"

    #Function that plots the results of the auxillary system
    def auxillary_plotter(self):
        try:
            fig = plt.figure()
            ax = plt.axes()
            ax.plot(self.p, self.q, linewidth=0.5)
            ax.set_xlabel('p')
            ax.set_ylabel('q')
            ax.set_title('Auxillary system plot')
        except: print("Try performing the chaos test first")

    #Function that plots Kc for different values of c
    def k_plotter(self):
        c_list = np.linspace(0,2*np.pi, 100)
        K_list = np.empty_like(c_list)
        for i in range(len(c_list)):
            self.test(print_result=False)
            K_list[i] = self.K

        fig = plt.figure()
        ax = plt.axes()
        ax.scatter(c_list, K_list)
        ax.set_xlabel('c')
        ax.set_ylabel('Kc')
        ax.set_title('Kc versus c')
        ax.set_xlim(c_list[0], c_list[-1])
        ax.hlines(1, c_list[0], c_list[-1], colors='red', label="K=1")
        ax.hlines(0, c_list[0], c_list[-1], colors='green', label="K=0")
        ax.legend()

#Function that runs the model for different examples
def run(example): 
    #Choose fixed c uniformly from (pi/5,4pi/5) -> these values have been taken from https://arxiv.org/pdf/0906.1418.pdf
    eps = 1e-5
    c = np.random.uniform(np.pi/5+eps,4*np.pi/5)

    if example == example1:
        #Example of chaotic L63 system
        k = zero_one_test(c=c, tau=0.1)
        k.solveL63()
        k.plotTimeseries(sample=True)
        k.plot3D()
        k.test()
        k.auxillary_plotter()
        k.k_plotter()
        plt.show()
    elif example == example2:
        #Example of regular L63 system
        k = zero_one_test(rho=14, c=c, tau=0.1)
        k.solveL63()
        k.plotTimeseries(sample=True)
        k.plot3D()
        k.test()
        k.auxillary_plotter()
        k.k_plotter()
        plt.show()
    elif example == example3:
        #Example of chaotic logistic map
        k = zero_one_test(c=c)
        k.solveLM()
        k.plotTimeseries()
        k.test()
        k.auxillary_plotter()
        k.k_plotter()
        plt.show()

    elif example == example4:
        #Example of regular logistic map
        k = zero_one_test(mu=3.55, c=c)
        k.solveLM()
        k.plotTimeseries()
        k.test()
        k.auxillary_plotter()
        k.k_plotter()
        plt.show()
    else: print("System not recognized")


#_____Code for running the model_____
example1 = 'chaotic L63'
example2 = 'regular L63'
example3 = 'chaotic LM'
example4 = 'regular LM'

run(example1)