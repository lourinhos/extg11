# Genetics Computer Simulation Exercise Computer Code
#
# EXTG11, EXTG15: Spring 2025
#
# Below I have provided several functions you will 
# need to complete your computer lab exercises.
#
# You should be able to run this code using whatever
# python interface you prefer (e.g., pasting directly
# into an interactive session using the terminal, or 
# using a Jupyter notebook). Either way, be sure that
# you have the correct packages installed. You may 
# need to create a python virtual environment to do so.

# Import necessary packages
from scipy import *
import matplotlib.pyplot as plt
import numpy as np
import random 
import os 
import sys


### Deterministic 1-locus selection model: recursion equation
# args:
# q   : Frequency of A2 allele  (numeric, between 0,1)
# w_1 : Relative fitness of A1A1 genotype (numeric, generally somewhere between 0 & 1.5)
# w_2 : Relative fitness of A1A2 genotype (numeric, generally somewhere between 0 & 1.5)
# w_3 : Relative fitness of A2A2 genotype (numeric, generally somewhere between 0 & 1.5)
# u   : Mutation rate from A1 --> A2 (numeric, generally quite small, e.g., 0.0001)
def qNext(q, w_1, w_2, w_3, u):
    q = q + (1 - q)*u
    return (q**2 * w_3 + (1 - q)*q*w_2) / ((1 - q)**2 * w_1 + 2*(1 - q)*q * w_2 + q**2 * w_3)


### Function to generate time series of deterministic allele frequencies
# args: passes same arguments to qNext(), with one additional
# q0   : Initial frequency of A2 allele (numeric, between 0,1)
# w_1  : Relative fitness of A1A1 genotype (numeric, generally somewhere between 0 & 1.5)
# w_2  : Relative fitness of A1A2 genotype (numeric, generally somewhere between 0 & 1.5)
# w_3  : Relative fitness of A2A2 genotype (numeric, generally somewhere between 0 & 1.5)
# u    : Mutation rate from A1 --> A2 (default value of 0, generally quite small, e.g., 0.0001)
# nGen : Number of generations to interate (positive integer value.)
def determTraj(q0, w_1, w_2, w_3, nGen, u=0):

    q = [0] * nGen # storage array
    q[0] = q0 # initial frequency
    # loop over time steps to calculate allele frequency trajectory
    for t in range(1,nGen):
        q[t] = qNext(q=q[t-1], w_1=w_1, w_2=w_2, w_3=w_3, u=u)

    # return allel freq. trajectory
    return np.array(q)


### Function to plot deterministic trajectory
# args: passes same arguments to determTraj(), with two additional
# q0   : Initial frequency of A2 allele (numeric, between 0,1)
# w_1  : Relative fitness of A1A1 genotype (numeric, generally somewhere between 0 & 1.5)
# w_2  : Relative fitness of A1A2 genotype (numeric, generally somewhere between 0 & 1.5)
# w_3  : Relative fitness of A2A2 genotype (numeric, generally somewhere between 0 & 1.5)
# u    : Mutation rate from A1 --> A2 (default value of 0, generally quite small, e.g., 0.0001)
# nGen : Number of generations to interate (positive integer value.)
# plotP: Plot the trajectory for A1 allele instead? (logical, True/False)
# yLim : Alternative y-axis limits (list with two values, default is [0,1])
def plotDetermTraj(q0, w_1, w_2, w_3, nGen, u=0, plotP=False, yLim=[0,1]):
    
    # make data
    q = determTraj(q0 = q0, w_1 = w_1, w_2 = w_2, w_3 = w_3, u = u, nGen = nGen)
    time = np.arange(start=0, stop=nGen, step=1)

    # plot options
    fig,ax=plt.subplots()
    ax.set_xlim([0, nGen])
    ax.set_ylim(yLim)
    ax.set_xlabel('Generations')
    ax.set_ylabel('Frequency of $A_2$ allele ($q$)')

    # make plot
    if plotP:
        plt.plot(time, 1-q, color='dodgerblue', ls='-', lw=2)
    else:
        plt.plot(time, q, color='dodgerblue', ls='-', lw=2)

    plt.show()


### Function to plot per-generation change in allele frequency
# args: passes same arguments to determTraj(), with two additional
# w_1  : Relative fitness of A1A1 genotype (numeric, generally somewhere between 0 & 1.5)
# w_2  : Relative fitness of A1A2 genotype (numeric, generally somewhere between 0 & 1.5)
# w_3  : Relative fitness of A2A2 genotype (numeric, generally somewhere between 0 & 1.5)
# u    : Mutation rate from A1 --> A2 (default value of 0, generally quite small, e.g., 0.0001)
# plotP: Plot the trajectory for A1 allele instead? (logical, True/False)
# yLim : Alternative y-axis limits (list with two values, default is [0,1])
def plotDeltaQ(w_1, w_2, w_3, u=0, yLim=[-0.1,0.1], plotP=False, plotPHat=True):
    
    # make data
    q = np.linspace(0, 1, 200)
    deltaQ = [0] * np.size(q)
    for i in range(0,np.size(q)):
        deltaQ[i] = qNext(q=q[i], w_1=w_1, w_2=w_2, w_3=w_3, u=u) - q[i]

    # Calculate equilibrium    
    if w_1 == w_3:
        pHat = 1/2
    else:
        h = (w_2-1)/(w_3-1)
        pHat = (h - 1)/(2*h - 1)

    # plot options
    fig,ax=plt.subplots()
    ax.set_xlim([0,1])
    ax.set_ylim(yLim)

    # make plot
    if plotP:
        deltaP = deltaQ
        deltaP = [x * -1 for x in deltaP]
        p = 1 - q
        ax.set_xlabel('Frequency of $A_1$ allele ($p$)')
        ax.set_ylabel('Per-generation change in frequency ($\\Delta p$)')
        plt.plot(p, deltaP, color='dodgerblue', ls='-', lw=2)
    else:
        ax.set_xlabel('Frequency of $A_2$ allele ($q$)')
        ax.set_ylabel('Per-generation change in frequency ($\\Delta q$)')
        plt.plot(q, deltaQ, color='dodgerblue', ls='-', lw=2)

    plt.axhline(y = 0, color = 'black', linestyle = '-')

    if(plotPHat):
        if(0< pHat < 1):
            plt.axvline(x = pHat, color = 'red', linestyle = '--') 
    
    plt.show()

### Plot multiple deterministic trajectories together - helps visualize equilibria
# args: Passes same arguments to plotDetermTraj(), with two additional
# q0   : Initial frequency of A2 allele (numeric, between 0,1)
# w_1  : Relative fitness of A1A1 genotype (numeric, generally somewhere between 0 & 1.5)
# w_2  : Relative fitness of A1A2 genotype (numeric, generally somewhere between 0 & 1.5)
# w_3  : Relative fitness of A2A2 genotype (numeric, generally somewhere between 0 & 1.5)
# u    : Mutation rate from A1 --> A2 (default value of 0, generally quite small, e.g., 0.0001)
# nGen : Number of generations to interate (positive integer value.)
# n_lines : Number of lines to plot (positive integer; 20 usually works well)
# clrmp : Colormap for lines (accepts matplotlib default named colormaps)
def plotDetermTrajMulti(w_1, w_2, w_3, nGen, u = 0, n_lines=20, clrmp = 'managua', plotP=False, plotPHat = False):
    cmap = plt.colormaps[clrmp]
    # Take colors at regular intervals spanning the colormap.
    colors = cmap(np.linspace(0, 1, n_lines))

    # Calculate equilibrium    
    if w_1 == w_3:
        pHat = 1/2
    else:
        h = (w_2-1)/(w_3-1)
        pHat = (h - 1)/(2*h - 1)
    
    fig,ax=plt.subplots()
    ax.set_xlim([0, nGen])
    ax.set_ylim([0, 1])
    ax.set_xlabel('Generations')
    if plotP:
        ax.set_ylabel('Frequency of $A_1$ allele ($p$)')
    else:
        ax.set_ylabel('Frequency of $A_2$ allele ($q$)')
            
    time = np.arange(start=0, stop=nGen, step=1)
    qInits = np.linspace(0, 1, n_lines)
    for i in range(n_lines):
        q = determTraj(q0 = qInits[i], w_1 = w_1, w_2 = w_2, w_3 = w_3, u = u, nGen = nGen)
        ax.plot(time, q, color=colors[i], ls='-', lw=2)
        #plt.pause(0.01)

    if(plotPHat):
        if(0< pHat < 1):
            plt.axhline(y = 1 - pHat, color = 'red', linestyle = '--') 
    
    plt.show()





### Wright-Fisher simulation for specified number of generations
##  Convenient for plotting, time axis is always the same
# args: 
# q0   : Initial frequency of A2 allele (numeric, between 0,1)
# w_1  : Relative fitness of A1A1 genotype (numeric, generally somewhere between 0 & 1.5)
# w_2  : Relative fitness of A1A2 genotype (numeric, generally somewhere between 0 & 1.5)
# w_3  : Relative fitness of A2A2 genotype (numeric, generally somewhere between 0 & 1.5)
# u    : Mutation rate from A1 --> A2 (default value of 0, generally quite small, e.g., 0.0001)
# nGen : Number of generations to interate (positive integer value.)
# N    : Population size (positive integer)
def WFSim(q0, w_1, w_2, w_3, u, nGen, N):

    q    = [0] * nGen
    q[0] = q0
    t    = np.arange(0,nGen,1)

    # Loop over number of generations
    for i in range(1,nGen):
        q[i-1] = q[i-1] + q[i-1]*u
        p = 1 - q[i-1]
        w_avg = (p**2 * w_1 + 2 * p * q[i-1] * w_2 + q[i-1]**2 * w_3)
        
        P_ij_det = np.array([p**2 * w_1 ,
                             2 * p * q[i-1]* w_2,
                             q[i-1]**2 * w_3]) / w_avg
        
        # Simulating multinomial distribution using numpy
        P_ij_drift = np.random.multinomial(int(round(N)), P_ij_det) / round(N)
        
        # Update q
        q[i] = P_ij_drift[2] + P_ij_drift[1] / 2

    return np.array([q,t])


# Function plotting many WF simulations. Good for visualizing outcomes, 
# strength of drift vs. selection
# args: Passes same arguments to WFSim(), with two additional
# q0   : Initial frequency of A2 allele (numeric, between 0,1)
# w_1  : Relative fitness of A1A1 genotype (numeric, generally somewhere between 0 & 1.5)
# w_2  : Relative fitness of A1A2 genotype (numeric, generally somewhere between 0 & 1.5)
# w_3  : Relative fitness of A2A2 genotype (numeric, generally somewhere between 0 & 1.5)
# u    : Mutation rate from A1 --> A2 (default value of 0, generally quite small, e.g., 0.0001)
# nGen : Number of generations to interate (positive integer value.)
# N    : Population size (positive integer)
# n_lines : Number of lines to plot (positive integer; 20 usually works well)
# clrmp : Colormap for lines (accepts matplotlib default named colormaps)
def plotWFSims(q0, w_1, w_2, w_3, nGen, N, n_lines, u = 0, clrmp = 'vanimo'):
    cmap = plt.colormaps[clrmp]
    # Take colors at regular intervals spanning the colormap.
    colors = cmap(np.linspace(0, 1, n_lines))
    
    fig,ax=plt.subplots()
    ax.set_xlim([0, nGen])
    ax.set_ylim([0, 1])
    ax.set_xlabel('Generations')
    ax.set_ylabel('Frequency of $A_2$ allele ($q$)')
    
    for _ in range(n_lines):
        tmpWF = WFSim(q0 = q0, w_1=w_1, w_2=w_2, w_3=w_3, u=u, nGen=nGen, N=N)
        if w_1 < w_2 < w_3:
            ind = int(np.round(np.max(tmpWF[0])*n_lines, 2)-1)
            plt.plot(tmpWF[0], color=colors[ind], ls='-', lw=2)
        else:
            plt.plot(tmpWF[0], color=colors[random.randrange(n_lines)], ls='-', lw=2)
        #plt.pause(0.05)
    
    plt.show()

# Analytic expression for expected heterozygosity decay over timed due to genetic drift
# H0: initial heterozygosity at time = 0
# N:  Population size
# t: number of generations to plot
def H_Hat(H0, N, t):
    return(H0*(1 - 1/(2*N))**t)

# Function to plot decay of heterozygosity over time, comparing the analytic expectation with simulated data
# q0   : Initial frequency of A2 allele (numeric, between 0,1)
# w_1  : Relative fitness of A1A1 genotype (numeric, generally somewhere between 0 & 1.5)
# w_2  : Relative fitness of A1A2 genotype (numeric, generally somewhere between 0 & 1.5)
# w_3  : Relative fitness of A2A2 genotype (numeric, generally somewhere between 0 & 1.5)
# u    : Mutation rate from A1 --> A2 (default value of 0, generally quite small, e.g., 0.0001)
# nGen : Number of generations to interate (positive integer value.)
# N    : Population size (positive integer)
# n_lines : Number of lines to plot (positive integer; 20 usually works well)
def H_decay_plot(q0, N, w_1, w_2, w_3, u, nGen, n_sims):
    # Data storage matrix
    dat = np.full((n_sims, nGen), np.nan)
    time = np.arange(nGen)
    
    # Loop over number of replicate simulations, calculate H over time
    for i in range(n_sims):
        # Assuming WFSim returns a dictionary with a key 'q' (similar to the R code)
        q_t = WFSim(q0=q0, w_1=w_1, w_2=w_2, w_3=w_3, u=u, nGen=nGen, N=N)[0]
        dat[i, :] = 2 * q_t * (1 - q_t)
    
    H_bar = np.mean(dat, axis=0)  # calculate mean Heterozygosity over time for all simulations
    H_exp = H_Hat(H0=2*q0*(1-q0), N=N, t=time)
    
    # Plot
    plt.figure(figsize=(8, 6))
    plt.plot(time, H_bar, color="red", linewidth=2, label="Mean Heterozygosity")
    plt.plot(time, H_exp, color="dodgerblue", linewidth=2, label="Expected Heterozygosity")
    plt.xlim(0, nGen)
    plt.ylim(0, 0.55)
    plt.xlabel("Generations")
    plt.ylabel("Heterozygosity")
    plt.legend()
    plt.show()


### Wright-Fisher simulation using a 'while' loop
##  Convenient for counting outcomes, not great for plotting
##  because time axis varies
# args: 
# q0   : Initial frequency of A2 allele (numeric, between 0,1)
# w_1  : Relative fitness of A1A1 genotype (numeric, generally somewhere between 0 & 1.5)
# w_2  : Relative fitness of A1A2 genotype (numeric, generally somewhere between 0 & 1.5)
# w_3  : Relative fitness of A2A2 genotype (numeric, generally somewhere between 0 & 1.5)
# N    : Population size (positive integer)
def WFSimToEnd(q0, w_1, w_2, w_3, N):

	# Define intial frequecy
    q = q0

    # Loop over number of generations
    t = 0
    while 0 < q*(1-q) < 1:
        p = 1 - q
        w_avg = (p**2 * w_1 + 2 * p * q * w_2 + q**2 * w_3)
        
        P_ij_det = np.array([p**2 * w_1 ,
                             2*p*q* w_2,
                             q**2 * w_3]) / w_avg
        
        # Simulating multinomial distribution using numpy
        P_ij_drift = np.random.multinomial(int(round(N)), P_ij_det) / round(N)
        
        # Update q, t
        q = P_ij_drift[2] + P_ij_drift[1] / 2
        t = t+1

    # return list of two numbers: final frequency of A2
	# and the duration of simulation
    return [int(q),t]


### Simple progress bar function to track progress in simulations
def progress(count, total, suffix=''):
    bar_len = 20
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', suffix))
    sys.stdout.flush()  # As suggested by Rom Ruben


### Wrapper function to run multiple WF simulations and count outcomes
# args: Passes same args to WFSimToEnd(), with one additional
# q0   : Initial frequency of A2 allele (numeric, between 0,1)
# w_1  : Relative fitness of A1A1 genotype (numeric, generally somewhere between 0 & 1.5)
# w_2  : Relative fitness of A1A2 genotype (numeric, generally somewhere between 0 & 1.5)
# w_3  : Relative fitness of A2A2 genotype (numeric, generally somewhere between 0 & 1.5)
# N    : Population size (positive integer)
# nSims: Number of replicate simulations to run (positive integer, usually >= 1000)
def WFSimCountFix(q0, w_1, w_2, w_3, N, nSims):
    fixations = int(0)
    doesA2Fix = []
    generations = []
    for i in range(nSims):
        simOut = WFSimToEnd(q0=q0, w_1=w_1, w_2=w_2, w_3=w_3, N=N)
        doesA2Fix.append(int(simOut[0]))
        fixations = fixations + int(simOut[0])
        generations.append(simOut[1])
        progress(i, nSims)

    print('\nReplicate simulations:')
    print(nSims)
    print('fixations:')
    print(fixations)
    print('Avg. time to extinction:')
    print(round(np.mean(np.array(generations)[np.array(doesA2Fix) == 0]),3))
    print('Avg. time to fixation:')
    print(round(np.mean(np.array(generations)[np.array(doesA2Fix) == 1]),3))



### Plotting function to plot decay of heterozygosity due to drift
# args: 
# H0   : Initial heterozygosity (numeric, between 0, 0.5)
# N    : Population size (positive integer)
# nGen : number of generations to plot
def plotDriftHeterozygosity(H0, N, nGen, yLim=[0,0.5]):
    
    # make data
    t  = np.arange(start=0, stop=nGen, step=1)
    Ht = H0*(1 - 1/(2*N))**t

    # plot options
    fig,ax=plt.subplots()
    ax.set_xlim([0, nGen])
    ax.set_ylim(yLim)
    ax.set_xlabel('Generations')
    ax.set_ylabel('Expected Heterozygosity ($H$)')

    # make plot
    plt.plot(t, Ht, color='dodgerblue', ls='-', lw=2)

    plt.show()
