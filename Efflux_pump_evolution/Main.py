"""
Evolutionary Model of Pseudomonas aeruginosa Efflux Pump Regulatory Network
under Antibiotic Treatment

This script models the evolution of efflux pump gene regulatory networks
in Pseudomonas aeruginosa under antibiotic selection pressure. It simulates
protein expression dynamics using ODEs and incorporates mutation events
(gene loss/duplication) that are selected based on fitness costs.

Author: Dr. Suvam Roy
Institution: Umeå University
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Command-line arguments for parallel execution
# source: iterates through 25 different sets of regulation rates
# target: 25 trials for each regulation rate set
source = sys.argv[-1]
target = sys.argv[-2]

def RK4(f, y0, t):
    """
    4th Order Runge-Kutta method for solving systems of ODEs
    
    Parameters:
    -----------
    f : function
        Right-hand side of the ODE system (dy/dt = f(y,t))
    y0 : array
        Initial conditions
    t : array
        Time points for integration
        
    Returns:
    --------
    y : array
        Solution array with dimensions [n_variables, n_timepoints]
    """
    n = len(y0)
    y = np.zeros((n, len(t)), dtype=float)
    y[:, 0] = y0
    
    for i in range(0, len(t) - 1):
        # Apply antibiotic treatment at t = 30 hours
        # Allow system to reach steady state before antibiotic exposure
        if i * h == 30:
            # Add antibiotics to the system (indices 31-33 represent antibiotics)
            y[n-3, i] = y[n-3, i] + 1000  # Tobramycin
            y[n-2, i] = y[n-2, i] + 1000  # Ciprofloxacin  
            y[n-1, i] = y[n-1, i] + 1000  # Meropenem
        
        # Runge-Kutta 4th order method
        F1 = h * f(y[:, i], t[i])
        F2 = h * f((y[:, i] + F1 / 2), (t[i] + h / 2))
        F3 = h * f((y[:, i] + F2 / 2), (t[i] + h / 2))
        F4 = h * f((y[:, i] + F3), (t[i] + h))
        y[:, i + 1] = y[:, i] + 1 / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        
        # Ensure concentrations cannot go negative
        for j in range(n):
            if y[j, i + 1] < 0:
                y[j, i + 1] = 0
                
    return y

def f(y, t):
    """
    ODE system representing the rate of change of protein and antibiotic concentrations
    
    Variables:
    ----------
    y[0:10] : Protein concentrations of 4 RND efflux pumps 
               (3 pumps with 3 genes each, 1 pump with 2 genes)
    y[11:30] : Protein concentrations of 20 regulatory genes
    y[31:33] : Antibiotic concentrations (tobramycin, ciprofloxacin, meropenem)
    
    Parameters:
    -----------
    d : Regulation matrix (activation/inhibition rates)
    g : Gene copy number states (1=normal, 0=deletion, 2=duplication)
    k : Production rate constant
    
    Returns:
    --------
    y_dot : Array of time derivatives for all variables
    """
    # Calculate negative regulation effects (inhibition)
    neg = np.multiply(np.multiply(-d * (d < 0), g), y[0:n-3]) + np.ones((n-3, n-3), dtype=float)
    
    # Calculate positive regulation effects (activation)
    pos1 = np.multiply(np.multiply(d * (d > 0), g), y[0:n-3]) + np.ones((n-3, n-3), dtype=float)
    pos2 = 2 * np.multiply(np.multiply(d * (d > 0), g), y[0:n-3]) + np.ones((n-3, n-3), dtype=float)
    
    # Initialize derivatives array
    y_dot = np.zeros((n,), dtype=float)
    
    # Protein production and degradation dynamics
    # Includes regulatory network effects and basal degradation rate (0.3)
    y_dot[0:n-3] = np.multiply(
        k * np.divide(
            np.divide(np.prod(pos2, axis=1), np.prod(pos1, axis=1)), 
            np.prod(neg, axis=1)
        ), g
    ) - 0.3 * y[0:n-3]
    
    # Antibiotic dynamics with efflux pump effects
    # Model accounts for antibiotic removal by efflux pumps
    if y[3] > 0:  # If MexXY genes are present
        y_dot[31] = y_dot[31] - 0.1 * y[3] * (y[0] / (y[0] + y[3]))  # Tobramycin efflux
        y_dot[32] = y_dot[32] - 0.1 * (y[0] * (y[0] / (y[0] + y[3])) + y[5] + y[8])  # Ciprofloxacin efflux
        y_dot[33] = y_dot[33] - 0.1 * y[0] * (y[0] / (y[0] + y[3]))  # Meropenem efflux
    else:
        y_dot[31] = y_dot[31]  # No change to tobramycin
        y_dot[32] = y_dot[32] - 0.1 * (y[0] + y[5] + y[8])  # Ciprofloxacin efflux
        y_dot[33] = y_dot[33] - 0.1 * y[0]  # Meropenem efflux
    
    # Ensure derivatives don't drive concentrations negative
    for i in range(n):
        if y[i] < 0:
            y_dot[i] = 0
            
    return y_dot

# System dimensions and initialization
n = 34  # Total number of variables in the system
# Variable indices:
# 0-10: Efflux pump proteins (11 components for 4 pumps)
# 11-30: Regulatory proteins (20 genes)
# 31-33: Antibiotics (tobramycin, ciprofloxacin, meropenem)

y0 = [0] * n  # Initial conditions (all concentrations start at 0)
g = np.transpose([1] * (n - 3))  # Gene copy numbers (1 = normal, 0 = deletion, 2 = duplication)
k = 60  # Protein production rate constant

# Time parameters for simulation
t0 = 0  # Start time (hours)
T = 48  # End time (hours)
t = np.linspace(t0, T, 4801)  # Time points (0.01 hour resolution)
h = t[2] - t[1]  # Time step size

# Load regulation network parameters
# d contains activation/inhibition rates between network components
d = np.load('nets/Nets/' + source + 'net.npy')

# Regulatory genes that can undergo mutation
gene = list(range(11, 31))  # Indices 11-30 represent regulatory proteins

# Names of regulatory genes for reference
regulatory_genes = ['armR', 'nalC', 'nalD', 'rocS2', 'mexR', 'brlR', 'mexT', 
                   'mexZ', 'armZ', 'rplU', 'suhB', 'parRS', 'amgRS', 'nfxB', 
                   'algU', 'vqsM', 'esrC', 'mexS', 'ampR', 'mvaT']

# Initial simulation to establish baseline
y = RK4(f, y0, t)

# Calculate fitness costs
# c1: Protein production cost (total protein concentration over a time period)
# c2: Antibiotic cost (internal antibiotic concentration over a time period)
# Integration from t=3000 (30 hours) to end captures post-antibiotic dynamics
c1 = (h * np.sum(y[0:11, 3000:])) + (h * np.sum(y[11:31, 3000:]))  # Protein cost
c2 = (h * np.sum(y[31:n, 3000:]))  # Antibiotic cost
c1o = c1  # Store initial protein cost for comparison

# Evolutionary simulation
Gen = [0]  # Generation counter
Cp = [c1]  # Protein cost trajectory
Ca = [c2]  # Antibiotic cost trajectory

for gen in range(0, 100):
    # Stop condition: if too many genes have mutated (deletion or duplication)
    if list(g).count(0) + list(g).count(2) >= 4:
        break

    # Run simulation with current genotype
    y = RK4(f, y0, t)
    
    # Calculate current costs
    c1 = (h * np.sum(y[0:11, 3000:])) + (h * np.sum(y[11:31, 3000:]))
    c2 = (h * np.sum(y[31:n, 3000:]))
    
    # Randomly select a regulatory gene for mutation
    mut = np.random.choice(gene, 1)[0]
    
    # Determine mutation type with equal probability
    if np.random.uniform(0, 1) < 0.5:
        g[mut] = 0  # Gene deletion
        mutation = 'loss'
    else:
        g[mut] = 2  # Gene duplication  
        mutation = 'gain'
    
    # Test fitness of mutant
    y = RK4(f, y0, t)
    
    # Calculate mutant costs
    c11 = (h * np.sum(y[0:11, 3000:])) + (h * np.sum(y[11:31, 3000:]))
    c21 = (h * np.sum(y[31:n, 3000:]))
    
    # Reset gene to normal state for evaluation
    g[mut] = 1
    
    # Fitness criterion: mutation is accepted if:
    # 1. Antibiotic cost decreases (c21 < 0.99*c2) AND
    # 2. Protein cost remains within 2x initial cost (c11 < 2*c1o)
    if c11 < 2 * c1o and c21 < 0.99 * c2:
        if mutation == 'loss':
            g[mut] = 0  # Accept deletion
        if mutation == 'gain':
            g[mut] = 2  # Accept duplication
            
        gene.remove(mut)  # Remove from pool of mutable genes
        Gen.append(gen + 1)  # Record generation
        Cp.append(c11)  # Record protein cost
        Ca.append(c21)  # Record antibiotic cost

# Save results
data = np.array([Gen, Cp, Ca])
final_target = str(int(source) * 25 + int(target))

# Save evolved genome and time series data
np.save('./tob/' + final_target + 'genome', g)
np.save('./tob/' + final_target + 'ts', data)
