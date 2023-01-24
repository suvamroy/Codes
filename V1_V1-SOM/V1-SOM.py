########### V1 inspired Self-Organizing Map #############

# Load modules:

import numpy as np
from math import pi,exp

np.random.seed(0)

# SOM grid size (nxn):

n=10 

# Load dataset (should be reshaped into 2D numpy array):

x=np.load('x.npy')

# Initialize SOM grid:

N,p=x.shape[0],x.shape[1]
w=np.random.uniform(0,1, size=(n,n,p))

# Set parameter values: 

eta_0,gamma_0,sigma_0=0.5,n/2,1.33
K_0,K_1,K_2=0.024,0.024,0.024


# Loop over epochs:

for t in range (0,100):

	eta_t=eta_0*exp(-K_0*t)
	gamma_t=gamma_0*exp(-K_1*t)
	sigma_t=sigma_0*exp(-K_2*t)
	
	# Loop over number of input data points:
	
	for I in range (0,N):
	
		# Winner node:
		
		distSq = (np.square(w - x[I])).sum(axis=2)
		r,s = np.unravel_index(np.argmin(distSq), distSq.shape)
		
		# Update weights:
		
		w_o=w
		for i in range (0,n):
			for j in range (0,n):
			
				# Gaussian term dependent on winner node:
				
				d2=exp(-((i-r)**2+(j-s)**2)/(2*gamma_t**2))
				
				# Neighborhood influence:
				
				d1_s=0
				m=np.zeros((p),dtype=float)
				
				for k in range (0,n):
					for l in range (0,n):
					
						d1=exp(-((i-k)**2+(j-l)**2)/(2*sigma_t**2))
						d1_s+=d1
						m=m+(x[I]-w_o[k,l])*d1
				
				# Update:		
						
				w[i,j]=w[i,j]+eta_t*d2*m/d1_s
	
	# Export updated weights:
				
	np.save('wf',w)
