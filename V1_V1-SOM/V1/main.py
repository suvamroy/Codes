###### Script to generate example of V1 action #######

# Load modules:

import numpy as np
import scipy.io
import os
from V1 import v1

# Load orientation selectivity matrices 'orr'
# and 2D Marr wavelet matrices 'wl':

orr1=np.load('orr_l4.npy')
w1l=np.load('w1l.npy')
orr2=np.load('orr_l23.npy')
w2l=np.load('w1l.npy')
orr3=np.load('orr_l5.npy')
w3l=np.load('w2l.npy')


# Load example LGN image from BSDS500 dataset:	
	
mat = scipy.io.loadmat('12003.mat')
mat=mat['groundTruth'][0][3]['Boundaries'][0][0]



N=max(mat.shape[0],mat.shape[1])	

lgn=np.zeros((N,N),dtype=float)

# Make images square shaped by padding

if mat.shape[0]>mat.shape[1]:
	lgn[:,80:N-80]=mat
else:
	lgn[80:N-80,:]=mat
del mat

# Crop image to even size
	
lgn=lgn[0:480,0:480]

# Export example LGN image

np.save('example_lgn',lgn)

# Apply V1 processing by calling the V1.py script
	
v1_o=v1(lgn,orr1,orr2,orr3,w1l,w2l,w3l)

# Export V1 output
	
np.save('example_v1',v1_o)

# Apply ReLu

for i in range (0,480):
	for j in range (0,480):
		if v1_o[i,j]<0.5:
			v1_o[i,j]=0
			
# Export V1 output (post ReLu)
			
np.save('example_v1_relu',v1_o)
	
	
	
	
