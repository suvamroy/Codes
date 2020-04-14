# Rate of evolution of TFT players in a finite population of ALLD

import numpy as np

x=0  					# times TFT gets fixed
N=100    				# Population size
w=1.0 					# Intensity of selection
Trials=1000				# Total trials
a,b,c,d=3.0,0.0,5.0,1.0			# 1 round Payoff matrix between ALLC and ALLD
m=10					# Number of rounds
for trial in xrange (0,Trials):   	
	y=1				# Initial number of TFT
    	q=[1]*(N-y)+[0]*y		# Population array (0 is ALLD, 1 is TFT)
    	t=0
	f=[0.0]*N			# Fitness array
    	while (0<y<N):

		# Payoff assignment after m round games :

		np.random.shuffle(q)
        	for n in xrange(0,N,2):
            		if (q[n]==0 and q[(n+1)]==0):
	        		f[n]=(1.0-w+w*m*a)
	       	 		f[n+1]=(1.0-w+w*m*a)
	    		if (q[n]==0 and q[(n+1)]==1):
	        		f[n]=(1.0-w+w*(b+(m-1)*d))
	        		f[n+1]=(1.0-w+w*(c+(m-1)*d))
	    		if (q[n]==1 and q[(n+1)]==0):
	        		f[n]=(1.0-w+w*(c+(m-1)*d))
	        		f[n+1]=(1.0-w+w*(b+(m-1)*d))  
	    		if (q[n]==1 and q[(n+1)]==1):
	        		f[n]=(1.0-w+w*m*d)
	        		f[n+1]=(1.0-w+w*m*d) 

		# Moran process :

		indices=range(0,N)
		i=np.random.choice(indices, 1, p=np.array(f)/sum(f))[0]
		del indices[i]
		j=np.random.choice(indices, 1)[0]
		indices=[]
		q[j]=q[i]
		y=(N-sum(q))  # Number of TFTs
       		t=t+1
    	if y==N:
		x+=1

print 'Rate of evolution, Np_TFT =',(1.0*x*N)/Trials
