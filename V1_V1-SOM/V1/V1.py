import numpy as np

def v1(lgn,orr1,orr2,orr3,w1l,w2l,w3l):

	# Add padding to LGN image
	
	K=7	# Spatial spread of L4 to LGN afferent weights
	N=lgn.shape[0]
	lgn1=lgn
	lgn=np.zeros((N+2*K,N+2*K),dtype=float)
	lgn[K:N+K,K:N+K]=lgn1
	
	
	### Layer 4 ###
	
	l4=np.zeros((N,N),dtype=float)
	
	# Activity from afferent connections
	
	for i in range (0,N):
		for j in range (0,N):
		
			n_iso_domains=int(N/40)	# Size of iso-orientation domain = (40x40)
			
			# Select orientation selective afferent weight matrix
		
			i1,j1=int(i/40),int(j/40)
			if i1 in list(range(0,n_iso_domains+4,2)):
				if j1 in list(range(0,n_iso_domains+4,2)):
					w1=orr1[0,:,:]
				else:
					w1=orr1[3,:,:]
			else:
				if j1 in list(range(0,n_iso_domains+4,2)):
					w1=orr1[1,:,:]
				else:
					w1=orr1[2,:,:]
					
			# Generate L4 node activity
				
			l4[i,j]=np.sum(np.multiply(lgn[i:i+2*K+1,j:j+2*K+1],w1))
		
	del lgn
	
	# Activity from lateral connections
	
	K=75	# padding
	
	l4l=np.zeros((N+2*K,N+2*K),dtype=float)
		
	for i in range (0,N):
		for j in range (0,N):
		
			# Generate L4 node activity 
			
			l4l[i:i+2*K+1,j:j+2*K+1]=l4l[i:i+2*K+1,j:j+2*K+1]+l4[i,j]*w1l
			
	# Crop and normalize before next layer
			
	l4=l4l[K-9:N+K+9,K-9:N+K+9]/np.max(l4l)
	del l4l
	
	
	
	### Layer 2/3 ###
	
	K=9	# Spatial spread of L2/3 to L4 afferent weights
	
	l23=np.zeros((N,N), dtype=float)
	
	# Activity from afferent connections
	
	for i in range (0,N):
		for j in range (0,N):
		
			n_iso_domains=int(N/40)	# Size of iso-orientation domain = (40x40)
			
			# Select orientation selective afferent weight matrix
			
			i1,j1=int(i/40),int(j/40)
			if i1 in list(range(0,n_iso_domains+4,2)):
				if j1 in list(range(0,n_iso_domains+4,2)):
					w2=orr2[0,:,:]
				else:
					w2=orr2[3,:,:]
			else:
				if j1 in list(range(0,n_iso_domains+4,2)):
					w2=orr2[1,:,:]
				else:
					w2=orr2[2,:,:]
					
			# Generate L2/3 node activity
						
			l23[i,j]=np.sum(np.multiply(l4[i:i+2*K+1,j:j+2*K+1],w2))
			
	del l4
	
	# Activity from lateral connections
	
	K=75	# padding
	
	l23l=np.zeros((N+2*K,N+2*K),dtype=float)
		
	for i in range (0,N):
		for j in range (0,N):
		
			# Generate L2/3 node activity
		
			l23l[i:i+2*K+1,j:j+2*K+1]=l23l[i:i+2*K+1,j:j+2*K+1]+l23[i,j]*w2l
			
	# Crop and normalize before next layer
			
	l23=l23l[K-12:N+K+12,K-12:N+K+12]/np.max(l23l)
	del l23l
	
	
	
	### Layer 5 ###
	
	N=int(N/2)
	K=12	# Spatial spread of L5 to L2/3 afferent weights
	
	l5=np.zeros((N,N),dtype=float)
	
	# Activity from afferent connections
	
	for i in range (0,N):
		for j in range (0,N):
		
			n_iso_domains=int(N/20)	# Size of iso-orientation domain = (20x20)
			
			# Select orientation selective afferent weight matrix
		
			i1,j1=int(i/20),int(j/20)
			if i1 in list(range(0,n_iso_domains+4,2)):
				if j1 in list(range(0,n_iso_domains+4,2)):
					w3=orr3[0,:,:]
				else:
					w3=orr3[3,:,:]
			else:
				if j1 in list(range(0,n_iso_domains+4,2)):
					w3=orr3[1,:,:]
				else:
					w3=orr3[2,:,:]
					
			# Generate L5 node activity

			l5[i,j]=np.sum(np.multiply(l23[i*2:i*2+2*K+1,j*2:j*2+2*K+1],w3))
			
	del l23
	
	# Activity from lateral connections
	
	K=37	# padding
	
	l5l=np.zeros((N+2*K,N+2*K),dtype=float)
	
	for i in range (0,N):
		for j in range (0,N):
		
			# Generate L5 node activity
		
			l5l[i:i+2*K+1,j:j+2*K+1]=l5l[i:i+2*K+1,j:j+2*K+1]+l5[i,j]*w3l
			
	# Crop and normalize before next layer
			
	l5=l5l[K-12:N+K+12,K-12:N+K+12]/np.max(l5l)
	del l5l
	
	
	
	### Layer 2/3 (back transition) ###
	
	N=N*2
	K=12
	l23=np.zeros((N,N),dtype=float)
	
	# Activity from afferent connections
	
	for i in range (0,N):
		for j in range (0,N):
		
			n_iso_domains=int(N/40)	# Size of iso-orientation domain = (40x40)
			
			# Select orientation selective afferent weight matrix
			
			i1,j1=int(i/40),int(j/40)
			if i1 in list(range(0,n_iso_domains+4,2)):
				if j1 in list(range(0,n_iso_domains+4,2)):
					w3=orr3[0,:,:]
				else:
					w3=orr3[3,:,:]
			else:
				if j1 in list(range(0,n_iso_domains+4,2)):
					w3=orr3[1,:,:]
				else:
					w3=orr3[2,:,:]
					
			# Generate L2/3 node activity

			l23[i,j]=np.sum(np.multiply(l5[int(i*0.5):int(i*0.5)+2*K+1,int(j*0.5):int(j*0.5)+2*K+1],w3))
	
	# Activity from lateral connections
	
	K=75	# padding
	
	l23l=np.zeros((N+2*K,N+2*K),dtype=float)
		
	for i in range (0,N):
		for j in range (0,N):
		
			# Generate L2/3 node activity
		
			l23l[i:i+2*K+1,j:j+2*K+1]=l23l[i:i+2*K+1,j:j+2*K+1]+l23[i,j]*w2l
			
	# Crop and normalize for output
			
	l23=l23l[K:N+K,K:N+K]/np.max(l23l)
	del l23l
	
	
	return l23
	
	
