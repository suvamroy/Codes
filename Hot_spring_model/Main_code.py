#.............................................................................................

# This code can generate the data files corresponding to Figure 2, 3 and 6 of the main text. #

# ...........................................................................................#

# Import python modules

import numpy as np
import random
from math import exp


######### Initial condition setup ###################	

# k_r0 and k_r correspond to the non-enzymatic replication rates (excluding the error part) for
# circ ssRNA to circ dsRNA and circ dsRNA to open-ended ssRNA creations (the 1st process requires 
# extending a 8-mer primer to 200 length, therefore effective length is 192 instead of 200)


N,k_r,k_r0,k_fast,h=30,exp(-3.22-0.005*200),exp(-3.22-0.005*192),0.362,0.0008
P_r,P_c,P_n,P_p,P_kill,P_hop=0.03,0.03,0.03,0.03,0.01,1
b=1.0

V=np.zeros((N,N), dtype=float)
S,S_max=np.zeros((N,N), dtype=float),np.zeros((N,N), dtype=float)
cell=np.empty((N,N), dtype=list)

for i in range (0,N):
	for j in range (0,N):
	
		# 1 circ-ssRNA per site initially (figure 2, 3)
		
		a=list(np.random.normal(0.35,0.0667, size=(1,1)).flatten())
		a=[k_r0*exp(-2.8*k) for k in a]
		cell[i][j]=[a,[],[],[],[],[],[]]
		
		# 5 circ-ssRNA at a central site initially (figure 6)
		'''
		if i==14 and j==14:
			a=list(np.random.normal(0.35,0.0667, size=(1,5)).flatten())
			a=[k_r0*exp(-2.8*k) for k in a]
			cell[i][j]=[a,[],[],[],[],[],[]]
		else:
			cell[i][j]=[[],[],[],[],[],[],[]]'''
			
			
		V[i][j]=100
		S[i][j]=80
		S_max[i][j]=80
		
# Matrices to track different types of ribozymes and total number of strands of each site

ribozymes=np.zeros((N,N), dtype=int)
replicase=np.zeros((N,N), dtype=int)
cyclase=np.zeros((N,N), dtype=int)
nucsy=np.zeros((N,N), dtype=int)
peptr=np.zeros((N,N), dtype=int)
numbers=np.zeros((N,N), dtype=int)



##################### Time evolution #########################

tot_rate=np.zeros((N,N), dtype=float)
T=8	# duration of each phase 8 hours

for day in range (0,1000):

	####### Dry phase #########

	##### Export data #########
	
	no=0
	for i in range (0,N):
		for i1 in range (0,N): 
			s,d,l,r,c,n,p=len(cell[i][i1][0]),len(cell[i][i1][1]),len(cell[i][i1][2]),len(cell[i][i1][3]),len(cell[i][i1][4]),len(cell[i][i1][5]),len(cell[i][i1][6])
			no+=n
			
			numbers[i][i1]=s+d+l+r+c+n+p
			replicase[i][i1]=r
			cyclase[i][i1]=c
			nucsy[i][i1]=n
			peptr[i][i1]=p
			
			# Make number and ribzyme matrices binary
			
			if s+d+l+r+c+n+p>0:
				numbers[i][i1]=1
			if r>0:
				replicase[i][i1]=1
			if c>0:
				cyclase[i][i1]=1
			if n>0:
				nucsy[i][i1]=1
			if p>0:
				peptr[i][i1]=1
				
	ribozymes=numbers+replicase+cyclase+nucsy+peptr	
	
	# Export data files for figure 2
	'''
	fl=open('D_1.28.txt','a+')		
	fl.write('%d' %day + '\t' + '%d' %(list(numbers.flatten()).count(5)) + '\n')
	fl.close()'''
	
	# Export data files for figure 3
	
	fl=open('D_1.28_P_kill_0.01.txt','a+')		
	fl.write('%d' %day + '\t' + '%d' %(list(numbers.flatten()).count(0)) + '\n')
	fl.close()
	
	# Export data files for figure 6
	'''
	if day==0:
		np.savetxt('Day_1.txt',ribozymes)
	if day==299:
		np.savetxt('Day_300.txt',ribozymes)
	if list(ribozymes.flatten()).count(5)/900>=0.9:
		np.savetxt('Day_90_percent.txt',ribozymes)
		break'''
		
		
		
		

	# Reactions:	
	
	# Calculate total rates of each cell and time step size:
	
	t=0.0
	while t<T:
	
		for i in range (0,N):
			for i1 in range (0,N):
				f=(S[i][i1]+b*len(cell[i][i1][5]))/S_max[i][i1]	# Rate reduction factor for finite monomer
				
				if sum([len(cell[i][i1][0]),len(cell[i][i1][1]),len(cell[i][i1][2]),len(cell[i][i1][3]),len(cell[i][i1][4]),len(cell[i][i1][5]),len(cell[i][i1][6])])<1000:
		
					a=[sum(cell[i][i1][0])*f,sum(cell[i][i1][1])*f,k_fast*len(cell[i][i1][3])*len(cell[i][i1][0])*f/100,k_fast*len(cell[i][i1][3])*len(cell[i][i1][1])*f/100,k_fast*len(cell[i][i1][4])*len(cell[i][i1][2])/100,h*(len(cell[i][i1][0])+len(cell[i][i1][1])+len(cell[i][i1][2])+len(cell[i][i1][3])+len(cell[i][i1][4])+len(cell[i][i1][5])+len(cell[i][i1][6]))]
				
				else:
					a=[sum(cell[i][i1][0])*f,0.0,k_fast*len(cell[i][i1][3])*len(cell[i][i1][0])*f/100,0.0,k_fast*len(cell[i][i1][4])*len(cell[i][i1][2])/100,h*(len(cell[i][i1][0])+len(cell[i][i1][1])+len(cell[i][i1][2])+len(cell[i][i1][3])+len(cell[i][i1][4])+len(cell[i][i1][5])+len(cell[i][i1][6]))]
					
				tot_rate[i][i1]=sum(a)
		
		max_rate=np.max(tot_rate)
		if max_rate<=0.0:
			break
		dt=1/max_rate
		if dt>=T:
			dt=T
		if t+dt<=T:
			t+=dt
		else:
			break
	
	
	
	# Reaction choosing inside each cell:
		
		g=0
		for i in range (0,N):
			for i1 in range (0,N):
				if tot_rate[i][i1]>0.0 and random.uniform(0,1)<tot_rate[i][i1]*dt:
		
					f=(S[i][i1]+b*len(cell[i][i1][5]))/S_max[i][i1]	# Rate reduction factor for finite monomer
					
					if sum([len(cell[i][i1][0]),len(cell[i][i1][1]),len(cell[i][i1][2]),len(cell[i][i1][3]),len(cell[i][i1][4]),len(cell[i][i1][5]),len(cell[i][i1][6])])<1000:
			
						a=[sum(cell[i][i1][0])*f,sum(cell[i][i1][1])*f,k_fast*len(cell[i][i1][3])*len(cell[i][i1][0])*f/100,k_fast*len(cell[i][i1][3])*len(cell[i][i1][1])*f/100,k_fast*len(cell[i][i1][4])*len(cell[i][i1][2])/100,h*(len(cell[i][i1][0])+len(cell[i][i1][1])+len(cell[i][i1][2])+len(cell[i][i1][3])+len(cell[i][i1][4])+len(cell[i][i1][5])+len(cell[i][i1][6]))]
						
					else:
						a=[sum(cell[i][i1][0])*f,0.0,k_fast*len(cell[i][i1][3])*len(cell[i][i1][0])*f/100,0.0,k_fast*len(cell[i][i1][4])*len(cell[i][i1][2])/100,h*(len(cell[i][i1][0])+len(cell[i][i1][1])+len(cell[i][i1][2])+len(cell[i][i1][3])+len(cell[i][i1][4])+len(cell[i][i1][5])+len(cell[i][i1][6]))]
			
					reaction=np.random.choice(['cs_cd','cd_oe','r_cs_cd','r_cd_oe','c_oe_cs','death'], 1, p=np.array(a)/tot_rate[i][i1])[0]
			
			
			# Circ ssRNA to Circ dsRNA:
			
					if reaction=='cs_cd':
						j=np.random.choice(range(0,len(cell[i][i1][0])), 1, p=np.array(cell[i][i1][0])/sum(cell[i][i1][0]))[0]
						cell[i][i1][1]+=[(cell[i][i1][0][j]*k_r)/k_r0]
						del cell[i][i1][0][j]
						g+=1
						
				
			# Circ dsRNA to ribozymes or non-catalytic open-ended RNA:
			
					if reaction=='cd_oe':
						p=random.uniform(0,1)
						if p<P_r:
							cell[i][i1][3]+=[k_fast]
						elif p>=P_r and p<P_r+P_c:
							cell[i][i1][4]+=[k_fast]
						elif p>=P_r+P_c and p<P_r+P_c+P_n:
							cell[i][i1][5]+=[0]
						elif p>=P_r+P_c+P_n and p<P_r+P_c+P_n+P_p:
							cell[i][i1][6]+=[0]
						else:
							cell[i][i1][2]+=[0]
						g+=1
						
				
			# Replicase catalyzed Circ ssRNA to Circ dsRNA:
					
					if reaction=='r_cs_cd':
						j=np.random.choice(range(0,len(cell[i][i1][0])), 1)[0]
						cell[i][i1][1]+=[(cell[i][i1][0][j]*k_r)/k_r0]
						del cell[i][i1][0][j]
						g+=1
				
				
			# Replicase catalyze Circ dsRNA to ribozymes or non-catalytic open-ended RNA:
				
					if reaction=='r_cd_oe':
						p=random.uniform(0,1)
						if p<P_r:
							cell[i][i1][3]+=[k_fast]
						elif p>=P_r and p<P_r+P_c:
							cell[i][i1][4]+=[k_fast]
						elif p>=P_r+P_c and p<P_r+P_c+P_n:
							cell[i][i1][5]+=[0]
						elif p>=P_r+P_c+P_n and p<P_r+P_c+P_n+P_p:
							cell[i][i1][6]+=[0]
						else:
							cell[i][i1][2]+=[0]
						g+=1
				
				
			# Cyclase catalyzed non-catalytic open-ended RNA to Circ ssRNA
					
					if reaction=='c_oe_cs':
						cell[i][i1][0]+=[k_r0*exp(-2.8*np.random.normal(0.35,0.0667))]
						del cell[i][i1][2][-1]
				
				
			# Degradation of strands:
						
					if reaction=='death':
						s,d,l,r,c,n,p=len(cell[i][i1][0]),len(cell[i][i1][1]),len(cell[i][i1][2]),len(cell[i][i1][3]),len(cell[i][i1][4]),len(cell[i][i1][5]),len(cell[i][i1][6])
						dead=np.random.choice(['s','d','l','r','c','n','p'], 1, p=np.array([s,d,l,r,c,n,p])/(s+d+l+r+c+n+p))[0]
				
						if dead=='s':
							j=np.random.choice(range(0,s), 1)[0]
							del cell[i][i1][0][j]
				
						if dead=='d':
							j=np.random.choice(range(0,d), 1)[0]
							del cell[i][i1][1][j]
					
						if dead=='l':
							j=np.random.choice(range(0,l), 1)[0]
							del cell[i][i1][2][j]
					
						if dead=='r':
							j=np.random.choice(range(0,r), 1)[0]
							del cell[i][i1][3][j]
					
						if dead=='c':
							j=np.random.choice(range(0,c), 1)[0]
							del cell[i][i1][4][j]
					
						if dead=='n':
							j=np.random.choice(range(0,n), 1)[0]
							del cell[i][i1][5][j]
					
						if dead=='p':
							j=np.random.choice(range(0,p), 1)[0]
							del cell[i][i1][6][j]
					
					
			# Threshold volume increase in presence of peptidyl transferase:
				
					if len(cell[i][i1][6])>0:
						V[i][i1]=100+20*len(cell[i][i1][6])
						
				
	# Change in S values after reactions:
				
		S_tot=np.sum(S)+(b*no-g)
		for i in range (0,N):
			for i1 in range (0,N):
				S[i][i1]=(S_tot/N**2)
				if S[i][i1]<0:
					S[i][i1]=0
				if S[i][i1]>S_max[i][i1]:
					S[i][i1]=S_max[i][i1]
	
	
		
	
	
	######## Wet phase ##############
	
	#### Protocell relocation #######
	'''
	for i in range (0,N**2-1):
		j=random.randint(i,N**2-1)
		i1,i2=int(i/N),(i%N)
		j1,j2=int(j/N),(j%N)
		cell[i1][i2],cell[j1][j2]=cell[j1][j2],cell[i1][i2]'''
		
	
	##### Protocell degration #####
	'''
	for i in range (0,N):
		for j in range (0,N):
			if random.uniform(0,1)<P_kill:
				cell[i][j]=[[],[],[],[],[],[],[]]'''
	
	
	##### Protocell division #####
	
	used=np.zeros((N,N), dtype=int)
	index=[]
	for i in range (0,N):
		for j in range (0,N):
			index+=[str(i)+' '+str(j)]
			
	np.random.shuffle(index)
	for k1 in range (0,N**2):	
		i,i1=int(index[k1].split()[0]),int(index[k1].split()[1])		
		s,d,l,r,c,n,p=len(cell[i][i1][0]),len(cell[i][i1][1]),len(cell[i][i1][2]),len(cell[i][i1][3]),len(cell[i][i1][4]),len(cell[i][i1][5]),len(cell[i][i1][6])	
		if s+d+l+r+c+n+p>=V[i][i1]:

			# Find neighbors

			# Apply periodic boundary condition
					
			if i==0:
				if i1==0:
					nei=[[N-1,N-1],[N-1,i1],[N-1,i1+1],[i,N-1],[i,i1+1],[i+1,N-1],[i+1,i1],[i+1,i1+1]]
				elif i1==N-1:
					nei=[[N-1,i1-1],[N-1,i1],[N-1,0],[i,i1-1],[i,0],[i+1,i1-1],[i+1,i1],[i+1,0]]
				else:
					nei=[[N-1,i1-1],[N-1,i1],[N-1,i1+1],[i,i1-1],[i,i1+1],[i+1,i1-1],[i+1,i1],[i+1,i1+1]]

			elif i==N-1:
				if i1==0:
					nei=[[i-1,N-1],[i-1,i1],[i-1,i1+1],[i,N-1],[i,i1+1],[0,N-1],[0,i1],[0,i1+1]]
				elif i1==N-1:
					nei=[[i-1,i1-1],[i-1,i1],[i-1,0],[i,i1-1],[i,0],[0,i1-1],[0,i1],[0,0]]
				else:
					nei=[[i-1,i1-1],[i-1,i1],[i-1,i1+1],[i,i1-1],[i,i1+1],[0,i1-1],[0,i1],[0,i1+1]]
			else:
				if i1==0:
					nei=[[i-1,N-1],[i-1,i1],[i-1,i1+1],[i,N-1],[i,i1+1],[i+1,N-1],[i+1,i1],[i+1,i1+1]]
				elif i1==N-1:
					nei=[[i-1,i1-1],[i-1,i1],[i-1,0],[i,i1-1],[i,0],[i+1,i1-1],[i+1,i1],[i+1,0]]
				else:
					nei=[[i-1,i1-1],[i-1,i1],[i-1,i1+1],[i,i1-1],[i,i1+1],[i+1,i1-1],[i+1,i1],[i+1,i1+1]]
					
			
			# If a neighbor already has been replaced by another dividing protocell, don't consider it
			
			nei_copy=nei	
			for k in range (0,8):
				if used[nei[k][0]][nei[k][1]]==1:
					nei[k]='X'
			nei=[k for k in nei if k!='X']
			
			# Assign weakness values to neighbors
				
			weakness=[]
			for k in range(0,len(nei)):
				j,j1=nei[k][0],nei[k][1]
				s1,d1,l1,r1,c1,n1,p1=len(cell[j][j1][0]),len(cell[j][j1][1]),len(cell[j][j1][2]),len(cell[j][j1][3]),len(cell[j][j1][4]),len(cell[j][j1][5]),len(cell[j][j1][6])
					
				if (s+d+l+r+c+n+p)>(s1+d1+l1+r1+c1+n1+p1):
					weakness+=[(s+d+l+r+c+n+p)-(s1+d1+l1+r1+c1+n1+p1)]
				else:
					weakness+=[0]
					
					
			# Choose neighbor to eliminate based on weakness values
			
			if sum(weakness)>0:			
				k=np.random.choice(list(range(0,len(nei))), 1, p=np.array(weakness)/sum(weakness))[0]
				j,j1=nei[k][0],nei[k][1]	
			else:	
				if len(nei)>0:
					k=np.random.choice(list(range(0,len(nei))), 1)[0]
					j,j1=nei[k][0],nei[k][1]
				else:
					k=np.random.choice(list(range(0,len(nei_copy))), 1)[0]
					j,j1=nei_copy[k][0],nei_copy[k][1]
					
			
			cell[j][j1]=[[],[],[],[],[],[],[]]
			used[i][i1],used[j][j1]=1,1
				
			
			# Random distribution of strands among two daughter cells:
			
			for k in range (0,s):
				if random.uniform(0,1)<0.5:
					cell[j][j1][0]+=[cell[i][i1][0][k]]
					cell[i][i1][0][k]='X'
			cell[i][i1][0]=[k for k in cell[i][i1][0] if k!='X']
			
			for k in range (0,d):
				if random.uniform(0,1)<0.5:
					cell[j][j1][1]+=[cell[i][i1][1][k]]
					cell[i][i1][1][k]='X'
			cell[i][i1][1]=[k for k in cell[i][i1][1] if k!='X']
			
			for k in range (0,l):
				if random.uniform(0,1)<0.5:
					cell[j][j1][2]+=[cell[i][i1][2][k]]
					cell[i][i1][2][k]='X'
			cell[i][i1][2]=[k for k in cell[i][i1][2] if k!='X']
			
			for k in range (0,r):
				if random.uniform(0,1)<0.5:
					cell[j][j1][3]+=[cell[i][i1][3][k]]
					cell[i][i1][3][k]='X'
			cell[i][i1][3]=[k for k in cell[i][i1][3] if k!='X']
			
			for k in range (0,c):
				if random.uniform(0,1)<0.5:
					cell[j][j1][4]+=[cell[i][i1][4][k]]
					cell[i][i1][4][k]='X'
			cell[i][i1][4]=[k for k in cell[i][i1][4] if k!='X']
			
			for k in range (0,n):
				if random.uniform(0,1)<0.5:
					cell[j][j1][5]+=[cell[i][i1][5][k]]
					cell[i][i1][5][k]='X'
			cell[i][i1][5]=[k for k in cell[i][i1][5] if k!='X']
			
			for k in range (0,p):
				if random.uniform(0,1)<0.5:
					cell[j][j1][6]+=[cell[i][i1][6][k]]
					cell[i][i1][6][k]='X'
			cell[i][i1][6]=[k for k in cell[i][i1][6] if k!='X']
			
			if j>=1 and j<=N and j1>=1 and j1<=N:
				V[j-1][j1-1]=100+20*len(cell[j][j1][6])
				V[i-1][i1-1]=100+20*len(cell[i][i1][6])
				
	
	
	
	
	######### Gel phase ##########
	
	###### Strand diffusion ######
	
	# Change number of steps and P_hop to vary diffusion coefficient
	
	for step in range (0,1000):
		index=[]
		for i in range (0,N):
			for j in range (0,N):
				index+=[str(i)+' '+str(j)]
			
		np.random.shuffle(index)
		for k in range (0,N**2):
			i,i1=int(index[k].split()[0]),int(index[k].split()[1])
			
			# Find neighbors
			
			# Apply periodic boundary
			
			if i==0:
				if i1==0:
					nei=[[N-1,N-1],[N-1,i1],[N-1,i1+1],[i,N-1],[i,i1+1],[i+1,N-1],[i+1,i1],[i+1,i1+1]]
				elif i1==N-1:
					nei=[[N-1,i1-1],[N-1,i1],[N-1,0],[i,i1-1],[i,0],[i+1,i1-1],[i+1,i1],[i+1,0]]
				else:
					nei=[[N-1,i1-1],[N-1,i1],[N-1,i1+1],[i,i1-1],[i,i1+1],[i+1,i1-1],[i+1,i1],[i+1,i1+1]]

			elif i==N-1:
				if i1==0:
					nei=[[i-1,N-1],[i-1,i1],[i-1,i1+1],[i,N-1],[i,i1+1],[0,N-1],[0,i1],[0,i1+1]]
				elif i1==N-1:
					nei=[[i-1,i1-1],[i-1,i1],[i-1,0],[i,i1-1],[i,0],[0,i1-1],[0,i1],[0,0]]
				else:
					nei=[[i-1,i1-1],[i-1,i1],[i-1,i1+1],[i,i1-1],[i,i1+1],[0,i1-1],[0,i1],[0,i1+1]]
			else:
				if i1==0:
					nei=[[i-1,N-1],[i-1,i1],[i-1,i1+1],[i,N-1],[i,i1+1],[i+1,N-1],[i+1,i1],[i+1,i1+1]]
				elif i1==N-1:
					nei=[[i-1,i1-1],[i-1,i1],[i-1,0],[i,i1-1],[i,0],[i+1,i1-1],[i+1,i1],[i+1,0]]
				else:
					nei=[[i-1,i1-1],[i-1,i1],[i-1,i1+1],[i,i1-1],[i,i1+1],[i+1,i1-1],[i+1,i1],[i+1,i1+1]]
			
			
			# Total number of molecules in site
				
			s,d,l,r,c,n,p=len(cell[i][i1][0]),len(cell[i][i1][1]),len(cell[i][i1][2]),len(cell[i][i1][3]),len(cell[i][i1][4]),len(cell[i][i1][5]),len(cell[i][i1][6])
		
			mols=s+d+l+r+c+n+p
		
			for m in range (0,mols):
			
				# Choose random neighbor
				
				j=np.random.choice(list(range(0,8)), 1)[0]
				j,j1=nei[j][0],nei[j][1]
				s1,d1,l1,r1,c1,n1,p1=len(cell[j][j1][0]),len(cell[j][j1][1]),len(cell[j][j1][2]),len(cell[j][j1][3]),len(cell[j][j1][4]),len(cell[j][j1][5]),len(cell[j][j1][6])
				
				# If neighbor contains fewer strands then move strand to it with probability P_hop
					
				if (s+d+l+r+c+n+p)>(s1+d1+l1+r1+c1+n1+p1) and random.uniform(0,1)<P_hop:
					mol=np.random.choice(['d','s','l','r','c','n','p'], 1, p=np.array([d,s,l,r,c,n,p])/(d+(s+l+r+c+n+p)))[0]	
			
					if mol=='s':
						k1=np.random.choice(list(range(0,s)), 1)[0]
						cell[j][j1][0]+=[cell[i][i1][0][k1]]
						s-=1
						del cell[i][i1][0][k1]
					
					if mol=='d':
						k1=np.random.choice(list(range(0,d)), 1)[0]
						cell[j][j1][1]+=[cell[i][i1][1][k1]]
						d-=1
						del cell[i][i1][1][k1]
					
					if mol=='l':
						cell[j][j1][2]+=[cell[i][i1][2][-1]]
						l-=1
						del cell[i][i1][2][-1]
					
					if mol=='r':
						cell[j][j1][3]+=[cell[i][i1][3][-1]]
						r-=1
						del cell[i][i1][3][-1]
					
					if mol=='c':
						cell[j][j1][4]+=[cell[i][i1][4][-1]]
						c-=1
						del cell[i][i1][4][-1]
					
					if mol=='n':
						cell[j][j1][5]+=[cell[i][i1][5][-1]]
						n-=1
						del cell[i][i1][5][-1]
					
					if mol=='p':
						cell[j][j1][6]+=[cell[i][i1][6][-1]]
						p-=1
						del cell[i][i1][6][-1]
		
					
				

######### End #####################
