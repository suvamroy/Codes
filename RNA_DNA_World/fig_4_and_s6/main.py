###### Load modules ###################

import os
import sys
import numpy as np
import random
from math import exp

target_directory=sys.argv[-1]
if not os.path.isdir(target_directory):
	os.makedirs(target_directory)
	
prob=float(sys.argv[1])			
ert=float(sys.argv[2])	

###################################################################################

# Initial condition

def f(x):	# The function theta in the manuscript
	if x>0:
		return 1
	else:
		return 0

N,V=400,1000

kr_rr,kr_rd,kr_dr,kr_dd,kf_rr,kf_rd,kf_dr,kf_dd=0.00255,0.001646,0.000823,0.00255,0.46447,0.226,0.1434,0

hr,hd=0.0008,0.00008

cell=[[]]*N
for i in range (0,N):
	cell[i]=[random.uniform(prob-prob*0.5,prob+prob*0.5),[1]*10,[],[],[],[],[],[]]
	


# Time evolution

t,tot_rate=0.0,[0]*N
T=30000
while t<T:

	# Export the time evolution data for the avg abundance of
	# replicase and all other type of strands per droplet
	
	tr,p,r,td,tdp,tdr=0,0,0,0,0,0
	for i in range (0,N):
		tr+=len(cell[i][1])
		p+=len(cell[i][2])
		r+=len(cell[i][3])
		td+=len(cell[i][4])
		tdp+=len(cell[i][5])
		tdr+=len(cell[i][6])
	
	if t>=25000:
		tr,p,r,td,tdp,tdr=tr/N,p/N,r/N,td/N,tdp/N,tdr/N
		m=open(target_directory+'data_1.txt','a+')
		m.write('%f' %t + '\t' + '%f' %r + '\t' + '%f' %(tr+p+td+tdp+tdr) + '\n')
		m.close()
		
	# Export the time evolution data for the percentage 
	# of droplets containing both R and T_dR
	
	both=0
	for i in range (0,N):
		if len(cell[i][3])>0 and len(cell[i][6])>0:
			both+=1
	
	if t>=25000:
		m=open(target_directory+'data_2.txt','a+')
		m.write('%f' %t + '\t' + '%f' %(both/N) + '\n')
		m.close()
	
	
	
	# Reactions:	
	
	# Calculate total rates of each protocell and time step size:
	
	re1,re2,re3,re4,re5,re6,re7,re8,re9=0,0,0,0,0,0,0,0,0
	
	for i in range (0,N):
		tr,p,r,td,tdp,tdr=len(cell[i][1]),len(cell[i][2]),len(cell[i][3]),len(cell[i][4]),len(cell[i][5]),len(cell[i][6])
		
		a=[kr_rr*tr,kr_rd*tr,kr_dr*(td+tdp+tdr),kr_dd*(td+tdp+tdr),kf_rr*f(tr)*r,kf_rd*f(tr)*r,kf_dr*f(td)*r,kf_dr*f(tdp)*r,kf_dr*f(tdr)*r,hr*(tr+p+r),hd*(td+tdp+tdr)]
		
		re1+=a[0]
		re2+=a[1]
		re3+=a[2]
		re4+=a[3]
		re5+=a[4]
		re6+=a[5]
		re7+=a[6]
		re8+=a[7]
		re9+=a[8]
		
		if tr+p+r+td+tdp+tdr>=V:
			a=[0,0,0,0,0,0,0,0,0,hr*(tr+p+r),hd*(td+tdp+tdr)]
			
		tot_rate[i]=sum(a)
		
		
	if t>=25000:
	
		# Export the time evolution data for the propensities 
		# of reaction-9 and all other polymerization reactiions
	
		re1,re2,re3,re4,re5,re6,re7,re8,re9=re1/N,re2/N,re3/N,re4/N,re5/N,re6/N,re7/N,re8/N,re9/N
		m=open(target_directory+'data_3.txt','a+')
		m.write('%f' %t + '\t' + '%f' %re9 + '\t' + '%f' %(re1+re2+re3+re4+re5+re6+re7+re8) + '\n')
		m.close()
		
		# Export the time evolution data for the propensities 
		# of reaction-5,6,7,8 and all polymerization reactiions
		
		m=open(target_directory+'data_4.txt','a+')
		m.write('%f' %(re5+re6+re7) + '\t' + '%f' %re8 + '\t' + '%f' %(re1+re2+re3+re4+re5+re6+re7+re8+re9) + '\n')
		m.close()
		
		
	# Calculate time step size
		
	max_rate=max(tot_rate)
	if max_rate<=0.0:
		break
	dt=1/max_rate
	if t+dt>T:
		break
	t+=dt
	
	
	
	# Reaction choosing inside each protocell:
	
	for i in range (0,N):
		if tot_rate[i]>0.0 and random.uniform(0,1)<tot_rate[i]/max_rate:
			
			tr,p,r,td,tdp,tdr=len(cell[i][1]),len(cell[i][2]),len(cell[i][3]),len(cell[i][4]),len(cell[i][5]),len(cell[i][6])
			
			if tr+p+r+td+tdp+tdr<V:
				a=[kr_rr*tr,kr_rd*tr,kr_dr*(td+tdp+tdr),kr_dd*(td+tdp+tdr),kf_rr*f(tr)*r,kf_rd*f(tr)*r,kf_dr*f(td)*r,kf_dr*f(tdp)*r,kf_dr*f(tdr)*r,hr*(tr+p+r),hd*(td+tdp+tdr)]
		
			else:
				a=[0,0,0,0,0,0,0,0,0,hr*(tr+p+r),hd*(td+tdp+tdr)]
		
		
			reaction=np.random.choice(range(1,12), 1, p=np.array(a)/tot_rate[i])[0]
			
			
			if reaction==1:
				if random.uniform(0,1)<cell[i][0]:
					cell[i][3]+=[1]
				else:
					cell[i][2]+=[1]
			
			if reaction==2:
				if random.uniform(0,1)<cell[i][0]:
					cell[i][6]+=[1]
				else:
					cell[i][5]+=[1]
				
			if reaction==3:
				if random.uniform(0,1)<cell[i][0]:
					cell[i][3]+=[1]
				else:
					cell[i][2]+=[1]
		
			if reaction==4:
				if random.uniform(0,1)<cell[i][0]:
					cell[i][6]+=[1]
				else:
					cell[i][5]+=[1]
		
			if reaction==5:
				if random.uniform(0,1)<1/(1+exp(8-ert)):
					cell[i][1]+=[1]
				else:
					cell[i][2]+=[1]
		
			if reaction==6:
				if random.uniform(0,1)<1/(1+exp(3-ert)):
					cell[i][4]+=[1]
				else:
					cell[i][5]+=[1]
		
			if reaction==7:
				if random.uniform(0,1)<1/(1+exp(3-ert)):
					cell[i][1]+=[1]
				else:
					cell[i][2]+=[1]
		
			if reaction==8:
				cell[i][2]+=[1]
				
			if reaction==9:
				if random.uniform(0,1)<1/(1+exp(3-ert)):
					cell[i][3]+=[1]
				else:
					cell[i][2]+=[1]
			
			if reaction==10:
				mol=np.random.choice(['tr','p','r'], 1, p=np.array([tr,p,r])/(tr+p+r))[0]
				if mol=='tr':
					del cell[i][1][-1]
				if mol=='p':
					del cell[i][2][-1]
				if mol=='r':
					del cell[i][3][-1]
			
			if reaction==11:
				mol=np.random.choice(['td','tdp','tdr'], 1, p=np.array([td,tdp,tdr])/(td+tdp+tdr))[0]
				if mol=='td':
					del cell[i][4][-1]
				if mol=='tdp':
					del cell[i][5][-1]
				if mol=='tdr':
					del cell[i][6][-1]
				
	
	
	
	
	# Protocell division:				
			
	vol=[]
	for i in range (0,N):
		tr,p,r,td,tdp,tdr=len(cell[i][1]),len(cell[i][2]),len(cell[i][3]),len(cell[i][4]),len(cell[i][5]),len(cell[i][6])					
		if tr+p+r+td+tdp+tdr>=V/5:
			vol+=[1]
		else:
			vol+=[0]
			
	# If all droplets contain >=V/5 number of 
	# strands then division round will happen
	
	if sum(vol)==N:
		index=list(range(0,N))
		np.random.shuffle(index)
		i1=0	
		while i1<N:
		
			# assign fitness and select winner droplet
		
			i,j=index[i1],index[i1+1]
			
			tr,p,r,td,tdp,tdr=len(cell[i][1]),len(cell[i][2]),len(cell[i][3]),len(cell[i][4]),len(cell[i][5]),len(cell[i][6])
			tr1,p1,r1,td1,tdp1,tdr1=len(cell[j][1]),len(cell[j][2]),len(cell[j][3]),len(cell[j][4]),len(cell[j][5]),len(cell[j][6])					
			
			f1,f2=r/(tr+p+td+tdp+tdr),r1/(tr1+p1+td1+tdp1+tdr1)
			
			if f1==0 and f2==0:
				f1,f2=(tr+p+td+tdp+tdr),(tr1+p1+td1+tdp1+tdr1)

			if random.uniform(0,1)<=f2/(f1+f2):
				i,j=j,i
				tr,p,r,td,tdp,tdr=len(cell[i][1]),len(cell[i][2]),len(cell[i][3]),len(cell[i][4]),len(cell[i][5]),len(cell[i][6])
		
		
			# Eliminate loser droplet

			cell[j]=[cell[i][0],[],[],[],[],[],[]]
			
			# Random distribution of strands among two daughter protocells:
			
			for k in range (0,tr):
				if random.uniform(0,1)<0.5:
					cell[j][1]+=[cell[i][1][k]]
					cell[i][1][k]='X'
			cell[i][1]=[k for k in cell[i][1] if k!='X']
			
			for k in range (0,p):
				if random.uniform(0,1)<0.5:
					cell[j][2]+=[cell[i][2][k]]
					cell[i][2][k]='X'
			cell[i][2]=[k for k in cell[i][2] if k!='X']
			
			for k in range (0,r):
				if random.uniform(0,1)<0.5:
					cell[j][3]+=[cell[i][3][k]]
					cell[i][3][k]='X'
			cell[i][3]=[k for k in cell[i][3] if k!='X']
			
			for k in range (0,td):
				if random.uniform(0,1)<0.5:
					cell[j][4]+=[cell[i][4][k]]
					cell[i][4][k]='X'
			cell[i][4]=[k for k in cell[i][4] if k!='X']
			
			for k in range (0,tdp):
				if random.uniform(0,1)<0.5:
					cell[j][5]+=[cell[i][5][k]]
					cell[i][5][k]='X'
			cell[i][5]=[k for k in cell[i][5] if k!='X']
			
			for k in range (0,tdr):
				if random.uniform(0,1)<0.5:
					cell[j][6]+=[cell[i][6][k]]
					cell[i][6][k]='X'
			cell[i][6]=[k for k in cell[i][6] if k!='X']
			
			
			i1+=+2
			



# End	
		
					

	
			
