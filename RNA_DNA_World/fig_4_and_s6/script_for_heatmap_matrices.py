import numpy as np

prob=0.04
ert=5

basic_directory='./data/'
a=np.array([0.04,0.02])
a=np.concatenate((a,np.linspace(0.01, 0.002, num=5)))
a=np.concatenate((a,np.linspace(0.001, 0.0002, num=5)))

for i in range (0,12):
	a[i]=round(a[i],4)

Ert=np.linspace(0,8,9)

m1,m2,m3=np.zeros((12,9),dtype=float),np.zeros((12,9),dtype=float),np.zeros((12,9),dtype=float)
m5,m6=np.zeros((12,9),dtype=float),np.zeros((12,9),dtype=float)

i1=0
for prob in a:
	j1=0
	for j in Ert:
		ert=int(j)
		target_directory=basic_directory+'p='+str(prob)+'/ert='+str(ert)+'/'
	
		data1=np.loadtxt(target_directory+'data_1.txt')
		r=data1[:,1]
		rest=data1[:,2]
	
		m=np.divide(r,(r+rest))
		m1[i1,j1]=np.mean(m)
	
	
		data2=np.loadtxt(target_directory+'data_2.txt')
		rtdr=data2[:,1]
		rtdr=np.mean(rtdr)
		m2[i1,j1]=rtdr
	
	
		data3=np.loadtxt(target_directory+'data_3.txt')
		re9=data3[:,1]
		re=data3[:,2]
	
		m=np.divide(re9,(re9+re))
		m3[i1,j1]=np.mean(m)
		
		
		data4=np.loadtxt(target_directory+'data_4.txt')
		re567=data4[:,1]
		re8=data4[:,2]
		re=data4[:,3]
		
		m=np.divide(re567,re)
		m5[i1,j1]=np.mean(m)
		
		m=np.divide(re8,re)
		m6[i1,j1]=np.mean(m)
		
		j1+=1
		
	i1+=1
	
m4=np.ones((12,9),dtype=float)-m3-m5-m6

np.save('m1',m1)
np.save('m2',m2)
np.save('m3',m3)
np.save('m4',m4)
np.save('m5',m5)
np.save('m6',m6)
