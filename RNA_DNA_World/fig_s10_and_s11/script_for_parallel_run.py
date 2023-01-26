import os
import numpy as np

prob=0.04
ert=5

basic_directory='./data/'
a=np.array([0.04,0.02])
a=np.concatenate((a,np.linspace(0.01, 0.002, num=5)))
a=np.concatenate((a,np.linspace(0.001, 0.0002, num=5)))

Ert=np.linspace(0,8,9)
for i in a:
	for j in Ert:
		prob=round(i,4)
		ert=int(j)
		target_directory=basic_directory+'p='+str(prob)+'/ert='+str(ert)+'/'
		cmd = 'python3 main.py '+str(prob)+' ' +str(ert) +' '+target_directory +' &'
		os.system(cmd)
		print(cmd)


	
