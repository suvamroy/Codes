import numpy as np
import matplotlib.pyplot as plt
from math import exp

x=np.linspace(0,8,81)
y1,y2=[],[]

for i in range (0,len(x)):
	y1+=[1/(1+exp(8-x[i]))]
	y2+=[1/(1+exp(3-x[i]))]
	
plt.plot(x,y2)
plt.plot(x,y1)
plt.legend(['RNA to DNA and DNA to RNA','RNA to RNA'])
plt.xlabel('Error threshold')
plt.ylabel('Probability of accurate replication')
plt.savefig('Fig_S1.eps',format='eps',dpi=300)
