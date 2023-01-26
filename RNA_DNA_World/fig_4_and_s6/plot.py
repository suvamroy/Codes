import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

prob=0.04
ert=5

basic_directory='./data/'
a=np.array([0.04,0.02])
a=np.concatenate((a,np.linspace(0.01, 0.002, num=5)))
a=np.concatenate((a,np.linspace(0.001, 0.0002, num=5)))

for i in range (0,12):
	a[i]=round(a[i],4)


m1=np.load('m1.npy')
m2=np.load('m2.npy')
m3=np.load('m3.npy')
m4=np.load('m4.npy')
m5=np.load('m5.npy')
m6=np.load('m6.npy')


fig=plt.figure(figsize=(12,3))

plt.subplot(1,3,1)
sns.heatmap(m1,cmap='RdPu_r',yticklabels=a)
plt.ylabel('Non-enzymatic probability')
plt.xlabel('Enzymatic error threshold')
plt.title('A')

plt.subplot(1,3,2)
sns.heatmap(m2,cmap='RdPu_r',yticklabels=a)
plt.ylabel('Non-enzymatic probability')
plt.xlabel('Enzymatic error threshold')
plt.title('B')

plt.subplot(1,3,3)
sns.heatmap(m3,cmap='RdPu_r',yticklabels=a)
plt.ylabel('Non-enzymatic probability')
plt.xlabel('Enzymatic error threshold')
plt.title('C')

plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=0.6, hspace=None)
plt.savefig('Fig_4.eps',format='eps',dpi=300)
plt.show()

plt.clf()

fig=plt.figure(figsize=(12,3))

plt.subplot(1,3,1)
sns.heatmap(m4,cmap='RdPu_r',yticklabels=a)
plt.ylabel('Non-enzymatic probability')
plt.xlabel('Enzymatic error threshold')
plt.title('A')

plt.subplot(1,3,2)
sns.heatmap(m5,cmap='RdPu_r',yticklabels=a)
plt.ylabel('Non-enzymatic probability')
plt.xlabel('Enzymatic error threshold')
plt.title('B')

plt.subplot(1,3,3)
sns.heatmap(m6,cmap='RdPu_r',yticklabels=a)
plt.ylabel('Non-enzymatic probability')
plt.xlabel('Enzymatic error threshold')
plt.title('C')



plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=0.6, hspace=None)
plt.savefig('Fig_S6.eps',format='eps',dpi=300)
plt.show()

