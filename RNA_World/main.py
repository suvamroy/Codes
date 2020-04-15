################################# RNA WORLD #####################################

import os
import sys
import numpy as np
import random
import RNA
from itertools import groupby
from math import exp

############################ Parallel computing #################################

target_directory=sys.argv[-1] # last argument is the directory
rates_directory=sys.argv[-1] # before last argument is the directory with rates

if not os.path.isdir(target_directory):
	os.makedirs(target_directory)
file_param = open(target_directory+'parameters.dat','w')
list_param=['k_hyd','k_con','dry','semi_wet','N','N_max','p_uu','p_up','p_pp','Day_max','Length_of_the_day']
for name_param,param in zip(list_param,sys.argv[1:]):
	file_param.write(name_param+': '+param+'\n')
file_param.close()

################################ Parameters #####################################

k_hyd=float(sys.argv[1])					# Rate of hydrolytic cleavage				
k_con=float(sys.argv[2])					# Rate of concatenation				
dry=float(sys.argv[3])						# Duration of dry phase
semi_wet=float(sys.argv[4])					# Duration of semi-wet phase
N=int(sys.argv[5])						# Number of initial strands
N_max=int(sys.argv[6])						# Maximum number of strands
p_uu=float(sys.argv[7])						# Hydrolysis Prefactor when both neighbor nucleotides are paired 
p_up=float(sys.argv[8])						# Hydrolysis Prefactor when either one neighbor nucleotide is paired
p_pp=float(sys.argv[9])						# Hydrolysis Prefactor when both neighbor nucleotides are unpaired
Day_max=int(sys.argv[10])					# Total number of days
Length_of_the_day = float(sys.argv[11])

# ~ k_hyd=0.04						# Rate of hydrolytic cleavage				
# ~ k_con=0.62						# Rate of concatenation				
# ~ dry=12.0						# Duration of dry phase
# ~ semi_wet=8.0					# Duration of semi-wet phase
# ~ N=100						# Number of initial strands
# ~ N_max=1000						# Maximum number of strands
# ~ p_uu=1.0						# Hydrolysis Prefactor when both neighbor nucleotides are paired 
# ~ p_up=0.1						# Hydrolysis Prefactor when either one neighbor nucleotide is paired
# ~ p_pp=0.01						# Hydrolysis Prefactor when both neighbor nucleotides are unpaired
# ~ Day_max=200						# Total number of days
# ~ Length_of_the_day = 24.0

#################################################################################

wet=Length_of_the_day-(dry+semi_wet)			# Duration of wet phase
monomers=['A','U','G','C']				# Types of monomers (free nucleotides)
complement = {'A': 'U', 'C': 'G', 'U': 'A', 'G': 'C'}	# Complements of nucleotides

# Loading primer extension rate matrices :

k=np.loadtxt(rates_directory+'rate0.txt')		# Primer extension rate table after a match (Niraja's paper)
k1=np.loadtxt(rates_directory+'rate1.txt')		# Primer extension rate table after a mismatch (Chen's paper)

# Reaction rates after a match :

def rate_match(primer,template,i):		# This function checks the templating base, and assigns incorporation rates from Niraja's paper (after match)
	l=len(primer[i])
	if template[i][l]=='A':
		a=k[0]
	if template[i][l]=='U':
		a=k[1]
	if template[i][l]=='G':
		a=k[2]
	if template[i][l]=='C':
		a=k[3]
	return a

# Reaction rates after a mismatch :

def rate_mismatch(primer,template,i):		# This function checks the templating base, and assigns incorporation rates from Chen's paper (after mismatch)
	l=len(primer[i])
	if template[i][l]=='A':
		a=k1[0]
	if template[i][l]=='U':
		a=k1[1]
	if template[i][l]=='G':
		a=k1[2]
	if template[i][l]=='C':
		a=k1[3]
	return a


############################### Subroutines #####################################

# Function to detect mismatches:

def mismatch(primer,template,i):
	l=len(primer[i])
	if complement[primer[i][l-1]]!=template[i][l-1]:
		return 1
	else:
		return 0
	

# Choosing templates from the strands based on their free energies and calculate abundance of secondary structures :

def template_choice_and_structure_check(strand):

	hairpin,double_hairpin,hammerhead,cloverleaf,more_complex=0.0,0.0,0.0,0.0,0.0
	
	indices,FE=[],[]
	for i in xrange (0,N):

		# Secondary Structure detection :
		
		(s, E)=RNA.fold(strand[i])
		s=s.replace('.','')
		s=[[c,len(list(cgen))] for c,cgen in groupby(s)]
		l=len(s)
		if l==4 and s[1][0]==')' and s[2][0]=='(' and s[0][1]>s[1][1] and s[3][1]>s[2][1] and (s[0][1]-s[1][1])==(s[3][1]-s[2][1]):
			hammerhead=hammerhead+1.0
		if l==4 and s[1][0]==')' and s[2][0]=='(' and s[0][1]==s[1][1] and s[3][1]==s[2][1]:
			double_hairpin=double_hairpin+1.0
		if l==2:
			hairpin=hairpin+1.0
		if l==6:
			cloverleaf=cloverleaf+1.0
		if l>6:
			more_complex=more_complex+1.0
		
		# Templating efficiency calculation and choosing templates based on templating efficiency :

		if random.uniform(0,1)<=exp(E/2.0) and len(strand[i])>=10:
			indices=indices+[i]
		FE=FE+[E]

	return indices,hairpin,double_hairpin,hammerhead,cloverleaf,more_complex,FE


# Generate initial primers :

def initial_primers(template):
	N=len(template)
	primer=['']*N
	for i in xrange (0,N):
		a=rate_match(primer,template,i)
		nucleotide=np.random.choice(monomers, 1, p=(a/sum(a)))[0]
		primer[i]=primer[i]+nucleotide
	return primer,N


# Template directed primer extension :

def primer_extension(primer,template,i):
	t=0.0
	while t<=semi_wet and len(primer[i])<len(template[i]):
		j=mismatch(primer,template,i)
		if j==0:
			a=rate_match(primer,template,i)
		if j==1:
			a=rate_mismatch(primer,template,i)
		a0=sum(a)
		dt=np.random.exponential(1.0/a0)
		t=t+dt
		if t<=semi_wet and len(primer[i])<len(template[i]):
			nucleotide=np.random.choice(monomers, 1, p=(a/a0))[0]
			primer[i]=primer[i]+nucleotide
	return primer[i]


# Hydrolytic cleavage of single strands :

def single_strand_hydrolysis(strand,t,i):

	# Assign cleaving probability to each bond depending on the secondary structure :
			
	(s, E)=RNA.fold(strand[i])
	s=s.replace('(','|')
	s=s.replace(')','|')
	indices=range(1,len(s))
	a=np.zeros(len(indices), dtype=float)
	broken=[]
	for j in indices:
		if (s[j-1]=='.' and s[j]=='|') or (s[j-1]=='|' and s[j]=='.'):
			a[j-1]=p_up
		if s[j-1]=='.' and s[j]=='.':
			a[j-1]=p_uu
		if s[j-1]=='|' and s[j]=='|':
			a[j-1]=p_pp
		
	while t[i]<=wet:
		a0=k_hyd*sum(a)
		if a0>0.0:
			dt=np.random.exponential(1.0/a0)
			t[i]=t[i]+dt
			if t[i]<=wet:

				# Choose cleaving bond :
					
				l=np.random.choice(indices, 1, p=(a/sum(a)))[0]
				broken=broken+[l]
				j=indices.index(l)
				del indices[j]
				a=np.delete(a, j)
		else:
			t[i]=wet+1.0

	# Break the strand at those bonds :

	l=len(broken)
	if l>0:
		broken.sort()
		for j in xrange (1,l):
			if len(strand[i][broken[j-1]:broken[j]])>1:
				strand=strand+[strand[i][broken[j-1]:broken[j]]]
		if len(strand[i][broken[l-1]:len(strand[i])])>1:
			strand=strand+[strand[i][broken[l-1]:len(strand[i])]]	
		if len(strand[i][0:broken[0]])>1:
			strand[i]=strand[i][0:broken[0]]
		else:
			strand[i]='X'

	return strand


# Hydrolysis of dangles of template-primer pairs :

def template_primer_hydrolysis(primer,template,strand,t,i):
	t1=0.0
	broken1,broken2=[],[]
	indices1,indices2=range(1,len(primer[i])),range(1,len(primer[i]))
	while t1<=wet:
		k_hyd_paired=p_pp*k_hyd*len(indices1+indices2)
		if len(template[i])==len(primer[i]):
			k_hyd_dangle=0.0
		else:
			k_hyd_dangle=p_uu*k_hyd*(len(template[i])-len(primer[i])-1)+p_up*k_hyd
		prop=np.array([k_hyd_paired,k_hyd_dangle])
		a0=sum(prop)
		if a0>0.0:
			dt=np.random.exponential(1.0/a0)
			t1=t1+dt
			if t1<=wet:
				reaction=np.random.choice(['paired','dangle'], 1, p=(prop/a0))

				# If hydrolysis happens within the paired region :

				if reaction=='paired':

					# If the broken bond is in the template :

					if random.uniform(0,1)<=(1.0*len(indices1)/len(indices1+indices2)):
						l=random.choice(indices1)
						broken1=broken1+[l]
						j=indices1.index(l)
						del indices1[j]

					# If the broken bond is in the primer :
	
					else:
						l=random.choice(indices2)
						broken2=broken2+[l]
						j=indices2.index(l)
						del indices2[j]

				# If hydrolysis happens in the dangle :
			
				if reaction=='dangle':
					indices=range(len(primer[i]),len(template[i]))
					a=np.array([p_uu]*len(indices))
					a[0]=p_up
					l=np.random.choice(indices, 1, p=(a/sum(a)))[0]
					if len(template[i][l:])>1:
						strand=strand+[template[i][l:]]
						t=t+[t1]	
					template[i]=template[i][:l]

		else:
			t1=wet+1.0	

	# If any bond in the paired region is hydrolyzed, break the template and primer at those bonds :
	
	l=len(broken1)
	if l>0:
		broken1.sort()
		for j in xrange (1,l):
			if len(template[i][broken1[j-1]:broken1[j]])>1:
				template=template+[template[i][broken1[j-1]:broken1[j]]]
		if len(template[i][broken1[l-1]:len(template[i])])>1:
			template=template+[template[i][broken1[l-1]:len(template[i])]]	
		if len(template[i][0:broken1[0]])>1:
			template[i]=template[i][0:broken1[0]]
		else:
			template[i]='X'

	l=len(broken2)
	if l>0:
		broken2.sort()
		for j in xrange (1,l):
			if len(primer[i][broken2[j-1]:broken2[j]])>1:
				primer=primer+[primer[i][broken2[j-1]:broken2[j]]]
		if len(primer[i][broken2[l-1]:len(primer[i])])>1:
			primer=primer+[primer[i][broken2[l-1]:len(primer[i])]]	
		if len(primer[i][0:broken2[0]])>1:
			primer[i]=primer[i][0:broken2[0]]
		else:
			primer[i]='X'

	return primer,template,strand,t


# Concatenation of strands and monomers :

def concatenation(strand):
	i=random.choice(range(0,len(strand)))
	nucleotide=random.choice(monomers)
	strand[i]=nucleotide+strand[i]
	return strand


# Sampling function to select strands for next dry phase :

def sample(i):
	return (1.0-1.0/len(i)**0.588)


############################## Initial strands ##################################

strand=['']*N						# Array containing RNA strands

# Stacking : Generating initial homogeneous strands of different lengths

for i in xrange (0,N):
	strand[i]=random.choice(['A','G'])
	length=int(round(np.random.normal(8,1)))
	for l in xrange (1,length):
		strand[i]=strand[i]+strand[i][0]


############################# Time evolution ####################################

for day in xrange (0,Day_max):

	N=len(strand)

	# Export Avg. Length and Max length  of strands :
	
	length=[1.0*len(i) for i in strand]
	n=open(target_directory+'data_1.txt','a+')
	n.write('%d' %day + '\t' + '%f' %(sum(length)/N) + '\t' + '%f' %max(length) + '\t')
	n.close()
	length=[]
	
	# Export average GC content of strands :

	gc=[]
	for i in xrange (0,N):
		gc=gc+[1.0*(strand[i].count('G')+strand[i].count('C'))/len(strand[i])]

	n=open(target_directory+'data_3.txt','a+')
	n.write('%d' %day + '\t' + '%f' %(sum(gc)/N) + '\t')
	n.close()
	gc=[]

	# Choosing templates from the strands based on their free energies and calculate abundance of secondary structures :
	
	indices,hairpin,double_hairpin,hammerhead,cloverleaf,more_complex,FE=template_choice_and_structure_check(strand)

	# Export Secondary structure abundances :
	
	m=open(target_directory+'data_2.txt','a+')		
	m.write('%d' %day + '\t' + '%f' %(hairpin/N) + '\t' + '%f' %(double_hairpin/N) + '\t' + '%f' %(hammerhead/N) + '\t' + '%f' %(cloverleaf/N) + '\t' + '%f' %(more_complex/N) + '\n')
	m.close()
	
	# Export avg. free energy of the strands :
	
	n=open(target_directory+'data_1.txt','a+')
	n.write('%f' %(sum(FE)/N) + '\t')
	n.close()
	FE=[]
	
	# Separate templates from the strands :

	template=[strand[i] for i in indices]

	# Remaining Non-template strands :	

	strand=[strand[i] for i in range(0,N) if i not in indices]
	indices=[]

	# Export fraction of templates among the strands :
	
	n=open(target_directory+'data_1.txt','a+')
	n.write('%f' %(1.0*len(template)/N) + '\t')
	n.close()

	# Export maximum length of templates :

	length=[1.0*len(i) for i in template]
	n=open(target_directory+'data_1.txt','a+')
	n.write('%f' %max(length) + '\t')
	n.close()
	length=[]
	


	############# Semi-wet phase ############################	

	# Template directed Primer extension :

	# Initialize primers of length 1 nucleotide :
	
	primer,N=initial_primers(template)
	
	# Primer extension :
	
	for i in xrange (0,N):
		primer[i]=primer_extension(primer,template,i)


	# Export Avg. Length and Max length  of primers :
	
	length=[1.0*len(i) for i in primer]
	n=open(target_directory+'data_1.txt','a+')
	n.write('%f' %(sum(length)/N) + '\t' + '%f' %max(length) + '\n')
	n.close()
	length=[]
	
	# Export average errors (or mismatches) in template-primer pairs :

	avg_error=[]
	for i in xrange (0,N):
		l=len(primer[i])
		if l>=10:
			error=0.0
			for j in xrange (0,l):
				if primer[i][j]!=complement[template[i][j]]:
					error=error+1.0
			avg_error=avg_error+[error/l]

	if len(avg_error)>0:
		n=open(target_directory+'data_3.txt','a+')
		n.write('%f' %(sum(avg_error)/len(avg_error)) + '\n')
		n.close()
		avg_error=[]
	else:
		n=open(target_directory+'data_3.txt','a+')
		n.write('%f' %(len(avg_error)) + '\n')
		n.close()

				
	################## Wet phase ############################	

	# Cleavage of dangles of template-primer pairs :
	
	t=[0.0]*len(strand)
	N=len(primer)
	for i in xrange (0,N):
		primer,template,strand,t=template_primer_hydrolysis(primer,template,strand,t,i)
	
	template=[i for i in template if i!='X']
	primer=[i for i in primer if i!='X']	


	# Hydrolytic cleavage of single strands :
	
	N=len(strand)
	for i in xrange (0,N):
		strand=single_strand_hydrolysis(strand,t,i)

	strand=[i for i in strand if i!='X']
	t=[]
	


	######################## Dry-phase ######################

	# Reject primers of length = 1 as they are just monomers and reverse their 5'-3' direction :
	
	N=len(primer)
	for i in xrange (0,N):
		if len(primer[i])==1:
			primer[i]='X'
		else:
			primer[i]=primer[i][::-1]

	primer=[i for i in primer if i!='X']

	strand=strand+template+primer

	# Sampling of strands based on their lengths :
		
	if len(strand)>N_max:
		a=np.array([sample(i) for i in strand])
		strand=list(np.random.choice(strand, N_max, p=(a/sum(a)), replace=False))

	# Concatenation and dimerization :

	t=0.0
	while t<=dry:
		a0=4.0*k_con*len(strand)
		dt=np.random.exponential(1.0/a0)
		t=t+dt
		if t<=dry:
			strand=concatenation(strand)

			# Dimer formation at each time step:		

			if random.uniform(0,1)<=16.0*k_con*dt:
				strand=strand+[''.join(np.random.choice(monomers, 2))]


			
#################################### End ########################################
