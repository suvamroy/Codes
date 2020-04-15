import os
import numpy as np


#basic_parameters

k_hyd=0.04						# Rate of hydrolytic cleavage				
k_con=0.62						# Rate of concatenation				
dry=12.0						# Duration of dry phase
semi_wet=8.0						# Duration of semi-wet phase
N=100							# Number of initial strands
N_max=1000						# Maximum number of strands
p_uu=1.0						# Hydrolysis Prefactor when both neighbor nucleotides are paired 
p_up=0.1						# Hydrolysis Prefactor when either one neighbor nucleotide is paired
p_pp=0.01						# Hydrolysis Prefactor when both neighbor nucleotides are unpaired
Day_max=200						# Total number of days
Length_of_the_day = 24.0

rates_directory='./'



#For example if you want to do two loops to vary both dry and semi wet ...

basic_directory='./data/'

#first we copy the code just for archive and then we do the loop over parameters ...
cmd='cp main.py '+basic_directory+'main.py.archive'
print(cmd)
os.system(cmd)

cmd='cp script_to_launch_multiple_runs.py '+basic_directory+'script_to_launch_multiple_runs.py.archive'
print(cmd)
os.system(cmd)


for dry in np.linspace(4,12,5):
	for semi_wet in np.linspace(4,12,5):
		if (dry==12 and semi_wet==12) or (dry==10 and semi_wet==12) or (dry==12 and semi_wet==10):
			pass
		else:
			target_directory=basic_directory+'dry='+str(dry)+'/semi_wet='+str(semi_wet)+'/'
			cmd = 'python2 main.py '+str(k_hyd)+' ' +str(k_con)+' ' +str(dry)+' ' +str(semi_wet)+' ' +str(N)+' ' +str(N_max)+' ' +str(p_uu)+' ' +str(p_up)+' ' +str(p_pp)+' ' +str(Day_max)+' ' +str(Length_of_the_day)+' ' +rates_directory+' '+target_directory +' &'
			print("starting a new job...")
			os.system(cmd)
			print(cmd)	
