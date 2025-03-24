# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 13:16:57 2023

@author:  annap
"""


from datetime import datetime
startTime = datetime.now()

from Functions import eccentricity_solver as ES
import numpy as np
import matplotlib.pyplot as plt 
import sys

file_directory = '........'
file_name = '.......'
file_type = 'txt'

exp_spec_nm = np.loadtxt(file_directory+'\\'+file_name+'.'+file_type)*10**(-9) # input in nm 

########################################################################################

##### Inputs ####

#a_range = np.array([2.25,2.75])*10**(-6) # estimated radius from microscopy image
b_dev = 50*10**(-9) # won't need changing
l = np.linspace(10,1000,991) # Mode number
n_cell_range = np.array([1.333,1.41]) # refractive index of the cell 
n_bead = 1.63 # refractive index of bead 
b_steps = 10000
gamma = 37*10**(-3) # Surface tension


peak_percent = 0.1 ##  This is the percentage of max-intensity that the peaks were analysed at.



optical_sys_res = 0.05  ## in nm resolution (FWHM) of the optical system calculated by the andor website


### Finding the optical resolution of the optical system to subtract this from each peak splitting
Gauss_1 = (2*np.sqrt(2*np.log(2)))**2
Gauss_2 = 2*(optical_sys_res**2)
Gauss_3 = -Gauss_2*np.log(peak_percent)
Gauss_4 = np.sqrt((Gauss_3)/(Gauss_1))
Optical_sys_res_FW_at_peak_percentage = (2*Gauss_4)*10**(-9) ## in nm
Optical_sys_res_FW_at_peak_percentage = 0 # when simulating data







# should be an even amount of pairs so that it is averaged evenly over TE and TM
if exp_spec_nm.ndim == 1:
    exp_spec_nm = exp_spec_nm[None,:]

### Getting an estimate for a_range using the FSR ####

### If there is not enough peaks
if len(exp_spec_nm[0,:])/2 < 4:
    sys.exit('Not enough peaks, need to manually enter an a_range.')
if len(exp_spec_nm[0,:])%2 != 0:
    sys.exit('Odd number of peaks, need to enter an even amount of TE and TM modes.')



fsr_all = []    
n_eff_range = [(n_bead*0.7+n_cell_range[0]*0.3),n_bead]
a_range = np.zeros((2,len(exp_spec_nm[:,0])))

for i in range(len(exp_spec_nm[:,0])):
    a_fsr_range = []
    for j in range(int(len(exp_spec_nm[0,:])/2)):
        
        
        fsr = exp_spec_nm[i,(j+4)] - exp_spec_nm[i,j] 
        fsr_all.append(fsr)
        lambda2 = (exp_spec_nm[i,(j+4)] - (Optical_sys_res_FW_at_peak_percentage/2))  * (exp_spec_nm[i,j]- (Optical_sys_res_FW_at_peak_percentage/2))
        
        r1 = lambda2/(2*np.pi*fsr*n_eff_range[0])   # r will = a
        r2 = lambda2/(2*np.pi*fsr*n_eff_range[1])   # r will = a
        
        a_fsr_range.append(r1)
        a_fsr_range.append(r2) 
       
    a_range[0,i] = np.min(a_fsr_range)
    a_range[1,i] = np.max(a_fsr_range)



##### MANUAL a_range for when not enough peaks are available


# l1= 626.48*1e-9
# l2= 638.02*1e-9
# l3= 632.29*1e-9
# l4= 620.89*1e-9

# r1min = l1*l2/(2*np.pi*(l2-l1)*n_eff_range[0])   # r will = a
# r2max = l3*l4/(2*np.pi*(l3-l4)*n_eff_range[1])   # r will = a

# a_range[1,:] = r1min*2
# a_range[0,:] = r2max*2

# a_range1 = [6.61863963,6.61863963,7.248605829,7.248605829] 
# a_range = np.zeros((2,2))
# a_range[0,:] = a_range1[0:2]
# a_range[1,:] = a_range1[2:]
# a_range = a_range*10**(-6)


# a_range= np.array([[8.07264007e-06, 8.05819975e-06, 8.07775313e-06, 8.07234671e-06,
#         8.07105221e-06, 8.07296764e-06, 8.06717738e-06, 8.06887046e-06,
#         8.05707437e-06, 8.07233948e-06, 8.05645337e-06, 8.07182799e-06,
#         8.07383210e-06, 8.07058701e-06, 8.07058701e-06, 8.07268486e-06,
#         8.07316719e-06, 8.06875189e-06, 8.06875189e-06, 8.07410737e-06,
#         8.05518178e-06, 8.04398148e-06, 8.04493462e-06, 8.05461752e-06,
#         8.06217572e-06, 8.05720380e-06, 8.06546573e-06, 8.06718581e-06],
#         [8.48410550e-06, 8.47040933e-06, 8.48410550e-06, 8.46899420e-06,
#         8.48333781e-06, 8.46861347e-06, 8.46750711e-06, 8.46750711e-06,
#         8.46807174e-06, 8.46804342e-06, 8.46600179e-06, 8.46632324e-06,
#         8.48160060e-06, 8.48577505e-06, 8.46503387e-06, 8.48643078e-06,
#         8.47900738e-06, 8.46323207e-06, 8.46804706e-06, 8.47966960e-06,
#         8.47698755e-06, 8.47698755e-06, 8.47337492e-06, 8.47337492e-06,
#         8.45881936e-06, 8.47208847e-06, 8.47208847e-06, 8.47275009e-06]])
    
#### EXPERIMENTAL SPECTRA ####

exp_spec_split = np.zeros((len(exp_spec_nm[:,0]),int((len(exp_spec_nm[0,:]))/2)))
exp_spec_sum = np.zeros((len(exp_spec_nm[:,0]),int((len(exp_spec_nm[0,:]))/2)))

for j in range(len(exp_spec_nm[:,0])):
    count = 0
    for i in range(int(len(exp_spec_nm[0,:])/2)):
        
        exp_spec_split[j,i] = -exp_spec_nm[j,count] + exp_spec_nm[j,count+1] - Optical_sys_res_FW_at_peak_percentage # Input split must always be positive!
        exp_spec_sum[j,i] = exp_spec_nm[j,count] + exp_spec_nm[j,count+1] - Optical_sys_res_FW_at_peak_percentage
        count = count+2

########################################################################################

e_2 = np.zeros((len(exp_spec_split[:,0]),len(exp_spec_split[0,:])))
e_error = np.zeros((len(exp_spec_split[:,0]),len(exp_spec_split[0,:])))
modes_L = np.zeros((len(exp_spec_split[:,0]),len(exp_spec_split[0,:]),2))
modes_U = np.zeros((len(exp_spec_split[:,0]),len(exp_spec_split[0,:]),2))
count = 0
for i in range(len(exp_spec_split[:,0])):
   for j in range(len(exp_spec_split[0,:])):
        count = count+1
        # result = ES(np.array([a_range[0,i],a_range[1,i]]),b_dev,l,n_cell_range,n_bead,b_steps,exp_spec_sum[i,j],exp_spec_split[i,j],count)
        result = ES(a_range[:,i],b_dev,l,n_cell_range,n_bead,b_steps,exp_spec_sum[i,j],exp_spec_split[i,j],count)
        e_2[i,j] = result[0]
        e_error[i,j] = result[1]
        modes_L[i,j,:] = result[2]
        modes_U[i,j,:] = result[3]
    
    
e_2_avg = np.mean(e_2,axis=1)
error_1 = np.square(e_error)
error_2 = np.sum(error_1,axis=1)
error_3 = np.sqrt(error_2)
e_2_error_avg = error_3/(np.sqrt(len(e_error[0,:])))


########################################################################################

#### a_estimation Calculation ####

## Error ith a abnd then no error in b allows us to have the same error in stress and force 
## Will need to look at the difference between the geometries and see if it is large?

a = np.zeros(len(a_range[0,:]))
a_error = np.zeros(len(a_range[0,:]))
for i in range(len(a_range[0,:])):
    
    a[i] = a_range[0,i] + ((a_range[1,i]-a_range[0,i])/2)
    a_error[i] = ((a_range[1,i]-a_range[0,i])/2)/a[i]

#### Absolute errors ####

a_error_abs = a*a_error
e_2_error_avg_abs = e_2_error_avg*e_2_avg

#### Force Measurment ####
x = 1-e_2_avg

## F prolate

f_pro = 3*np.pi*gamma*a*( (x**(-1/2)) - (x**(-1/6)))

dfp_da = 3*np.pi*gamma*((x**(-1/2)) - (x**(-1/6)))
dfp_de = -np.pi*gamma*np.sqrt(e_2_avg)*a*((x**(1/3)) - 3)*(x**(-3/2))

df_pro_abs = np.sqrt( (dfp_da * a_error_abs)**2 + (dfp_de * e_2_error_avg_abs)**2   ) # absolute error in force oblate
df_pro = df_pro_abs/f_pro


## F oblate

f_ob = 3*np.pi*gamma*a*( (x**(1/6)) - (x**(1/2)))

dfo_da = 3*np.pi*gamma*((x**(1/6)) - (x**(1/2)))
dfo_de = np.pi*gamma*np.sqrt(e_2_avg)*a*((3*(x**(1/3))) - 1)*(x**(-5/6))

df_ob_abs = np.sqrt( (dfo_da * a_error_abs)**2 + (dfo_de * e_2_error_avg_abs)**2   ) # absolute error in force prolate
df_ob = df_ob_abs/f_ob

## Force_total and force_total error

F = (f_pro+f_ob)/2

F_error_abs  = np.zeros(len(F))
for i in range(len(F)):
    F_error_abs[i]  = (np.max([f_pro[i]+df_pro_abs[i],f_ob[i]+df_ob_abs[i]]) - np.min([f_pro[i]-df_pro_abs[i],f_ob[i]-df_ob_abs[i]]))/2
F_error = F_error_abs/F



#### Calculating r and db ####

b_ob = ((a**2) * (1-e_2_avg))**(1/2)
b_pro = ((a**2) / (1-e_2_avg))**(1/2)

r_ob = ((a**2)*b_ob)**(1/3)
r_pro = ((a**2)*b_pro)**(1/3)
r = (r_ob+r_pro) /2

delta_b_ob = abs(r_ob-b_ob)
delta_b_pro = abs(b_pro-r_pro)
delta_b = (delta_b_ob +delta_b_pro)/2





print('----------------------------------')
print('Force (nN) = ',F*10**(9))
#print('----------------------------------')
print('% Error = ',F_error*100)
print('Abs Error (N) = ',F_error_abs)

print('----------------------------------')
print('Eccentricity^2  = ',e_2_avg**(1/2))
#print('----------------------------------')
print('% Error = ',e_2_error_avg*0.5*100)
print('Abs Error (um) = ',(e_2_error_avg*0.5*(e_2_avg**(1/2))))

print('----------------------------------')
print('a axis (um) = ',a*10**(6))
#print('----------------------------------')
print('% Error = ',a_error*100)
print('Abs Error (um) = ',a_error*a*10**(6))

exp_spec_split_avg = np.mean(exp_spec_split,axis = 1)
print('----------------------------------')
print('Split (m) = ',exp_spec_split_avg)
print('----------------------------------')
print('Gamma (N/m) = ', gamma )

print('----------------------------------')
print('Delta b (m) = ',delta_b)

print('----------------------------------')
print('Radius (m) = ',r)



np.savetxt(file_directory+'\\'+file_name+'_Force.txt',F)
np.savetxt(file_directory+'\\'+file_name+'_Force_error_abs.txt',F_error_abs)
np.savetxt(file_directory+'\\'+file_name+'_Force_error.txt',F_error)
np.savetxt(file_directory+'\\'+file_name+'_Split_average.txt',exp_spec_split_avg)
np.savetxt(file_directory+'\\'+file_name+'_Split.txt',exp_spec_split)
np.savetxt(file_directory+'\\'+file_name+'_a.txt',np.array([a]))
np.savetxt(file_directory+'\\'+file_name+'_a_error.txt',np.array([a_error]))
np.savetxt(file_directory+'\\'+file_name+'_Gamma.txt',np.array([gamma]))
np.savetxt(file_directory+'\\'+file_name+'_Delta-b.txt',np.array([delta_b]))
np.savetxt(file_directory+'\\'+file_name+'_Radius.txt',np.array([r]))
np.savetxt(file_directory+'\\'+file_name+'_eccentricity.txt',np.array([e_2_avg]))


########################################################################################
Time = datetime.now() - startTime
print('----------------------------------')
print('Time = ', Time)
print('----------------------------------')
















































##### Simulated Spectra #####

#from TE_SvE import TE_P0_Pl
#from TM_SvE import TM_P0_Pl
#
#a_sim = 7.1*10**(-6)
#b_sim = np.array([7.103*10**(-6)])
#l_sim = [110,109,108]
#n_cell_sim = 1.36
#n_bead = 1.67
#
#exp_e = 1-(a_sim**2/b_sim**2)
#
#exp_spec_sum = np.zeros((1,2*len(l_sim))) 
#exp_spec_split = np.zeros((1,2*len(l_sim))) 
#
#count = 0
#for i in range(len(l_sim)):
#    
#    W_TE_p0_pl = TE_P0_Pl(a_sim,b_sim,l_sim[i],n_cell_sim,n_bead)
#    W_TM_p0_pl = TM_P0_Pl(a_sim,b_sim,l_sim[i],n_cell_sim,n_bead)
#    
#    exp_spec_sum[0,i+count] = W_TM_p0_pl[2] 
#    exp_spec_sum[0,i+1+count] = W_TE_p0_pl[2]    
#    exp_spec_split[0,i+count] = W_TM_p0_pl[0] 
#    exp_spec_split[0,i+1+count] = W_TE_p0_pl[0] 
#    count = count+1
#    
#
#
###### User Inputs #####
#
#a_range = np.array([7,7.5])*10**(-6)
#b_dev = 50*10**(-9) 
#l = np.linspace(10,300,291)
#n_cell_range = np.array([1.34,1.39])
#n_bead = 1.67
#b_steps = 10000
#
#exp_sum = np.abs(exp_spec_sum[0,0]) # l = 110
#exp_split = np.abs(exp_spec_split[0,0])



























