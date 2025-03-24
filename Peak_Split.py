# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 18:16:10 2023

@author:  annap
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import os

file_directory = './'
file_name = 'Exp_Data'
file_type = 'txt'
new_folder = file_name+'_results'

# %% Preliminary conditions

w_L = 639 # lower wavelength in nm
peak_divs = [2]## defines how many peaks you want to use, min number of peaks is 2 to average TE and TM. If you want 2 peaks you would have one value in peak_divs that is inbetween the peaks.
w_U = 641 # upper wavelength in nm
Threshold = 5 #noise lvl
Threshold_2 = 5 # the amount of counts the peaks need to be above to after signal processing to analyse
Threshold_3 = 1 #nm (110pm) this is the threshold of the max difference of peak splits -  it is the resolution of the spectrometer
peak_percent = 0.1 # The percentage at which you want to analyse the peak as a decimal from 0-1
optical_sys_res = 0.05  ## in nm resolution (FWHM) of the optical system calculated by the andor website
Graph_dpi = 600

# %% Finding the optical resolution of the optical system to subtract this from each peak splitting

Gauss_1 = (2*np.sqrt(2*np.log(2)))**2
Gauss_2 = 2*(optical_sys_res**2)
Gauss_3 = -Gauss_2*np.log(peak_percent)
Gauss_4 = np.sqrt((Gauss_3)/(Gauss_1))
Optical_sys_res_FW_at_peak_percentage = 2*Gauss_4 ## in nm
# Optical_sys_res_FW_at_peak_percentage = 0 ## in nm


exp_spec_load1 = np.loadtxt(file_directory+'\\'+file_name + '.'+file_type)[:,:]
exp_spec_load = exp_spec_load1[:,:]
## checkpoint whole spectra plotting
# =============================================================================
# plt.plot(exp_spec_load[:,0],exp_spec_load[:,1:]) ## raw data; X-first column, Y-data
# plt.xlabel('Wavelength ($nm$)')
# plt.ylabel('Counts') 
# plt.savefig(file_directory+'\\'+file_name+'_Graph1.png',dpi = Graph_dpi,bbox_inches='tight')
# plt.clf()
# =============================================================================

w_L_index = np.argmin(np.abs(exp_spec_load[:,0]-w_L)) #index of the w_L in the main file
w_U_index = np.argmin(np.abs(exp_spec_load[:,0]-w_U)) #index of the w_U in the main

exp_spec = exp_spec_load[w_L_index:w_U_index,:] #select only specified wavelengths

#%%  checkpoint if spectra was properly selected

# =============================================================================
# plt.plot(exp_spec[:,0],exp_spec[:,1:]) ## trimming data
# plt.xlabel('Wavelength ($nm$)')
# plt.ylabel('Counts') 
# plt.savefig(file_directory+'\\'+file_name+'_Graph2.png',dpi = Graph_dpi,bbox_inches='tight')
# plt.clf()
# =============================================================================

# %% Spectra selection

#Frames only with signal above threshold

Frames_to_keep = []
for i in range(len(exp_spec[0,:])):
    if np.any(exp_spec[:,i] > Threshold):
        Frames_to_keep.append(i)
    else:
        pass

#Manual selection
# Frames_to_keep = []
# number_of_frames = 7
# ftkaa = np.linspace(0,number_of_frames,number_of_frames+1)
# for i in range(len(ftkaa)):
#     Frames_to_keep.append(int(ftkaa[i]))

exp_spec_final = np.zeros((len(exp_spec[:,0]),len(Frames_to_keep)))
for i in range(len(Frames_to_keep)):
    exp_spec_final[:,i] = exp_spec[:,Frames_to_keep[i]]
    
peak_divs_index = [0]*len(peak_divs)
for i in range(len(peak_divs)):
    peak_divs_index[i] = np.argmin(np.abs(exp_spec_final[:,0]-peak_divs[i]))    
div = [w_L_index] + peak_divs_index + [w_U_index]

## checkpoint
# =============================================================================
# plt.plot(exp_spec_final[:,0],exp_spec_final[:,1:]) ## removal of non lasing data
# plt.xlabel('Wavelength ($nm$)')
# plt.ylabel('Counts') 
# plt.savefig(file_directory+'\\'+file_name+'_Graph3.png',dpi = Graph_dpi,bbox_inches='tight')
# plt.clf()
# =============================================================================

#%% Fitting

FtR_all = []
spec_fit = np.zeros((40,(len(exp_spec_final[0,:]))))
spec_split_w = np.zeros((len(exp_spec_final[0,:])-1,(len(div)-1)*2))
Sav_Go_Filter = np.zeros((len(exp_spec_final[:,0]),len(exp_spec_final[0,:])))
Sav_Go_Filter[:,0] = exp_spec_final[:,0] #keeps the wavelengths (x axis)

back_sub = np.zeros((len(exp_spec_final[:,0]),len(exp_spec_final[0,:])))
back_sub[:,0] = exp_spec_final[:,0]

fitted_lines = np.zeros((4,4)) ## only for 1 spectrum 
##########################################################################################################################

for i in range(len(exp_spec_final[0,:])):
    if i == 0:
        spec_fit[0:20,i] = exp_spec_final[0:20,i]
        spec_fit[20:,i] = exp_spec_final[len(exp_spec_final[:,0])-20:,i]
        
    else:
        
        ##### Background removal 
        
        spec_fit[0:20,i] = exp_spec_final[0:20,i]
        spec_fit[20:,i] = exp_spec_final[len(exp_spec_final[:,0])-20:,i]
        
        fit = np.polyfit(spec_fit[:,0],spec_fit[:,i],1) #[m,c]
        exp_spec_final[:,i] = (exp_spec_final[:,i] - fit[1] - (fit[0]*exp_spec_final[:,0]))
        back_sub[:,i] = exp_spec_final[:,i]

        #### Filtering
        ## Savitzky-Golay filter to smooth the data remove higher frequencies
        exp_spec_final[:,i] = savgol_filter(exp_spec_final[:,i], 3, 1) # The Savitsky-golay filter is applicable for split peaks and doesn't affect the shape but for very sharp undeformed peaks it isn't suitable it deforms the peak shape.
        Sav_Go_Filter[:,i] = exp_spec_final[:,i]

        #####
        count = 0
        for j in range(len(div)-1):
            
            split_peak_I = exp_spec_final[div[j]:div[j+1],i]
            split_peak_I_r = split_peak_I[::-1]
            
            #### I_ref using percentage of max intensity 
            I_max = np.max(split_peak_I)
            I_max_index = np.argmax(split_peak_I)
            I_ref = I_max*peak_percent
            
            #### I_ref from background level
#            
#            I_bgrd1 = np.max(split_peak_I[0:50])
#            I_bgrd2 = np.max(split_peak_I[-50:])
#            if I_bgrd1 > I_bgrd2:
#                I_ref = I_bgrd1*2
#            else:
#                I_ref = I_bgrd2*2
#            
        #####
            
            k = 0
            x = split_peak_I[k] - I_ref
            while x < 0:
                k = k + 1
                x = split_peak_I[k] - I_ref

            #### Linear spacing
            w1 = exp_spec_final[k+div[j],0]
            w2 = exp_spec_final[k+1+div[j],0]
            w_space1 = np.linspace(w1,w2,1000)
            I1 = exp_spec_final[k+div[j],i]
            I2 = exp_spec_final[k+1+div[j],i]
            I_space1 = np.linspace(I1,I2,1000)
            I_index1 = np.argmin(np.abs(I_space1 - I_ref))
            spec_split_w[i-1,count] = w_space1[I_index1]
            a = w_space1[I_index1]
            II1 = I_space1[I_index1]
#            
            ####
            
            count = count + 1
        #####
            l = 0
            y = split_peak_I_r[l] - I_ref
            while y < 0:
                l = l + 1
                y = split_peak_I_r[l] - I_ref

            #### Linear spacing
            w3 = exp_spec_final[len(split_peak_I_r) - (l + 1) + div[j],0]
            w4 = exp_spec_final[len(split_peak_I_r) - l + div[j],0]
            w_space2 = np.linspace(w3,w4,10000)
            I3 = exp_spec_final[len(split_peak_I_r) - (l + 1) + div[j],i]
            I4 = exp_spec_final[len(split_peak_I_r) - l + div[j],i]
            I_space2 = np.linspace(I3,I4,10000)
            I_index2 = np.argmin(np.abs(I_space2 - I_ref))
            spec_split_w[i-1,count] = w_space2[I_index2]
            b = w_space2[I_index2]
            II2 = I_space2[I_index2]
            ####    
            
        
            count = count + 1
        #####
        ## Resolving if all peaks in a spectra have a intensity greater than Threshold 2
            if I_max < Threshold_2:
                FtR_all.append(i)
                
        
        ##### 
        ## visulising the fitted peak splits.
#            W_fit_x = np.array([a,b])
#            W_fit_y = np.array([II1,II2])            
            W_fit_x = np.array([a,b])
            W_fit_y = np.array([I_ref,I_ref])
#            W_fit_x = np.array([ exp_spec_final[index_l,0],exp_spec_final[index_u,0] ])
#            W_fit_y = np.array([exp_spec_final[index_l,i],exp_spec_final[index_u,i]])
            
            
        ######   
            fitted_lines[j,0:2] = W_fit_x
            fitted_lines[j,2:] = W_fit_y
            ## here we would need to plot the individual graphs for 'peak fit'.
            plt.plot(exp_spec_final[div[j]:div[j+1],0],exp_spec_final[div[j]:div[j+1],i],color='r') # Raw data 
            plt.plot(W_fit_x , W_fit_y,marker = '.',color='b') # percentage of peak line

        plt.xlabel('Wavelength ($nm$)')
        plt.ylabel('Counts') 
        plt.savefig(file_directory+'\\'+file_name+'_Frame_'+str(Frames_to_keep[i])+'.png',dpi = Graph_dpi,bbox_inches='tight')
        plt.clf()
        
        
        ######
        
        
np.savetxt(file_directory+'\\'+file_name+'_fitted_lines_nm.txt',fitted_lines)     
        
### plotting the background subtraction        

# plt.plot(back_sub[:,0],back_sub[:,1:]) ## background removal
# plt.xlabel('Wavelength ($nm$)')
# plt.ylabel('Counts') 
# plt.savefig(file_directory+'\\'+file_name+'_Graph4.png',dpi = Graph_dpi,bbox_inches='tight')
# plt.clf()
# np.savetxt(file_directory+'\\'+file_name+'_BackgroundSub_nm.txt',back_sub)     
   
## plotting the savitsky golay filter
# plt.plot(Sav_Go_Filter[:,0],Sav_Go_Filter[:,1:]) ## Smoothing of data
# plt.xlabel('Wavelength ($nm$)')
# plt.ylabel('Counts') 
# plt.savefig(file_directory+'\\'+file_name+'_Graph5.png',dpi = Graph_dpi,bbox_inches='tight')
# plt.clf()
# np.savetxt(file_directory+'\\'+file_name+'_Sav-GoFilter_nm.txt',Sav_Go_Filter)


#
####
## Removing the frames that we don't want, for the Frames_to_keep and for the matrix have the peaks in.

FtR = []            
for i in FtR_all: 
    if i not in FtR: 
        FtR.append(i)

FtK_indices =range(len(Frames_to_keep))
FtK_indices = np.delete(FtK_indices,FtR)
Frames_to_keep = np.delete(Frames_to_keep,FtR)
 
         
exp_spec_final_1 = np.zeros((len(exp_spec_final[:,0]),len(Frames_to_keep))) 
spec_split_w_1 = np.zeros((len(FtK_indices)-1,len(spec_split_w[0,:]))) 
for i in range(len(FtK_indices)-1):
        spec_split_w_1[i,:] = spec_split_w[FtK_indices[i],:]
for i in range(len(FtK_indices)):
        exp_spec_final_1[:,i] = exp_spec_final[:,FtK_indices[i]]

# del FtR,FtK_indices,exp_spec_final,spec_split_w,Threshold_2,FtR_all,W_fit_x,W_fit_y
#####

spec_split_abs_minOpRes = np.zeros((len(exp_spec_final_1[0,:])-1, int(np.ceil(len(spec_split_w_1[0,:])/2))))
spec_split_abs = np.zeros((len(exp_spec_final_1[0,:])-1, int(np.ceil(len(spec_split_w_1[0,:])/2))))
spec_split_avg = np.zeros(len(exp_spec_final_1[0,:])-1)
spec_split_avg_minOpRes = np.zeros(len(exp_spec_final_1[0,:])-1)
spec_split_T3 = np.zeros(((len(exp_spec_final_1[0,:])-1),2))

for i in range(len(spec_split_w_1[:,0])):
    count = 0
    for j in range(len(spec_split_abs[0,:])):
        k = j + count
        spec_split_abs[i,j] = - spec_split_w_1[i,k] +  spec_split_w_1[i,k+1]
        spec_split_abs_minOpRes[i,j] = - spec_split_w_1[i,k] +  spec_split_w_1[i,k+1]  - Optical_sys_res_FW_at_peak_percentage
        count = count + 1
 
    spec_split_avg[i] = np.mean(spec_split_abs[i,:])
    spec_split_avg_minOpRes[i] = np.mean(spec_split_abs[i,:]) - Optical_sys_res_FW_at_peak_percentage
    
    spec_split_T3[i,0] = np.max(spec_split_abs[i,0::2]) - np.min(spec_split_abs[i,0::2])#even # need to change this to be the std or another analytical value for TE and TM peaks
    
    if peak_divs == []:
        continue
    else:        
        spec_split_T3[i,1] = np.max(spec_split_abs[i,1::2]) - np.min(spec_split_abs[i,1::2])#odd # need to change this to be the std or another analytical value for TE and TM peaks
    
del i,j,spec_fit,fit,div,I_max,I_max_index,I_ref,split_peak_I, split_peak_I_r,count,w1,w2,w3,w4,w_space1,w_space2,I_space1,I_space2,I1,I2,I3,I4,x,y,I_index1,I_index2

#### standard deviation filter

Frames_to_keep = Frames_to_keep[1:]

delete = []
for i in range(len(spec_split_T3[:,0])):
    if spec_split_T3[i,0] > Threshold_3:
        delete.append(i)
    elif spec_split_T3[i,1] > Threshold_3:
        delete.append(i)
    

del Threshold_3,a,b,i,k,l

spec_split_w_1 = np.delete(spec_split_w_1,delete,axis = 0)
spec_split_abs = np.delete(spec_split_abs,delete,axis = 0)
spec_split_avg = np.delete(spec_split_avg,delete)
spec_split_T3 = np.delete(spec_split_T3,delete,axis = 0)
Frames_to_keep = np.delete(Frames_to_keep,delete)

del delete

np.savetxt(file_directory+'\\'+file_name+'_pk_nm_'+str(peak_percent)+'-'+str(optical_sys_res)+'.txt',spec_split_w_1) # spectrum of the peak widths defined 
np.savetxt(file_directory+'\\'+file_name+'_pk_abs_nm_'+str(peak_percent)+'-'+str(optical_sys_res)+'.txt',spec_split_abs) # THe splitting for each peak 
np.savetxt(file_directory+'\\'+file_name+'_pk_abs_minOpRes_nm_'+str(peak_percent)+'-'+str(optical_sys_res)+'.txt',spec_split_abs_minOpRes) # the real splitting for each peak minus the optical resolution of the system
np.savetxt(file_directory+'\\'+file_name+'_pk_avg_nm_'+str(peak_percent)+'-'+str(optical_sys_res)+'.txt',spec_split_avg) # The average splitting across all peaks in the spectrum for every time point
np.savetxt(file_directory+'\\'+file_name+'_pk_avg_minOpRes_nm_'+str(peak_percent)+'-'+str(optical_sys_res)+'.txt',spec_split_avg_minOpRes) # The average splitting across all peaks in the spectrum for every time point
# np.savetxt(file_directory+'\\'+file_name+'_peaks_std_nm_pp_'+str(peak_percent)+'-'+str(optical_sys_res)+'.txt',spec_split_T3)
np.savetxt(file_directory+'\\'+file_name+'_pk_Time_'+str(peak_percent)+'-'+str(optical_sys_res)+'.txt',Frames_to_keep) ## This lets you know what frames have been discarded because the haven't met the thresholds set at teh start/

# plt.plot(exp_spec_final_1[:,0],exp_spec_final_1[:,1:]) ## After second lasing threshold and a standard deviation filter
# plt.xlabel('Wavelength ($nm$)')
# plt.ylabel('Counts') 
# plt.savefig(file_directory+'\\'+file_name+'_Graph6.png',dpi = 250,bbox_inches='tight')





#############################################
Time = datetime.now() - startTime
print('----------------------------------')
print('Time = ', Time)
print('----------------------------------')






# np.savetxt(file_directory+'\\'+file_name+'_peak4_nm.txt',spec_split_w_1)
# np.savetxt(file_directory+'\\'+file_name+'_peak4_abs_nm.txt',spec_split_abs)
# np.savetxt(file_directory+'\\'+file_name+'_peak4_avg_nm.txt',spec_split_avg)
# # np.savetxt(file_directory+'\\'+file_name+'_peak1_std_nm.txt',spec_split_T3)
# np.savetxt(file_directory+'\\'+file_name+'_peak4_Time.txt',Frames_to_keep)



