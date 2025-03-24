# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:43:56 2020

@author: Joe1
"""

#%%
#####################################################################################################

def eccentricity_solver(a_range,b_dev,l,n_cell_range,n_bead,b_steps,exp_sum,exp_split,count):
       
    from Functions import TE_P0_Pl as TE
    from Functions import TM_P0_Pl as TM
    from Functions import Mode_Refining
    import numpy as np
    import matplotlib.pyplot as plt 
    #####################################################################################################
    
    
    
    #### Mode Refining ####
    
    modes = Mode_Refining(a_range,b_dev,l,n_cell_range,n_bead,exp_sum)
    l_ranges = modes[0]
    ##### a_lower
    modes_L = l_ranges[0]
    b_L = np.linspace(a_range[0]-b_dev,a_range[0]+b_dev,b_steps)
    
    ##### a_upper
    modes_U = l_ranges[1]
    b_U = np.linspace(a_range[1]-b_dev,a_range[1]+b_dev,b_steps)
        
    
    #####################################################################################################
    
    ##########################################################################################
#    #### Need to plot all the variable possibilities here to check the extremes then will hard
#    #### code these to improve speed
#    # can do te and tm for each step and also b_u and b_l in the same loop
#    modes_test = [modes_U,modes_L]
#    split_test = []
#    e_2_test = []
#    constraints_test = []
#    for i in range(len(a_range)):
#        for j in range(len(n_cell_range)):    
#            for k in range(len(modes_L)):
#                    if i == 0:
#                        TE_result = TE(a_range[i],b_L,modes_L[k],n_cell_range[j],n_bead)
#                        split_test.append(TE_result[0])
#                        e_2_test.append(TE_result[1])
#                        constraints_test.append([a_range[i],b_L,modes_L[k],n_cell_range[j],n_bead,'TE'])
#                        
#                        TE_result = TE(a_range[i],b_L,modes_U[k],n_cell_range[j],n_bead)
#                        split_test.append(TE_result[0])
#                        e_2_test.append(TE_result[1])
#                        constraints_test.append([a_range[i],b_L,modes_U[k],n_cell_range[j],n_bead,'TE'])
#                        
#                        TM_result = TM(a_range[i],b_L,modes_L[k],n_cell_range[j],n_bead)
#                        split_test.append(TM_result[0])
#                        e_2_test.append(TM_result[1])
#                        constraints_test.append([a_range[i],b_L,modes_L[k],n_cell_range[j],n_bead,'TM'])
#                        
#                        TM_result = TM(a_range[i],b_L,modes_U[k],n_cell_range[j],n_bead)
#                        split_test.append(TM_result[0])
#                        e_2_test.append(TM_result[1])
#                        constraints_test.append([a_range[i],b_L,modes_U[k],n_cell_range[j],n_bead,'TM'])
#                    else:
#                        TE_result = TE(a_range[i],b_U, modes_L[k],n_cell_range[j],n_bead)
#                        split_test.append(TE_result[0])
#                        e_2_test.append(TE_result[1])
#                        constraints_test.append([a_range[i],b_U,modes_L[k],n_cell_range[j],n_bead,'TE'])
#                        
#                        TE_result = TE(a_range[i],b_U,modes_U[k],n_cell_range[j],n_bead)
#                        split_test.append(TE_result[0])
#                        e_2_test.append(TE_result[1])
#                        constraints_test.append([a_range[i],b_U,modes_U[k],n_cell_range[j],n_bead,'TE'])
#                        
#                        TM_result = TM(a_range[i],b_U,modes_L[k],n_cell_range[j],n_bead)
#                        split_test.append(TM_result[0])
#                        e_2_test.append(TM_result[1])
#                        constraints_test.append([a_range[i],b_U,modes_L[k],n_cell_range[j],n_bead,'TM'])
#            
#                        TM_result = TM(a_range[i],b_U,modes_U[k],n_cell_range[j],n_bead)
#                        split_test.append(TM_result[0])
#                        e_2_test.append(TM_result[1])
#                        constraints_test.append([a_range[i],b_U,modes_U[k],n_cell_range[j],n_bead,'TM'])
#            
#    
#    for i in range(len(split_test[:])):
#   
#        plt.plot(split_test[i],e_2_test[i],linewidth = 0.1)               
#        plt.xlabel('Spliting ($nm$)')            
#        plt.ylabel('$e^2$')
#    plt.savefig('test_1.pdf')
#    plt.clf()
#    np.save('split_test'+str(count)+'.npy',split_test)
#    np.save('e_2_test'+str(count)+'.npy',e_2_test)
#    np.save('constraints_test'+str(count)+'.npy',constraints_test)
#    np.save('modes_test'+str(count)+'.npy',modes_test)
#    
#    
#    
    ################################################################################
    
    
    ##### e_2_lower

    split_L = np.zeros(len(b_L))
    e_2_L = np.zeros(len(b_L))
#    constraints_L = []
    
    TE_result = TE(a_range[1],b_U,modes_L[0],n_cell_range[1],n_bead)
#    constraints_L.append([a_range[0],b_L,modes_L[0],n_cell_range[1],n_bead,'TE'])
    split_L[:] = TE_result[0]
    e_2_L[:] = TE_result[1]
    
    ##### e_2_upper

    split_U = np.zeros(len(b_L))
    e_2_U = np.zeros(len(b_L))
#    constraints_U = []
    
    TM_result = TM(a_range[0],b_L,modes_U[1],n_cell_range[0],n_bead)
#    constraints_U.append([a_range[1],b_U,modes_U[0],n_cell_range[1],n_bead,'TM'])
    split_U[:] = TM_result[0]
    e_2_U[:] = TM_result[1]
    
    
    ######################################################################################################
#    ##### Plotting e_2 against split.
#    
#    #plt.plot(split_L[1,:]*10**9,e_2_L[1,:],label = 'Lower_Upper', linewidth = 0.25)
#    plt.plot(split_L[0,:]*10**9,e_2_L[0,:], label = 'Lower_Lower',linewidth = 0.25)
#    #plt.plot(split_U[0,:]*10**9,e_2_U[0,:],label = 'Upper_lower',linewidth = 0.25)
#    plt.plot(split_U[1,:]*10**9,e_2_U[1,:],label = 'Upper_upper', linewidth = 0.25)
##    plt.legend()
#    plt.xlabel('Split (nm)')
#    plt.ylabel('e**2')
#    plt.ylim(0,0.005)
#    plt.xlim(-0.8,0.8)
#    plt.tick_params(labelright=True,right=True)
#    plt.savefig('e_2 vs Split.pdf')
#        
    ######################################################################################################
    
    ##### Finding the Eccentricity ####
    
    
    x_steps = int(np.round(b_steps/2))
    y_steps = int(np.round(b_steps/2))
    
    # split into 2 becasue the for oblate and prolate geometries the e_2 is neg and pos
    
    
    ##### lower
    split_diff_y_min =   abs(abs(split_L[y_steps:]) - exp_split)
    split_diff_x_min =   abs(split_L[0:x_steps] - exp_split)
    
    ##### upper
    split_diff_y_max =   abs(abs(split_U[y_steps:]) - exp_split)
    split_diff_x_max =   abs(split_U[0:x_steps] - exp_split)
    
    ##### Eccentricity Indices #####
    
    ##### lower
    split_diff_y_min_index = int(y_steps + (np.argmin(split_diff_y_min)))   
    split_diff_x_min_index = np.argmin(split_diff_x_min)
    
    e_2_y_min = e_2_L[split_diff_y_min_index]
    e_2_x_min = e_2_L[split_diff_x_min_index]
    
    ##### upper
    split_diff_y_max_index = int(y_steps + (np.argmin(split_diff_y_max)))
    split_diff_x_max_index = np.argmin(split_diff_x_max)
  
    e_2_y_max = e_2_U[split_diff_y_max_index]
    e_2_x_max = e_2_U[split_diff_x_max_index]
    
    ##### min and max eccentriciteies ####
    
    e_2_range = [e_2_y_max,e_2_x_max,e_2_y_min,e_2_x_min]
    
    e_2_min = np.min(e_2_range)
    e_2_max = np.max(e_2_range)
    
    
    
    ##### Final e_2 ####
    
    e_2_final = e_2_min+((e_2_max - e_2_min)/2)
    error = ((e_2_max - e_2_min)/2)/e_2_final


    return [e_2_final,error,modes_L,modes_U]
    
#%%

#%%

def TM_P0_Pl(a,b,l,n_cell,n_bead):
    
    import numpy as np
    
    #################################################################################################################
        
    #def Wavelength_matrix(a,b,l,n_cell,n_bead,p):
#    a = 8*10**(-6)
#    b = np.linspace(5,13,4801)*10**(-6)
#    l = 133
#    n_cell = 1.34
#    n_bead = 1.67
    
    b_T = b
    l_T = l
    
    m = n_bead / n_cell
    z = -2.338107410459774 # We assume we work in the first radial mode q=1
    u = ((a**2) - (b_T**2))/3
    #    p = 0 # Because l=m asume the equitorial mode this would be represented by the redest most peaks of a spectrum.
    
    #################################################################################################################
    p = 0
    
    chi = 1/(m**2)
    
    #    lambda_te =  (m*2*np.pi*a)/ \
    lambda_tm =  (n_bead*2*np.pi*a)/ \
        ( \
        (l_T) - \
        ((z * (l_T**(1/3))) / ((2**(1/3)))) + \
        (((2 * p *(a - b_T)) + a) /(2 * b_T)) - \
        ((chi *m)/ (np.sqrt((m**2) - 1))) +\
        ((3 * (l_T**(-1/3)) * (z**2)) / (20* (2**(-1/3)))) -\
        ((z *(l_T**(-2/3))*(2 * p * ((a**3) - (b_T**3)) + (a**3))) / (12 * (b_T**(3)) * (2**(-2/3)))) - \
        ((z * (l_T**(-2/3)) * 2 * (m**3) * chi * ((2 * (chi**2)) - 3 )) / (12 * (2**(-2/3)) * (((m**2) - 1)**(3/2)))) + \
        ((((z**3)+10)*(l_T**(-1.)))/(1400*(2**(-1.)))) + \
        ((((2*p + 1)**2)  * (a**2) * ((b_T**2) * (1+(3*u)) - (a**2)))* (l_T**(-1.)) / (32 * (b_T**4) * (2**(-1.)))) \
        )
        
        
        
    p = l_T
    
    #    lambda_te =  (m*2*np.pi*a)/ \
    lambda_tm_l =  (n_bead*2*np.pi*a)/ \
        ( \
        (l_T) - \
        ((z * (l_T**(1/3))) / ((2**(1/3)))) + \
        (((2 * p *(a - b_T)) + a) /(2 * b_T)) - \
        ((chi *m)/ (np.sqrt((m**2) - 1))) +\
        ((3 * (l_T**(-1/3)) * (z**2)) / (20* (2**(-1/3)))) -\
        ((z *(l_T**(-2/3))*(2 * p * ((a**3) - (b_T**3)) + (a**3))) / (12 * (b_T**(3)) * (2**(-2/3)))) - \
        ((z * (l_T**(-2/3)) * 2 * (m**3) * chi * ((2 * (chi**2)) - 3 )) / (12 * (2**(-2/3)) * (((m**2) - 1)**(3/2)))) + \
        ((((z**3)+10)*(l_T**(-1.)))/(1400*(2**(-1.)))) + \
        ((((2*p + 1)**2)  * (a**2) * ((b_T**2) * (1+(3*u)) - (a**2)))* (l_T**(-1.)) / (32 * (b_T**4) * (2**(-1.)))) \
        )
        
        
    
    #################################################################################################################
    split = np.zeros((len(b)))
    summ = np.zeros((len(b)))
    e_2 = np.zeros((len(b)))
    geometry = np.zeros((len(b)))
    for i in range(len(b)):
        
        if a < b[i]:
            # positive split for prolate
            geometry[i] = -1
            e_2[i] = 1-((a**2)/(b[i]**2))
#            eccentricity[i] = 1-((a**2)/(b[i]**2))
            split[i] = lambda_tm[i] - lambda_tm_l[i]
            summ[i] = lambda_tm[i] + lambda_tm_l[i]
            
        elif a > b[i]:
            # negative split for oblate
            geometry[i] = 1
            e_2[i] = 1-((b[i]**2)/(a**2))
#            eccentricity[i] = 1-((b[i]**2)/(a**2))
            split[i] = lambda_tm[i] - lambda_tm_l[i]
            summ[i] = lambda_tm[i] + lambda_tm_l[i]
        else:
            # zero split for sphere
            geometry[i] = 0
            e_2[i] = 0
            split[i] = lambda_tm[i] - lambda_tm_l[i]
            summ[i] = lambda_tm[i] + lambda_tm_l[i]
    result = [split,e_2,summ,geometry,lambda_tm,lambda_tm_l]
    
    
    return result

#%%
#%%

def TE_P0_Pl(a,b,l,n_cell,n_bead):
    
    import numpy as np
    
    #################################################################################################################
        
    #def Wavelength_matrix(a,b,l,n_cell,n_bead,p):
#    a = 8*10**(-6)
#    b = np.linspace(5,13,4801)*10**(-6)
#    l = 133
#    n_cell = 1.34
#    n_bead = 1.67
#    
    b_T = b
    l_T = l
    
    m = n_bead / n_cell
    z = -2.338107410459774 # We assume we work in the first radial mode q=1
    u = ((a**2) - (b_T**2))/3
    #    p = 0 # Because l=m asume the equitorial mode this would be represented by the redest most peaks of a spectrum.
    
    #################################################################################################################
    p = 0
    
    chi= 1
    
    #    lambda_te =  (m*2*np.pi*a)/ \
    lambda_te =  (n_bead*2*np.pi*a)/ \
        ( \
        (l_T) - \
        ((z * (l_T**(1/3))) / ((2**(1/3)))) + \
        (((2 * p *(a - b_T)) + a) /(2 * b_T)) - \
        ((chi *m)/ (np.sqrt((m**2) - 1))) +\
        ((3 * (l_T**(-1/3)) * (z**2)) / (20* (2**(-1/3)))) -\
        ((z *(l_T**(-2/3))*(2 * p * ((a**3) - (b_T**3)) + (a**3))) / (12 * (b_T**(3)) * (2**(-2/3)))) - \
        ((z * (l_T**(-2/3)) * 2 * (m**3) * chi * ((2 * (chi**2)) - 3 )) / (12 * (2**(-2/3)) * (((m**2) - 1)**(3/2)))) + \
        ((((z**3)+10)*(l_T**(-1.)))/(1400*(2**(-1.)))) + \
        ((((2*p + 1)**2)  * (a**2) * ((b_T**2) * (1+(3*u)) - (a**2)))* (l_T**(-1.)) / (32 * (b_T**4) * (2**(-1.)))) \
        )
        
        
        
    p = l_T
    
    #    lambda_te =  (m*2*np.pi*a)/ \
    lambda_te_l =  (n_bead*2*np.pi*a)/ \
        ( \
        (l_T) - \
        ((z * (l_T**(1/3))) / ((2**(1/3)))) + \
        (((2 * p *(a - b_T)) + a) /(2 * b_T)) - \
        ((chi *m)/ (np.sqrt((m**2) - 1))) +\
        ((3 * (l_T**(-1/3)) * (z**2)) / (20* (2**(-1/3)))) -\
        ((z *(l_T**(-2/3))*(2 * p * ((a**3) - (b_T**3)) + (a**3))) / (12 * (b_T**(3)) * (2**(-2/3)))) - \
        ((z * (l_T**(-2/3)) * 2 * (m**3) * chi * ((2 * (chi**2)) - 3 )) / (12 * (2**(-2/3)) * (((m**2) - 1)**(3/2)))) + \
        ((((z**3)+10)*(l_T**(-1.)))/(1400*(2**(-1.)))) + \
        ((((2*p + 1)**2)  * (a**2) * ((b_T**2) * (1+(3*u)) - (a**2)))* (l_T**(-1.)) / (32 * (b_T**4) * (2**(-1.)))) \
        )
        
        
    
    #################################################################################################################
    split = np.zeros((len(b)))
    summ = np.zeros((len(b)))
    e_2 = np.zeros((len(b)))
    geometry = np.zeros((len(b)))
    for i in range(len(b)):
        
        if a < b[i]:
            # positive split for prolate
            geometry[i] = -1
            e_2[i] = 1-((a**2)/(b[i]**2))
#            eccentricity[i] = 1-((a**2)/(b[i]**2))
            split[i] = lambda_te[i] - lambda_te_l[i]
            summ[i] = lambda_te[i] + lambda_te_l[i]
            
        elif a > b[i]:
            # negative split for oblate
            geometry[i] = 1
            e_2[i] = 1-((b[i]**2)/(a**2))
#            eccentricity[i] = 1-((b[i]**2)/(a**2))
            split[i] = lambda_te[i] - lambda_te_l[i]
            summ[i] = lambda_te[i] + lambda_te_l[i]
        else:
            # zero split for sphere
            geometry[i] = 0
            e_2[i] = 0
            split[i] = lambda_te[i] - lambda_te_l[i]
            summ[i] = lambda_te[i] + lambda_te_l[i]
    result = [split,e_2,summ,geometry,lambda_te,lambda_te_l]
    
    
    
    
    return result

#%%

#%%
###############################################################################
def Mode_Refining(a_range,b_dev,l,n_cell_range,n_bead,exp_sum):

    
    
    from Functions import TE_MR 
    from Functions import TM_MR
    import numpy as np
    
    
    ##### mode refining #####
    
    b_range_L = np.array([a_range[0]-b_dev,a_range[0]+b_dev])
    b_range_U = np.array([a_range[1]-b_dev,a_range[1]+b_dev])
    ##### a_lower
    
    x = (len(b_range_L)*len(n_cell_range)*2)
    
    wave_combined_L = np.zeros((x,len(l)))
    constraints_L = []
    
    count = 0
    a = a_range[0]
    for j in range(len(b_range_L)):
        b = b_range_L[j]
        for k in range(len(n_cell_range)):
            n_cell = n_cell_range[k]
    
            count1 = count+1
            
            TE_result = TE_MR(a,b,l,n_cell,n_bead)
            TM_result = TM_MR(a,b,l,n_cell,n_bead)
            constraints_L.append([a,b,l,n_cell,n_bead,'TE'])
            constraints_L.append([a,b,l,n_cell,n_bead,'TM'])
            wave_combined_L[count,:] = TE_result[0]
            wave_combined_L[count+1,:] = TM_result[0]
    
            count = count1 + 1
    
    ##### a_upper
            
    x = (len(b_range_U)*len(n_cell_range)*2)
    
    wave_combined_U = np.zeros((x,len(l)))
    constraints_U = []
    
    count = 0
    a = a_range[1]
    for j in range(len(b_range_U)):
        b = b_range_U[j]
        for k in range(len(n_cell_range)):
            n_cell = n_cell_range[k]
    
            count1 = count+1
            
            TE_result = TE_MR(a,b,l,n_cell,n_bead)
            TM_result = TM_MR(a,b,l,n_cell,n_bead)
            constraints_U.append([a,b,l,n_cell,n_bead,'TE'])
            constraints_U.append([a,b,l,n_cell,n_bead,'TM'])
            wave_combined_U[count,:] = TE_result[0]
            wave_combined_U[count+1,:] = TM_result[0]
    
            count = count1 + 1
    
       
    ##### Finding the best fitting mode number #####
    
    #####a_lower
    modes_L = np.zeros(len(wave_combined_L[:,0])) 
    
    #####a_upper
    modes_U = np.zeros(len(wave_combined_L[:,0])) 
    
    for i in range(len(wave_combined_L[:,0])): 
        
        ######a_lower
        mode_sum_index_L = np.argmin(abs(wave_combined_L[i,:] - exp_sum))
        modes_L[i] = l[mode_sum_index_L]
    
        #####a_upper
        mode_sum_index_U = np.argmin(abs(wave_combined_U[i,:] - exp_sum))
        modes_U[i] = l[mode_sum_index_U]
    
    
    ##### Finding the mode ranges for a_upper and a_lower #####
    
    ##### a_lower
    mode_range_L = np.array([np.min(modes_L),np.max(modes_L)])    
    
    ##### a_upper 
    mode_range_U = np.array([np.min(modes_U),np.max(modes_U)])    
    
    
    modes_LU = [mode_range_L,mode_range_U]
    return [modes_LU]
    



#%%

#%%
def TE_MR(a,b,l,n_cell,n_bead):
    
    import numpy as np
    
    #################################################################################################################
        
    #def Wavelength_matrix(a,b,l,n_cell,n_bead,p):
#    a = 8*10**(-6)
#    b = np.linspace(5,13,4801)*10**(-6)
#    l = 133
#    n_cell = 1.34
#    n_bead = 1.67
#    
    b_T = b
    l_T = l
    
    m = n_bead / n_cell
    z = -2.338107410459774 # We assume we work in the first radial mode q=1
    u = ((a**2) - (b_T**2))/3
    #    p = 0 # Because l=m asume the equitorial mode this would be represented by the redest most peaks of a spectrum.
    
    #################################################################################################################
    p = 0
    
    chi= 1
    
    #    lambda_te =  (m*2*np.pi*a)/ \
    lambda_te =  (n_bead*2*np.pi*a)/ \
        ( \
        (l_T) - \
        ((z * (l_T**(1/3))) / ((2**(1/3)))) + \
        (((2 * p *(a - b_T)) + a) /(2 * b_T)) - \
        ((chi *m)/ (np.sqrt((m**2) - 1))) +\
        ((3 * (l_T**(-1/3)) * (z**2)) / (20* (2**(-1/3)))) -\
        ((z *(l_T**(-2/3))*(2 * p * ((a**3) - (b_T**3)) + (a**3))) / (12 * (b_T**(3)) * (2**(-2/3)))) - \
        ((z * (l_T**(-2/3)) * 2 * (m**3) * chi * ((2 * (chi**2)) - 3 )) / (12 * (2**(-2/3)) * (((m**2) - 1)**(3/2)))) + \
        ((((z**3)+10)*(l_T**(-1.)))/(1400*(2**(-1.)))) + \
        ((((2*p + 1)**2)  * (a**2) * ((b_T**2) * (1+(3*u)) - (a**2)))* (l_T**(-1.)) / (32 * (b_T**4) * (2**(-1.)))) \
        )
        
        
        
    p = l_T
    
    #    lambda_te =  (m*2*np.pi*a)/ \
    lambda_te_l =  (n_bead*2*np.pi*a)/ \
        ( \
        (l_T) - \
        ((z * (l_T**(1/3))) / ((2**(1/3)))) + \
        (((2 * p *(a - b_T)) + a) /(2 * b_T)) - \
        ((chi *m)/ (np.sqrt((m**2) - 1))) +\
        ((3 * (l_T**(-1/3)) * (z**2)) / (20* (2**(-1/3)))) -\
        ((z *(l_T**(-2/3))*(2 * p * ((a**3) - (b_T**3)) + (a**3))) / (12 * (b_T**(3)) * (2**(-2/3)))) - \
        ((z * (l_T**(-2/3)) * 2 * (m**3) * chi * ((2 * (chi**2)) - 3 )) / (12 * (2**(-2/3)) * (((m**2) - 1)**(3/2)))) + \
        ((((z**3)+10)*(l_T**(-1.)))/(1400*(2**(-1.)))) + \
        ((((2*p + 1)**2)  * (a**2) * ((b_T**2) * (1+(3*u)) - (a**2)))* (l_T**(-1.)) / (32 * (b_T**4) * (2**(-1.)))) \
        )
        
        
    
    #################################################################################################################
    Wave_Combined = np.zeros((len(l)))
    Wave_Split = np.zeros((len(l)))
#    eccentricity = np.zeros((len(b)))
#    geometry = np.zeros((len(b)))
    for i in range(len(l)):
        
        if a < b:
            
#            geometry[i] = -1
#            eccentricity[i] = np.sqrt(1-((a**2)/(b[i]**2)))
            Wave_Combined[i] = lambda_te[i] + lambda_te_l[i]
            Wave_Split[i] = abs(lambda_te[i] - lambda_te_l[i])
            
        elif a > b:
            
#            geometry[i] = 1
#            eccentricity[i] = np.sqrt(1-((b[i]**2)/(a**2)))
            Wave_Combined[i] = lambda_te[i] + lambda_te_l[i]
            Wave_Split[i] = abs(lambda_te[i] - lambda_te_l[i])

        else:
#            geometry[i] = 0
#            eccentricity[i] = 0
            Wave_Combined[i] = lambda_te[i] + lambda_te_l[i]
            Wave_Split[i] = abs(lambda_te[i] - lambda_te_l[i])

    result = [Wave_Combined,Wave_Split]
    
    return result

#%%
#%%
def TM_MR(a,b,l,n_cell,n_bead):
    
    import numpy as np
    
    #################################################################################################################
        
    #def Wavelength_matrix(a,b,l,n_cell,n_bead,p):
#    a = 8*10**(-6)
#    b = np.linspace(5,13,4801)*10**(-6)
#    l = 133
#    n_cell = 1.34
#    n_bead = 1.67
    
    b_T = b
    l_T = l
    
    m = n_bead / n_cell
    z = -2.338107410459774 # We assume we work in the first radial mode q=1
    u = ((a**2) - (b_T**2))/3
    #    p = 0 # Because l=m asume the equitorial mode this would be represented by the redest most peaks of a spectrum.
    
    #################################################################################################################
    p = 0
    
    chi = 1/(m**2)
    
    #    lambda_te =  (m*2*np.pi*a)/ \
    lambda_tm =  (n_bead*2*np.pi*a)/ \
        ( \
        (l_T) - \
        ((z * (l_T**(1/3))) / ((2**(1/3)))) + \
        (((2 * p *(a - b_T)) + a) /(2 * b_T)) - \
        ((chi *m)/ (np.sqrt((m**2) - 1))) +\
        ((3 * (l_T**(-1/3)) * (z**2)) / (20* (2**(-1/3)))) -\
        ((z *(l_T**(-2/3))*(2 * p * ((a**3) - (b_T**3)) + (a**3))) / (12 * (b_T**(3)) * (2**(-2/3)))) - \
        ((z * (l_T**(-2/3)) * 2 * (m**3) * chi * ((2 * (chi**2)) - 3 )) / (12 * (2**(-2/3)) * (((m**2) - 1)**(3/2)))) + \
        ((((z**3)+10)*(l_T**(-1.)))/(1400*(2**(-1.)))) + \
        ((((2*p + 1)**2)  * (a**2) * ((b_T**2) * (1+(3*u)) - (a**2)))* (l_T**(-1.)) / (32 * (b_T**4) * (2**(-1.)))) \
        )
        
        
        
    p = l_T
    
    #    lambda_te =  (m*2*np.pi*a)/ \
    lambda_tm_l =  (n_bead*2*np.pi*a)/ \
        ( \
        (l_T) - \
        ((z * (l_T**(1/3))) / ((2**(1/3)))) + \
        (((2 * p *(a - b_T)) + a) /(2 * b_T)) - \
        ((chi *m)/ (np.sqrt((m**2) - 1))) +\
        ((3 * (l_T**(-1/3)) * (z**2)) / (20* (2**(-1/3)))) -\
        ((z *(l_T**(-2/3))*(2 * p * ((a**3) - (b_T**3)) + (a**3))) / (12 * (b_T**(3)) * (2**(-2/3)))) - \
        ((z * (l_T**(-2/3)) * 2 * (m**3) * chi * ((2 * (chi**2)) - 3 )) / (12 * (2**(-2/3)) * (((m**2) - 1)**(3/2)))) + \
        ((((z**3)+10)*(l_T**(-1.)))/(1400*(2**(-1.)))) + \
        ((((2*p + 1)**2)  * (a**2) * ((b_T**2) * (1+(3*u)) - (a**2)))* (l_T**(-1.)) / (32 * (b_T**4) * (2**(-1.)))) \
        )
        
        
    
    #################################################################################################################
    Wave_Combined = np.zeros((len(l)))
    Wave_Split = np.zeros((len(l)))
#    eccentricity = np.zeros((len(b)))
#    geometry = np.zeros((len(b)))
    for i in range(len(l)):
        
        if a < b:
            
#            geometry[i] = -1
#            eccentricity[i] = np.sqrt(1-((a**2)/(b[i]**2)))
            Wave_Combined[i] = lambda_tm[i] + lambda_tm_l[i]
            Wave_Split[i] = abs(lambda_tm[i] - lambda_tm_l[i])
            
        elif a > b:
            
#            geometry[i] = 1
#            eccentricity[i] = np.sqrt(1-((b[i]**2)/(a**2)))
            Wave_Combined[i] = lambda_tm[i] + lambda_tm_l[i]
            Wave_Split[i] = abs(lambda_tm[i] - lambda_tm_l[i])

        else:
#            geometry[i] = 0
#            eccentricity[i] = 0
            Wave_Combined[i] = lambda_tm[i] + lambda_tm_l[i]
            Wave_Split[i] = abs(lambda_tm[i] - lambda_tm_l[i])

    result = [Wave_Combined,Wave_Split]
    
    return result

#%%