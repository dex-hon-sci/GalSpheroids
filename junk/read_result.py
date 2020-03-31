a#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 15:43:15 2019

@author: dexter
"""

import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from numpy.linalg import inv


###########initial parameters###############
c, H0=    299792.458, 68.0        #speed of light in km/s, Hubble constant in(km/s)/Mpc
M_sun = 4.53
############################################

#%%
def pd_read(filename,check_equvi):
    ####################################################################################################
    ## Reading profiler (by Bogdan) output log file by pandas
    ## It generate a python list of the parameters of  
    ## 
    ####################################################################################################
    ## filename: The ASCII log file from profiler
    ####################################################################################################  
    A = pd.read_table(filename,sep='\n')  #seperating the string by \n
    #print(filename)
    Gal_list=[]
    Gal_list.append(filename[2:9])  #recording the name of the galaxy
    #equvi = str(check_equvi)
    i=0 
    for i in range(A.shape[0]):
        if A.iloc[i][0] == "-----------------------------Detailed fit report:-------------------------------- ":
                ## record the rms for the residual as the second element of the list
                Gal_list.append(float(A.iloc[i-1][0].split()[3]))

        elif A.iloc[i][0][0:17] == "Sersic component:":
                ########################################################################################
                ## record the fitting parameters for the Sersic function 
                ## ['Sersic', np.array([mu_e. r_e, n])], np.array([mag]])
                ########################################################################################

                mu_e = float(A.iloc[i+1][0].split()[2])
                r_e = float(A.iloc[i+2][0].split()[2])
                n = float(A.iloc[i+3][0].split()[2])
                Sersic_mag = float(A.iloc[i+4][0].split()[4]) if check_equvi else 0
                
                Gal_list.append('Sersic')
                Gal_list.append(np.array([mu_e,r_e,n]))
                Gal_list.append(Sersic_mag)

        elif A.iloc[i][0][0:22] == "Core-Sersic component:":
                ########################################################################################
                ## record the fitting parameters for the Core Sersic function 
                ## ['CoreSersic', np.array([mu_p. r_e, r_b, n, alpha, gamma])], np.array([mag]])
                ########################################################################################
                mu_p = float(A.iloc[i+1][0].split()[2])
                r_e = float(A.iloc[i+2][0].split()[2])
                r_b = float(A.iloc[i+3][0].split()[2])
                n = float(A.iloc[i+4][0].split()[2])
                alpha = float(A.iloc[i+5][0].split()[2])
                gamma = float(A.iloc[i+6][0].split()[2])
                CoreSersic_mag = float(A.iloc[i+7][0].split()[4]) if check_equvi else 0

                
                Gal_list.append('CoreSersic')
                Gal_list.append(np.array([mu_p,r_e,r_b,n,alpha,gamma]))
                Gal_list.append(CoreSersic_mag)
                
        elif A.iloc[i][0][0:22] == "Exponential component:":
                ########################################################################################
                ## record the fitting parameters for the Exponential function 
                ## ['Exp', np.array([mu_0. h])], np.array([mag]])
                ########################################################################################
                mu_0 = float(A.iloc[i+1][0].split()[2])
                h = float(A.iloc[i+2][0].split()[2])
                Exp_mag = float(A.iloc[i+3][0].split()[4]) if check_equvi else 0
            
                Gal_list.append('Exp')                
                Gal_list.append(np.array([mu_0,h]))
                Gal_list.append(Exp_mag)


        elif A.iloc[i][0][0:25] == "Truncated disc component:":
                ########################################################################################
                ## record the fitting parameters for the Borken Exponential function 
                ## ['BrokenExp', np.array([mu_0td. r_b , h1, h2])], np.array([BrokenExp_mag])
                ########################################################################################
                mu_0td = float(A.iloc[i+1][0].split()[2])
                r_b = float(A.iloc[i+2][0].split()[2])
                h1 = float(A.iloc[i+3][0].split()[2])
                h2 = float(A.iloc[i+4][0].split()[2])
                BrokenExp_mag = float(A.iloc[i+5][0].split()[4]) if check_equvi else 0

                Gal_list.append('BrokenExp')                                
                Gal_list.append(np.array([mu_0td,r_b,h1,h2]))
                Gal_list.append(BrokenExp_mag)

        elif A.iloc[i][0][0:17] == "Ferrer component:":
                ########################################################################################
                ## record the fitting parameters for the Ferrer function 
                ## ['Ferrer', np.array([mu_0,r_0,alpha,beta])], np.array([Ferrer_mag]])
                ########################################################################################           
                mu_0 = float(A.iloc[i+1][0].split()[2])
                r_0 = float(A.iloc[i+2][0].split()[2])
                alpha = float(A.iloc[i+3][0].split()[2])
                beta = float(A.iloc[i+4][0].split()[2])
                Ferrer_mag = float(A.iloc[i+6][0].split()[4]) if check_equvi else 0
                
                Gal_list.append('Ferrer')                                
                Gal_list.append(np.array([mu_0,r_0,alpha,beta]))
                Gal_list.append(Ferrer_mag)
                
        elif A.iloc[i][0][0:24] == "Inclined disc component:":
                ########################################################################################
                ## record the fitting parameters for the Gaussian function 
                ## ['InclExp', np.array([mu_0z. z0])], np.array([InclExp_mag]])
                ########################################################################################
                mu_0z = float(A.iloc[i+1][0].split()[2])
                z_0 = float(A.iloc[i+2][0].split()[2])
                InclExp_mag = float(A.iloc[i+3][0].split()[4]) if check_equvi else 0
                
                Gal_list.append('InclExp')                                
                Gal_list.append(np.array([mu_0z,z_0]))
                Gal_list.append(InclExp_mag)
                   
        elif A.iloc[i][0][0:19] == "Gaussian component:":
                ########################################################################################
                ## record the fitting parameters for the Gaussian function 
                ## ['Gauss', np.array([mu_0. r_0, FWHM])], np.array([Gauss_mag]])
                ########################################################################################
                mu_0 = float(A.iloc[i+1][0].split()[2])
                r_0 = float(A.iloc[i+2][0].split()[2])
                FWHM = float(A.iloc[i+3][0].split()[2])          
                Gauss_mag = float(A.iloc[i+4][0].split()[4]) if check_equvi else 0
                
                Gal_list.append('Gauss')                                
                Gal_list.append(np.array([mu_0,r_0,FWHM]))
                Gal_list.append(Gauss_mag)
                
        elif A.iloc[i][0][0:24] == "PSF [nuclear] component:":
                ########################################################################################
                ## record the fitting parameters for the Gaussian function 
                ## ['PSF', np.array([mu_0])], np.array([PSF_mag]])
                ########################################################################################
                mu_0 = float(A.iloc[i+1][0].split()[2])
                PSF_mag = float(A.iloc[i+2][0].split()[4]) if check_equvi else 0
                
                Gal_list.append('PSF')                                                
                Gal_list.append(np.array([mu_0]))
                Gal_list.append(PSF_mag)
        elif A.iloc[i][0][0:23] == "Total galaxy magnitude:":
                ########################################################################################
                ## record the total magnitude of the galaxy 
                ## ['Total_mag',  total_mag]
                ########################################################################################
                total_mag = float(A.iloc[i][0].split()[3])
                
                Gal_list.append('Total_mag')                                                
                Gal_list.append(total_mag)            
    return Gal_list

#%%    
def read_list(input_file):
    with open(input_file, 'rb') as f:
        mylist = pickle.load(f)
    return mylist

#%%
def run_list(input_list,output_name,check_equvi):
    ####################################################################################################
    ## Put the log file information into a python readable format
    ## All the list of individual galaxy will be bundled into one python list 
    ## 
    ####################################################################################################
    ## input_list: A ASCII file containing the profiler log at each row 
    ## output_name: The filename for saving the list bundle
    ####################################################################################################  
    input_gal = pd.read_table(input_list,sep='\n')
    Gal_bundle = []
    j=0 
    for j in range(input_gal.shape[0]):
        Gal_list = pd_read(input_gal.iloc[j][0],check_equvi) #, check_equvi)
        Gal_bundle.append(Gal_list)
        #print(Gal_list)
        
    with open(output_name, 'wb') as f:
        pickle.dump(Gal_bundle, f)
        
    return Gal_bundle


#%%
def compare_para(list_input1,list_input2, func_index, element_index, para_name):  
    ####################################################################################################
    ## Compare the parameters of the fitting components
    ## Export an ASCII file containing the  
    ####################################################################################################
    ## list_input1, 2: The python lists containing the information of ALL galaxy
    ## func_index: The function you are inteested in 
    ## element_index: The parameters of interest in the function 
    ## para_name: the parameter name
    ####################################################################################################
    fig, ax = plt.subplots()
    
    list1, list2 = read_list(list_input1), read_list(list_input2)
    delta_para, index, gal_name = [], [], []
    i = 0 
    for i in range(len(list1)):
        delta_para.append(list1[i][func_index][element_index] - list2[i][func_index][element_index])
        index.append(i)
        gal_name.append(list1[i][0])
    #delta_para = np.array(delta_para)

    ax.plot(index, delta_para, 'bo')
    
    j=0
    for j in range(len(list1)):
       ax.text(index[j],delta_para[j], gal_name[j],fontsize=12)
        
    #plt.xlabel(para_name,fontsize=16)
    plt.ylabel("$\Delta$ %s" %para_name,fontsize=16)
    plt.hlines(1.0, 0, 150, linestyle="dashed",linewidth=3, color='b' )
    plt.hlines(-1.0, 0, 150, linestyle="dashed",linewidth=3, color='b' )
    return delta_para

# new criteria 
# look at the number of Sersic/Ferrer/Exp profile -> check which one are bulge/extended disk-> select
# what to compare major vs equvi
# Bulge, extended Disk, nuclear Disk
# Bulge = Sersic, with the higher sersic index 
## bulge_can = [array1,array2] (if  label = "sersic" or label = "CoreSersic"), bulge = array1 (if array1_n = max(array_n))
    # maybe Sph not bulge
# nuclear Disk = exp that has lowest h, usually h<7
    ## nucDisk_can = [array1,array2] (if label = "Exp", h<7 ), nucDisk = 
# intermediate Disk = exp lower h than the extended disk
    ## intDisk_can = array1 (if label = "Exp", "incl_Exp", "Broken_Exp", h<7 ) if at R_tot/2, Sersic(x) > Exp (x)
# Extended Disk = exp with the largest h 
    ##  extDisk = array1 (if label = "Exp", "incl_Exp", "Broken_Exp", h<7 ) if at R_tot/2, Sersic(x) < Exp (x) 
## primary Bar = bridge between disk and bulge; secondary Bar = internal bar form out of instability    
    ##  Bar_can = [array1, array2] (if label = "Ferrer" h<7 ), prim_bar = array1, sec_bar = array2 if r_out1 > r_out2
# Anse = the accumulation of stars builded by 
    ## Anse_can = array1 (if label = "Gauss", )
# rings/ spiral arm = 
    ## Ring_can = [array1] (if label = "Gauss", on the disk)
# 
    
# For Bars and other stuff,create a checklist of the components markers [NGC1234, Bar = 2]
# Loop through the marker table and cross matches the feature to the parameters
# override list [NGC1234, Sersic == disk], check override list first then 
    #output cpt label replace func label
    
#cpt_label = [Spheroid, NucDisk, IntDisk , ExtDisk, PrimBar, SecBar] 
    
    
def compare_para2(list_input1,list_input2, feature, element_index, para_name):  
    ####################################################################################################
    ## Compare the parameters of the fitting components
    ## Export an ASCII file containing the  
    ## 
    ####################################################################################################
    ## list_input1, 2: The python lists containing the information of ALL galaxy
    ## func_index: The function you are inteested in 
    ## feature: The component you are interested in, acceptable input: 
    ## 1) Bulge, 2)Core-depleted Bulge 3) Extended Disk, 4) Nuclear Disk, 5) Intermediate Disk, 
    ## 6) outer dearth (typeII truncation), 7) outer Surplus (typeIII truncation) 
    ## 8) Primary Bar, 9) Secondary Bar 
    ## Extra label to override the selection scheme
    #################################################################################################### 
    fig, ax = plt.subplots()
    
    list1, list2 = read_list(list_input1), read_list(list_input2)
    delta_para, index, gal_name = [], [], []
    
    #i = 0 
    for i in range(len(list1)):  #i=row
        j=0
        for j in range(len(list1[i])): #j=index and element
            
            if list1[i][j] == feature and list2[i][j] == feature:
                delta_para.append(list1[i][j+1][element_index] - list2[i][j+1][element_index])
                index.append(i)
                gal_name.append(list1[i][0])
            else:
                pass
    #delta_para = np.array(delta_para)

    #ax.plot(index, delta_para, 'bo')
    ax.plot(index, delta_para, 'bo')

    k=0
    for k in range(len(index)):
       ax.text(index[k],delta_para[k], gal_name[k],fontsize=12)
       
    print(feature, len(index))
        
    #plt.xlabel(para_name,fontsize=16)
    plt.ylabel("$\Delta$ %s" %para_name,fontsize=16)
    plt.hlines(1.0, 0, 150, linestyle="dashed",linewidth=3, color='b')
    plt.hlines(-1.0, 0, 150, linestyle="dashed",linewidth=3, color='b')
    return delta_para


def compare_para3(list_input1,list_input2, feature, element_index, para_name1,para_name2):  
    ####################################################################################################
    ## Compare the parameters of the fitting components
    ## Export an ASCII file containing the  
    ## 
    ####################################################################################################
    ## list_input1, 2: The python lists containing the information of ALL galaxy
    ## func_index: The function you are inteested in 
    ## feature: The component you are interested in, acceptable input: 
    ## 1) Bulge, 2)Core-depleted Bulge 3) Extended Disk, 4) Nuclear Disk, 5) Intermediate Disk, 
    ## 6) outer dearth (typeII truncation), 7) outer Surplus (typeIII truncation) 
    ## 8) Primary Bar, 9) Secondary Bar 
    ## Extra label to override the selection scheme
    #################################################################################################### 
    fig, ax = plt.subplots()
    
    list1, list2 = read_list(list_input1), read_list(list_input2)
    delta_para, index, gal_name = [], [], []
    
    #i = 0 
    for i in range(len(list1)):  #i=row
        j=0
        for j in range(len(list1[i])): #j=index and element
            
            if list1[i][j] == feature and list2[i][j] == feature:
                delta_para.append(list1[i][j+1][element_index]/list2[i][j+1][element_index])
                index.append(i)
                gal_name.append(list1[i][0])
            else:
                pass
    #delta_para = np.array(delta_para)

    #ax.plot(index, delta_para, 'bo')
    ax.plot(index, delta_para, 'bo')

    #k=0
    #for k in range(len(index)):
    #   ax.text(index[k],delta_para[k], gal_name[k],fontsize=12)
       
    print(feature, len(index))
        
    plt.xlabel(feature,fontsize=16)
    plt.ylabel("$log(%s/%s)$" %(para_name1,para_name2),fontsize=16)
    plt.yscale( 'log' )
    plt.xticks([1,50,100], [])
    
    return delta_para


#%%
def compare_mag(list_input1,list_input2, element_index, para_name):  
    ####################################################################################################
    ## Compare the magnitude of the fitting components
    ## 
    ####################################################################################################
    ## list_input1, 2: The python lists containing the information of ALL galaxy
    ## func_index: The function you are inteested in 
    ## element_index: The parameters of interest in the function 
    ## para_name: the parameter name
    #################################################################################################### 
    list1, list2 = read_list(list_input1), read_list(list_input2)
    delta_para, index, gal_name = [], [], []
    i = 0 
    for i in range(len(list1)):
        delta_para.append(list1[i][element_index] - list2[i][element_index])
        index.append(i)
        gal_name.append(list1[i][0])

    ax.plot(index, delta_para, 'bo')
    
    j=0
    for j in range(len(list1)):
       ax.text(index[j],delta_para[j], gal_name[j],fontsize=12)
        
    #plt.xlabel(para_name,fontsize=16)


    plt.ylabel("$\Delta$ %s" %para_name,fontsize=16)
    return delta_para
#%%
def cpt_seperator_demo(input_list_name): #, override_instruction):
    ####################################################################################################
    ## replacing analytic function label to component label
    ## 
    ####################################################################################################
    ## input_list_name:
    ## output_list_name:
    ## override_instruction
    #################################################################################################### 
    sample = read_list(input_list_name) #sample package
    
    master_Bulge_can, master_Bulge_can_index = [], []
    master_CoreBulge_can, master_CoreBulge_can_index = [], []

    master_Disk_can1, master_Disk_can2, master_Disk_can3 = [],[],[]
    master_Disk_can1_index, master_Disk_can2_index, master_Disk_can3_index = [],[],[]
    master_Bar_can, master_Bar_can_index = [], []
    master_Ring_can, master_Ring_can_index = [],[]
    
    master_Bulge_can_row = []
    master_CoreBulge_can_row = []
    master_Disk_can1_row = []
    master_Disk_can2_row = []
    master_Disk_can3_row = []
    master_Bar_can_row = []
    master_Ring_can_row = []
    
    for number in range(len(sample)):
        
        #print(sample[number])
        sample_indi = sample[number]    #indivdual sample
        
        Bulge_can,Bulge_can_index = [], []
        CoreBulge_can,CoreBulge_can_index = [],[]
        Disk_can1, Disk_can2, Disk_can3, Disk_can1_index, Disk_can2_index, Disk_can3_index = [],[],[],[],[],[]
        Bar_can, Bar_can_index = [], []
        Ring_can, Ring_can_index = [],[]
        
        Bulge_can_row = []
        CoreBulge_can_row = []

        Disk_can1_row = []
        Disk_can2_row = []
        Disk_can3_row = []
        Bar_can_row = []
        Ring_can_row = []
        
        for index in range(len(sample_indi)):
            ###############################
            # Bulge and Lens

            if sample_indi[index] == "Sersic":
                Bulge_can.append(sample_indi[index+1])
                Bulge_can_index.append(index)
                Bulge_can_row.append(number)
                #print(sample_indi[0],"Sersic")
            ###############################
            # CoreBulge
            if sample_indi[index] == "CoreSersic":
                CoreBulge_can.append(sample_indi[index+1])
                CoreBulge_can_index.append(index)
                CoreBulge_can_row.append(number)
                #print(sample_indi[0],"CoreSersic")
            ###############################       
            # nucDisk, intDisk, extDisk            
            elif sample_indi[index] == "Exp":
                Disk_can1.append(sample_indi[index+1])
                Disk_can1_index.append(index)
                Disk_can1_row.append(number)
                #print(sample_indi[0],"Exp")

            elif sample_indi[index] == "BrokenExp":
                Disk_can2.append(sample_indi[index+1])
                Disk_can2_index.append(index)
                Disk_can2_row.append(number)
                #print(sample_indi[0],"BrokenExp")

            elif sample_indi[index] == "InclExp":
                Disk_can3.append(sample_indi[index+1])   
                Disk_can3_index.append(index)
                Disk_can3_row.append(number)
                #print(sample_indi[0],"InclExp")

            ################################
            # primBar, secBar
            elif sample_indi[index] == "Ferrer":
                Bar_can.append(sample_indi[index+1])
                Bar_can_index.append(index)
                Bar_can_row.append(number)
                #print(sample_indi[0],"Ferrer")

            ################################
            # Anse, Rings
            elif sample_indi[index] == "Gauss":
                Ring_can.append(sample_indi[index+1])
                Ring_can_index.append(index)  
                Ring_can_row.append(number)
                #print(sample_indi[0],"Gauss")

            ################################
        master_Bulge_can.append(Bulge_can), master_Bulge_can_index.append(Bulge_can_index) 
        master_CoreBulge_can.append(CoreBulge_can), master_CoreBulge_can_index.append(CoreBulge_can_index) 
        master_Disk_can1.append(Disk_can1), master_Disk_can1_index.append(Disk_can1_index)
        master_Disk_can2.append(Disk_can2), master_Disk_can2_index.append(Disk_can2_index)
        master_Disk_can3.append(Disk_can3), master_Disk_can3_index.append(Disk_can3_index)
        master_Bar_can.append(Bar_can), master_Bar_can_index.append(Bar_can_index)
        master_Ring_can.append(Ring_can), master_Ring_can_index.append(Ring_can_index)
        
        master_Bulge_can_row.append(Bulge_can_row)
        master_CoreBulge_can_row.append(CoreBulge_can_row)
        master_Disk_can1_row.append(Disk_can1_row)
        master_Disk_can2_row.append(Disk_can2_row)
        master_Disk_can3_row.append(Disk_can3_row)
        master_Bar_can_row.append(Bar_can_row)
        master_Ring_can_row.append(Ring_can_row)
        
    return{"master_Bulge_can":master_Bulge_can, "master_Bulge_can_index":master_Bulge_can_index, "master_Bulge_can_row":master_Bulge_can_row,\
           "master_CoreBulge_can":master_CoreBulge_can, "master_CoreBulge_can_index":master_CoreBulge_can_index, "master_CoreBulge_can_row":master_CoreBulge_can_row,\
           "master_Disk_can1":master_Disk_can1, "master_Disk_can1_index":master_Disk_can1_index, "master_Disk_can1_row":master_Disk_can1_row,\
           "master_Disk_can2":master_Disk_can2, "master_Disk_can2_index":master_Disk_can2_index,"master_Disk_can2_row":master_Disk_can2_row,\
           "master_Disk_can3":master_Disk_can3, "master_Disk_can3_index":master_Disk_can3_index,"master_Disk_can3_row":master_Disk_can3_row,\
           "master_Bar_can":master_Bar_can, "master_Bar_can_index":master_Bar_can_index,"master_Bar_can_row":master_Bar_can_row,\
           "master_Ring_can":master_Ring_can, "master_Ring_can_index":master_Ring_can_index, "master_Ring_can_row":master_Ring_can_row}
            
                    
#%%
def cpt_classifier_demo(input_list_name, input_sep_dict ,output_list_name, override_instruction):
    
    sample = read_list(input_list_name) #sample package
    sep_dict = input_sep_dict
    
    ###Bulge######
    for number in range(len(sep_dict["master_Bulge_can"])):
        
        Bulge_can = sep_dict["master_Bulge_can"][number]
        Bulge_can_index = sep_dict["master_Bulge_can_index"][number]# the cpt marker
        Bulge_can_row = sep_dict["master_Bulge_can_row"][number] #the galaxy
        
        #print("Gal",sample[number][0], number, "Sersic")

        if len(Bulge_can) > 1:
            n_list = [Bulge_can[i][2] for i in range(len(Bulge_can))] # Sersic index n
            
            lens_index = Bulge_can_index[n_list.index(min(n_list))]
            lens_row = Bulge_can_row[n_list.index(min(n_list))]
            
            bulge_index = Bulge_can_index[n_list.index(max(n_list))]
            bulge_row = Bulge_can_row[n_list.index(max(n_list))]
            
            sample[lens_row][lens_index] = "Lens"
            sample[bulge_row][bulge_index] = "Bulge"
                
            pass
        elif len(Bulge_can) == 1:
            row = Bulge_can_row[0]
            index = Bulge_can_index[0]
            
            sample[row][index]= "Bulge"
            pass
        else: 
            pass    
        
    #####CoreBulge###
    for number in range(len(sep_dict["master_CoreBulge_can"])):
        CoreBulge_can = sep_dict["master_CoreBulge_can"][number]
        CoreBulge_can_index = sep_dict["master_CoreBulge_can_index"][number]# the cpt marker
        CoreBulge_can_row = sep_dict["master_CoreBulge_can_row"][number] #the galaxy
        
        #print("Gal",sample[number][0], number,"CoreBulge")
        if CoreBulge_can ==[]:
            pass
        else:
            row = CoreBulge_can_row[0]
            index = CoreBulge_can_index[0]
            
            sample[row][index] = "CoreBulge"
            pass
    ####Disk1######
    for number in range(len(sep_dict["master_Disk_can1"])):
        Disk_can1 = sep_dict["master_Disk_can1"][number]
        Disk_can1_index = sep_dict["master_Disk_can1_index"][number]# the cpt marker
        Disk_can1_row = sep_dict["master_Disk_can1_row"][number] #the galaxy
        
       # print("Gal",sample[number][0], number,"Exp")
        if Disk_can1 ==[]:
            pass
        else:
            row = Disk_can1_row[0]
            index = Disk_can1_index[0]
            if (Disk_can1[0][1] < 4):  # h < 4
                sample[row][index] = "nucDisk"
                pass
            elif ((Disk_can1[0][1]) > 4):   # h > 4
                sample[row][index] = "Disk"
                pass
            else: 
                pass
            
    #####Disk2####
    for number in range(len(sep_dict["master_Disk_can2"])):
        Disk_can2 = sep_dict["master_Disk_can2"][number]
        Disk_can2_index = sep_dict["master_Disk_can2_index"][number]# the cpt marker
        Disk_can2_row = sep_dict["master_Disk_can2_row"][number] #the galaxy
                
        #print("Gal",sample[number][0], number,"BrokenExp")
        if Disk_can2 ==[]:
            pass
        else:
            row = Disk_can2_row[0]
            index = Disk_can2_index[0]
            
            if (Disk_can2[0][2] < 4):  # h < 4
                sample[row][index] = "nucDisk"
                pass
            elif ((Disk_can2[0][2]) > 4):   # h > 4
                sample[row][index] = "Disk"
                pass
            else: 
                pass    
    #####Disk3###
    for number in range(len(sep_dict["master_Disk_can3"])):
        Disk_can3 = sep_dict["master_Disk_can3"][number]
        Disk_can3_index = sep_dict["master_Disk_can3_index"][number]# the cpt marker
        Disk_can3_row = sep_dict["master_Disk_can3_row"][number] #the galaxy
        

        #print(Disk_can3,Disk_can3_index,Disk_can3_row)
        
        #print("Gal",sample[number][0], number,"InclExp")
        if Disk_can3 ==[]:
            pass
        else:

            row = Disk_can3_row[0]
            index = Disk_can3_index[0]
            
            sample[row][index] = "Disk"
            pass
 
    ####Bar######
    for number in range(len(sep_dict["master_Bar_can"])):
        Bar_can = sep_dict["master_Bar_can"][number]
        Bar_can_index = sep_dict["master_Bar_can_index"][number]# the cpt marker
        Bar_can_row = sep_dict["master_Bar_can_row"][number] #the galaxy
        
        #print("Gal",sample[number][0], number, "Ferrer")
        if Bar_can == []:
            pass
        else:
            if len(Bar_can) > 1:
                r_out_list = [Bar_can[i][2] for i in range(len(Bar_can))] # Sersic index n
                
                Bar1_index = Bar_can_index[r_out_list.index(max(r_out_list))]
                Bar1_row = Bar_can_row[r_out_list.index(max(r_out_list))]
            
                Bar2_index = Bar_can_index[r_out_list.index(min(r_out_list))]
                Bar2_row = Bar_can_row[r_out_list.index(min(r_out_list))]
            
                sample[Bar1_row][Bar1_index] = "PrimBar"
                sample[Bar2_row][Bar2_index] = "SecBar"
                
                pass
            elif len(Bar_can) == 1:
                
                row = Bar_can_row[0]
                index = Bar_can_index[0]
            
                sample[row][index]= "PrimBar"
                pass
            
            else:
                pass
            
    ####Ring#####
    for number in range(len(sep_dict["master_Ring_can"])):
        Ring_can = sep_dict["master_Ring_can"][number]
        Ring_can_index = sep_dict["master_Ring_can_index"][number]# the cpt marker
        Ring_can_row = sep_dict["master_Ring_can_row"][number] #the galaxy
        
        #print("Gal",sample[number][0], number,"Gauss")
        if Ring_can ==[]:
            pass
        
        else:
            if len(Ring_can) > 1:
                for i in range(len(Ring_can_row)):
                    
                    row, index = Ring_can_row[i], Ring_can_index[i] 
                    sample[row][index] = "Ring"
                    pass
            elif len(Ring_can) == 1:

                row = Ring_can_row[0]
                index = Ring_can_index[0]
            
                sample[row][index] = "Ring"
                pass
            else: 
                pass
    ###override##
    if len(override_instruction) == None or len(override_instruction) == []:
        pass
    else:
        for number in range(len(override_instruction)):   
            for i in range(len(sample)):
                if sample[i][0] == override_instruction[number]:
                    cpt_index = override_instruction[number+1]
                    sample[i][cpt_index] = override_instruction[number+2]
                    #print("replaced!", sample[i][cpt_index], override_instruction[number+2])
                else:
                    pass
    #############
    with open(output_list_name, 'wb') as f:
        pickle.dump(sample, f)
        
    return sample
#%%               
 ## plotting ###    
def plot_compare_rms(input_list1,input_list2):
    fig, ax = plt.subplots()
    # same size
    input_list1 = read_list(input_list1)
    input_list2 = read_list(input_list2)
    index, rms1, rms2, gal_name= [], [], [], []
    
    for i in range(len(input_list1)):        
        if input_list1[i][0] == input_list2[i][0]: #check name consistency
            index.append(i)
            rms1.append(input_list1[i][1])
            rms2.append(input_list2[i][1])
            gal_name.append(input_list1[i][0])
            
        
    plt.hlines(np.average(rms1), 0, 104, linestyle="solid",linewidth=3, color='#7e491b')
    plt.hlines(np.average(rms2), 0, 104, linestyle="solid",linewidth=3, color='#244b87')
               
    plt.hlines(np.average(rms1)+np.std(rms1), 0, 104, linestyle="dashed",linewidth=3, color='#7e491b')
    plt.hlines(np.average(rms2)+np.std(rms2), 0, 104, linestyle="dashed",linewidth=3, color='#244b87')   
               
    plt.hlines(np.average(rms1)-np.std(rms1), 0, 104, linestyle="dashed",linewidth=3, color='#7e491b')
    plt.hlines(np.average(rms2)-np.std(rms2), 0, 104, linestyle="dashed",linewidth=3, color='#244b87')  
    
    s_rms1, s_rms2 = np.std(rms1), np.std(rms2)
    avg_rms1, avg_rms2 = np.average(rms1), np.average(rms2)
    
    #print(avg_rms1, avg_rms2,s_rms1, s_rms2)
    
    x_edge1,y_edge1= [0,0,104,104], [avg_rms1-s_rms1,avg_rms1+s_rms1,avg_rms1+s_rms1,avg_rms1-s_rms1]
    x_edge2,y_edge2= [0,0,104,104], [avg_rms2-s_rms2,avg_rms2+s_rms2,avg_rms2+s_rms2,avg_rms2-s_rms2]
           
    ax.fill(x_edge1, y_edge1, alpha=0.2, color='#e87964')
    ax.fill(x_edge2, y_edge2, alpha=0.2, color='#41a9c4')
          
    ax.plot(index, rms1, 'o', color='#d5721d', label= "multi-cpt")
    ax.plot(index, rms2, 'o', color='#1d64d5', label="B+D")

    #j=0
    #for j in range(len(input_list1)):
    #   ax.text(index[j],rms1[j], gal_name[j], fontsize=12)
       
    ax.text(50, np.average(rms1)+0.005, r"$\langle \Delta rms \rangle$ = %s $\pm$ %s" %(round(avg_rms1,3),round(s_rms1,2)),fontsize=14, color="#453816")
    ax.text(50, np.average(rms2)+0.005, r"$\langle \Delta rms \rangle$ = %s $\pm$ %s"  %(round(avg_rms2,3),round(s_rms2,2)),fontsize=14, color="#162945")    

    #plt.xticks(index,gal_name,rotation='vertical', fontsize)
    plt.ylabel("$\Delta$ rms", fontsize=16)
    plt.legend()   
    plt.xticks([1,50,100], [])
    return fig, ax

#%%
#Ctable_Into_all.txt
#Dtable_Into_all_compare.txt
#Dtable_Into_BD_all_compare.txt


    
def plot_hist_percentage(input_list, dist):
    #category_names = ["Bulge", "CoreBulge", "PrimBar","SecBar","Lens","nucDisk","Disk","Ring"]   
    category = {
                "nucDisk" : "#165180",
                "SecBar" : "#a88e25",
                "Bulge" :"#ee2b10",
                "CoreBulge" : "#b62c19",
                "Lens": "#d67a6d",
                "IntDisk": "#6d88d6",
                "PrimBar": "#dad035",
                "Disk" : "#3535da", 
                "Ring" : "#12e3e7", 
                "Point Source": "#3cf017"}
                
    #Q = create_table("Dtable_Into_all_compare.txt", "Ctable_Into_all.txt", 4.53,103)
    #P = create_table("Dtable_Into_BD_all_compare.txt", "Ctable_Into_all.txt", 4.53,103)

    
    category_names = list(category.keys())
    catagory_colour = list(category.values())
        
    input_list = read_list(input_list)
    
    gal_dict, gal_cpt = {},{} #gal_dict with gal name as key, array of mag as data, gal_cpt same but with cpt
    
    for row in range(len(input_list)):
        label, data = [], []
        #zdist = Q[8][row]
        #D = Q[11][row]
        #print(input_list[row][0],zdist,D)
        D = dist[row]
        for item in range(len(category_names)):  #re-shuffle
            #print(input_list[row][0])
            for index in range(len(input_list[row])):
                if input_list[row][index] == category_names[item]:  #matching cpt
                    
                    #print(input_list[row][index])
                    m = input_list[row][index+2]
                    M = m - 25 -5*np.log10(D)

                    Lum = 10**((M-(4.53))/(-2.5))
                    #Lum = input_list[row][index+2]
                    label.append(category_names[item])
                    data.append(Lum)
                else:
                    pass
        gal_dict[input_list[row][0]], gal_cpt[input_list[row][0]] = data, label
        
    #print("gal_dict",gal_dict,"gal_cpt",gal_cpt)
    
    #shuffle-to-order-then-plot
    #################################################  
    #gal_labels = list(gal_dict.keys())
    #gal_data = np.array(list(gal_dict.values()))
    #gal_data_cum = gal_data.cumsum()
    
    gal_labels = list(gal_dict.keys())
    gal_data = np.array(list(gal_dict.values()))
    gal_data_cum = gal_data.cumsum()

    
    fig, ax = plt.subplots(figsize=(10, 20))
    
    ax.invert_yaxis()
    #ax.xaxis.set_visible(False)
    #ax.set_xlim(0, np.sum(gal_data, axis=1).max())
    
    for i in range(len(input_list)): # loop through the galaxy list
        gal_name = input_list[i][0]
        #print(gal_name)
        widths, starts = 0,0  
        
        gal_data_cum = np.array(gal_data[i]).cumsum()

        for j in range(len(gal_data[i])): # loop through the property & cpt list
            
            widths = gal_data[i][j]
            
            starts = gal_data_cum[j] - widths
            
            #print(gal_cpt[gal_name][j], category[gal_cpt[gal_name][j]],widths,starts)
            ax.barh(gal_labels[i], widths, left=starts, height=0.9, label=gal_cpt[gal_name][j], color=category[gal_cpt[gal_name][j]])

    #mass bin 
    #morph Bin
    #specific feature bin
    
    #ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),loc='lower left', fontsize='small')
    plt.xlabel( '$L_*/L_{\odot}$' ,fontsize=16)
    #plt.xscale( 'log' )
    return fig, ax

#%%
def lookup_bundle(gal_bundle,gal_name):
    # to lookup singular entry of the bundle
    GN= read_list(gal_bundle)
    for row in range(len(GN)):
        if GN[row][0]==gal_name:
            print(GN[row])
        else:
            pass
        
#%%
def read_table(table):
    D=np.genfromtxt(table, dtype='float')
    return D
#%%
def grab_dist(dir_dist_list,dist_list):
    dir_dist_list1 ,dir_dist_list2 = np.genfromtxt(dir_dist_list, dtype='str'), np.genfromtxt(dir_dist_list, dtype='float')
    
    dist_list1,dist_list2 = np.genfromtxt(dist_list, dtype='str'), np.genfromtxt(dist_list, dtype='float')
    
    ##calculate v/H = dist
    ## replace dir_dist and calculate scale
    ## store in numpy array
    
    name = dist_list1[:,0]
    vel , vel_err= dist_list2[:,3], dist_list2[:,4]
    scale = dist_list2[:,5]
    
    dist = vel / H0 
    dist_err = vel_err / H0
    
    #####
    dir_name = dir_dist_list1[:,0]
    dir_dist, dir_dist_err = dir_dist_list2[:,1], dir_dist_list2[:,2]
    dir_method = dir_dist_list1[:,5]
    dir_method_flag = dir_dist_list2[:,7]
    #####
    dir_scale = dir_dist* ((1e3) / 206264.806) #turn Mpc/rad to kpc/arcsec 
    
    ####replace####
    for row in range(len(name)):
        for row2 in range(len(dir_name)):
            if name[row] == dir_name[row2]:
                dist[row] = dir_dist[row2]
                scale[row] = dir_scale[row2]
            else:
                pass
                
    return {"Gal_name": name, "Dist": dist, "Scale": scale}

#%%
def grab_parameter(filename, keyword, number):
    table = read_list(filename)
    storage = []
    #generate a list of singular parameter
    for row in range(len(table)):
        for item in range(len(table[row])):
            if table[row][item] in (keyword):
                storage.append(table[row][item+1][number])
                
    storage = np.array(storage)       
    return storage

def grab_mag(filename, keyword):
    table = read_list(filename)
    storage = []
    #generate a list of singular parameter
    for row in range(len(table)):
        for item in range(len(table[row])):
            if table[row][item] in (keyword):
                storage.append(table[row][item+2])
                
    storage = np.array(storage)       
    return storage

#%%
def plot_sizeMass_sub(input_list):
    return None

def plot_compare_Sph_shift(input_list):
    return None
#%%
##execution area###        
    
override_list_maj = \
["NGC2862",2,"Disk","NGC2872",5,"Point Source", "NGC3805",5,"IntDisk","NGC3805",8,"Point Source","NGC3812",2,"Background",\
 "NGC3872",5,"Halo","NGC3940",5,"Point Source", "NGC4065",5,"Point Source", "NGC4555",5,"Point Source", "NGC4570", 2, "Disk",\
 "UGC8736",2,"Disk" ]

override_list_equ = \
["NGC2862",2,"Disk","NGC2872",5,"Point Source", "NGC3805",5,"IntDisk","NGC3805",8,"Point Source","NGC3812",2,"Background",\
 "NGC3872",5,"Halo","NGC3940",5,"Point Source", "NGC4065",5,"Point Source", "NGC4555",5,"Point Source",\
 "UGC8736",2,"Disk" ]
#run_list("equvi_list.txt", "Gal_bundle",True)
#A = cpt_seperator_demo("Gal_bundle")
#cpt_classifier_demo('Gal_bundle',A,'Gal_bundle_cpt',override_list)

name = grab_dist("./distance/dir_dist_list.txt","./distance/dist_list.txt")["Gal_name"]
DD = grab_dist("./distance/dir_dist_list.txt","./distance/dist_list.txt")["Dist"]
scale = grab_dist("./distance/dir_dist_list.txt","./distance/dist_list.txt")["Scale"]

R_e = grab_parameter("Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"], 1) #get Re

mag = grab_mag("Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"])


CC2 = np.genfromtxt("Ctable_Into_all.txt" , dtype='float')

AbsMag_g = CC2[:,6]
AbsMag_i = CC2[:,8]


from plot_execute import ML_relation_Iband
from plot_execute import cal_Mass
from plot_execute import create_table
from plot_execute import Mass_Re_plot
from plot_execute import plot_cut
Mass_bulge_Into = cal_Mass(mag,DD,ML_relation_Iband(AbsMag_g,AbsMag_i).Into13_MassRatio,M_sun)
Re_kpc = R_e * scale

Mag = mag - 25 -5*np.log10(DD)
Lum = 10**((Mag-(4.53))/(-2.5))

Ie = Lum/(np.pi*Re_kpc**2)

import matplotlib.gridspec as gridspec
import matplotlib.patches as patches


#plot_cut(S)
fig, ax = plt.subplots()

Mass_Re_plot(Mass_bulge_Into,Re_kpc,name,'bo','Bulge Mass estimated by us using stellar M/L ratio from Into et al.2013',1.0)
plt.show()


def fun_plane(Re,K):
    return K*Re**(-1/(0.83))

S = np.linspace(50,1e5,300)

plt.plot(S,fun_plane(S,2e7))
plt.plot(Re_kpc*1e3,Ie/((1e3)**2),'bo')
plt.xlabel("$R_e$ (pc)",fontsize=16)
plt.ylabel("$I_{total}()$",fontsize=16)
plt.xscale( 'log' )
plt.yscale( 'log' )
#Mass_Re_plot(Re_kpc,Ie,name,'bo','Bulge Mass estimated by us using stellar M/L ratio from Into et al.2013',1.0)
plt.show()
#Q = create_table("Dtable_Into_all_compare.txt", "Ctable_Into_all.txt", 4.53,103)


#cal_Mass(m_gal,DD,ML_ratio,M_sun)

##############
#un_list("BD_equvi_list.txt","Gal_bundle_BD_equvi", True)

run_list("equvi_list_bin2.txt","Gal_bundle_equvi_bin2",True)
run_list("equvi_list_bin3.txt","Gal_bundle_equvi_bin3",True)
run_list("equvi_list_bin4.txt","Gal_bundle_equvi_bin4",True)

C2,C3,C4 = cpt_seperator_demo('Gal_bundle_equvi_bin2'),cpt_seperator_demo('Gal_bundle_equvi_bin3'),cpt_seperator_demo('Gal_bundle_equvi_bin4')
#cpt_classifier_demo('Gal_bundle_equvi',C,'Gal_bundle_equvi_cpt',override_list_equ)

cpt_classifier_demo('Gal_bundle_equvi_bin2',C2,'Gal_bundle_equvi_bin2_cpt',override_list_equ)
cpt_classifier_demo('Gal_bundle_equvi_bin3',C3,'Gal_bundle_equvi_bin3_cpt',override_list_equ)
cpt_classifier_demo('Gal_bundle_equvi_bin4',C4,'Gal_bundle_equvi_bin4_cpt',override_list_equ)



name2 = grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin2.txt")["Gal_name"]
DD2 = grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin2.txt")["Dist"]
scale2 = grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin2.txt")["Scale"]

R_e_2 = grab_parameter("Gal_bundle_equvi_bin2_cpt", ["Bulge","CoreBulge"], 1) #get Re
mag2 = grab_mag("Gal_bundle_equvi_bin2_cpt", ["Bulge","CoreBulge"])

name3 = grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin3.txt")["Gal_name"]
DD3 = grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin3.txt")["Dist"]
scale3 = grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin3.txt")["Scale"]

R_e_3 = grab_parameter("Gal_bundle_equvi_bin3_cpt", ["Bulge","CoreBulge"], 1) #get Re
mag3 = grab_mag("Gal_bundle_equvi_bin3_cpt", ["Bulge","CoreBulge"])

name4 = grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin4.txt")["Gal_name"]
DD4 = grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin4.txt")["Dist"]
scale4 = grab_dist("./distance/dir_dist_list.txt","./distance/dist_list_bin4.txt")["Scale"]

R_e_4 = grab_parameter("Gal_bundle_equvi_bin4_cpt", ["Bulge","CoreBulge"], 1) #get Re
mag4 = grab_mag("Gal_bundle_equvi_bin4_cpt", ["Bulge","CoreBulge"])

CCC2 = np.genfromtxt("Ctable_Into_bin2_final.txt" , dtype='float')
CCC3 = np.genfromtxt("Ctable_Into_bin3_final.txt" , dtype='float')
CCC4 = np.genfromtxt("Ctable_Into_bin4_final.txt" , dtype='float')

AbsMag_g_bin2 = CCC2[:,6]
AbsMag_i_bin2 = CCC2[:,8]

AbsMag_g_bin3 = CCC3[:,6]
AbsMag_i_bin3 = CCC3[:,8]

AbsMag_g_bin4 = CCC4[:,6]
AbsMag_i_bin4 = CCC4[:,8]


##Taylor11_MassRatio

Mass_bulge_Into2 = cal_Mass(mag2,DD2,ML_relation_Iband(AbsMag_g_bin2,AbsMag_i_bin2).Into13_MassRatio,M_sun)
Re_kpc2 = R_e_2 * scale2

Mag2 = mag2 - 25 -5*np.log10(DD2)
Lum2 = 10**((Mag2-(4.53))/(-2.5))

Ie2 = Lum2/(np.pi*Re_kpc2**2)

Mass_bulge_Into3 = cal_Mass(mag3,DD3,ML_relation_Iband(AbsMag_g_bin3,AbsMag_i_bin3).Into13_MassRatio,M_sun)
Re_kpc3 = R_e_3 * scale3

Mag3 = mag3 - 25 -5*np.log10(DD3)
Lum3 = 10**((Mag3-(4.53))/(-2.5))

Ie3 = Lum3/(np.pi*Re_kpc3**2)

Mass_bulge_Into4 = cal_Mass(mag4,DD4,ML_relation_Iband(AbsMag_g_bin4,AbsMag_i_bin4).Into13_MassRatio,M_sun)
Re_kpc4 = R_e_4 * scale4

Mag4 = mag4 - 25 -5*np.log10(DD4)
Lum4 = 10**((Mag4-(4.53))/(-2.5))

Ie4 = Lum4/(np.pi*Re_kpc4**2)

SS=np.linspace(2e8,1e13,2000)


from plot_execute import Barro13_cut
from plot_execute import vDokkum15_cut
from plot_execute import vdWel14_cut
x_edge,y_edge= [0,0,2,2], [7e10,90e11,90e11,7e10]


plt.plot(Mass_bulge_Into2,Re_kpc2,'ro')
plt.plot(Mass_bulge_Into3,Re_kpc3,'bo')
plt.plot(Mass_bulge_Into4,Re_kpc4,'ko')

for i in range(np.size(name2)):
    plt.text(Mass_bulge_Into2[i],Re_kpc2[i], name2[i],fontsize=6)
for i in range(np.size(name3)):
    plt.text(Mass_bulge_Into3[i],Re_kpc3[i], name3[i],fontsize=6)
for i in range(np.size(name4)):
    plt.text(Mass_bulge_Into4[i],Re_kpc4[i], name4[i],fontsize=6)



plt.plot(SS,Barro13_cut(SS),"g--" , linewidth=3,label="Barro et al.2013" )
plt.vlines(1e10, 0, 10**((np.log10(1e10)-10.3)/1.5), linestyle="dashed", linewidth=3, color='g' )

plt.plot(SS,vDokkum15_cut(SS),"y--" , linewidth=3, label="van Dokkum et al.2015" )
plt.vlines(10**10.6, 0, 10**(np.log10(10**10.6)-10.7), linestyle="dashed", linewidth=3, color='y' )

plt.plot(SS,vdWel14_cut(SS),"b--" , linewidth=3, label="van der Wel et al.2014" )
plt.vlines(10**10.7, 0, 2.5*(((10**10.7)/1e11)**0.75), linestyle="dashed",linewidth=3, color='b' )

plt.vlines(7e10,0,2,color='k', linestyle="dashed", linewidth=3)
plt.hlines(2,7e10,90e11, color='k', linestyle="dashed", linewidth=3)
plt.fill(y_edge,x_edge, alpha=0.3, color='green')

plt.text(2e11, 30, "Bin2",fontsize=6)
plt.text(1e11,30, "Bin3",fontsize=6)
plt.text(5e10,30, "Bin4",fontsize=6)



plt.vlines(9e11,0,200,color='k', linestyle="dashed", linewidth=1)
plt.vlines(4e11,0,200,color='k', linestyle="dashed", linewidth=1)
plt.vlines(2e11,0,200,color='k', linestyle="dashed", linewidth=1)
plt.vlines(1e11,0,200,color='k', linestyle="dashed", linewidth=1)
plt.vlines(5e10,0,200,color='k', linestyle="dashed", linewidth=1)


plt.xscale( 'log' )
plt.yscale( 'log' )
plt.show()
    #Mass_Re_plot(Mass_bulge_Into2,Re_kpc2,name,'ro','Bulge Mass estimated by us using stellar M/L ratio from Into et al.2013',1.0)
#Mass_Re_plot(Mass_bulge_Into3,Re_kpc3,name,'bo','Bulge Mass estimated by us using stellar M/L ratio from Into et al.2013',1.0)
#Mass_Re_plot(Mass_bulge_Into4,Re_kpc4,name,'ko','Bulge Mass estimated by us using stellar M/L ratio from Into et al.2013',1.0)



##############
#run_list("BD_major_list.txt","Gal_bundle_BD_major", False)
#run_list("BD_equvi_list.txt","Gal_bundle_BD_equvi", True)

#run_list("major_list_2.txt","Gal_bundle_major",False)
#run_list("equvi_list_2.txt","Gal_bundle_equvi",True)
#
#A1, A2 = cpt_seperator_demo('Gal_bundle_BD_major'), cpt_seperator_demo('Gal_bundle_BD_equvi')
#B, C = cpt_seperator_demo('Gal_bundle_major'), cpt_seperator_demo('Gal_bundle_equvi')
#
#cpt_classifier_demo('Gal_bundle_BD_major',A1,'Gal_bundle_BD_major_cpt',[])
#cpt_classifier_demo('Gal_bundle_BD_equvi',A2,'Gal_bundle_BD_equvi_cpt',[])
#
#cpt_classifier_demo('Gal_bundle_major',B,'Gal_bundle_major_cpt',override_list_maj)
#cpt_classifier_demo('Gal_bundle_equvi',C,'Gal_bundle_equvi_cpt',override_list_equ)
#
lookup_bundle("Gal_bundle_equvi","UGC8736")
lookup_bundle("Gal_bundle_equvi_cpt","UGC8736")

lookup_bundle("Gal_bundle_equvi","NGC4826")
lookup_bundle("Gal_bundle_equvi_cpt","NGC4826")

#
#############

#print('-----')

#print("master",A["master_Disk_can1"], A["master_Disk_can1_row"], A["master_Disk_can1_index"])
#print('A',A["master_Disk_can1"][0][0][1])
#print('-----')

#cpt_classifier_demo('Gal_bundle_BD_major',A,'Gal_bundle_BD_major_cpt',[])


#print('-----')
#print('cpt',read_list('Gal_bundle_major_cpt'))
#print('-----')

#print('cpt',read_list('Gal_bundle_equvi_cpt'))

###############
plot_compare_rms('Gal_bundle_major_cpt','Gal_bundle_BD_major_cpt')
plt.show()

#compare_para3('Gal_bundle_major_cpt','Gal_bundle_equvi_cpt', "Bulge", 2, "n_{major}","n_{equvi}")
#plt.show()
#
#compare_para3('Gal_bundle_major_cpt','Gal_bundle_equvi_cpt', "PrimBar", 1, "r_{out,major}","r_{out,equvi}")
#plt.show()

#print(read_list("Gal_bundle_cpt"))

#print(read_list("Gal_bundle_equvi_cpt"))

plot_hist_percentage("Gal_bundle_equvi_cpt",DD)
plt.show()


plot_hist_percentage("Gal_bundle_equvi_bin2_cpt",DD2)
plt.show()

plot_hist_percentage("Gal_bundle_equvi_bin3_cpt",DD3)
plt.show()
plot_hist_percentage("Gal_bundle_equvi_bin4_cpt",DD4)
plt.show()
##############
#fig, ax = plt.subplots()
#compare_para('Gal_bundle_BD_major','Gal_bundle_BD_equvi', 6, 0, '$\mu_{0(Exp)}$')
#plt.show()
#
#fig, ax = plt.subplots()
#compare_para('Gal_bundle_BD_major','Gal_bundle_BD_equvi', 3, 0, '$\mu_{0(Sersic)}$')
#plt.show()
#
#
#fig, ax = plt.subplots()
#compare_para('Gal_bundle_BD_major','Gal_bundle_BD_equvi', 3, 1, '$R_{e(Sersic)}$')
#plt.show()
#
#fig, ax = plt.subplots()
#compare_mag('Gal_bundle_BD_major','Gal_bundle_BD_equvi', -1, '$Mag_{total}$')
#plt.show()
###############
