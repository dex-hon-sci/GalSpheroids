#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:46:49 2020

@author: dexter
"""

import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from numpy.linalg import inv

from SphRead import read_list

#%% tested
def cpt_seperator_demo(input_list_name): 
    """
    A demo function for seperating galaxy components to 
    ...

    Parameters
    ----------
    input_list_name : str
        The file name of the galaxy bundle.

    Return
    -------

    
    """
    ####################################################################################################
    ## replacing analytic function label to component label
    ## 
    ####################################################################################################
    ## input_list_name:
    ## output_list_name:
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
            
                    
#%% tested
def cpt_classifier_demo(input_list_name, input_sep_dict ,output_list_name, override_instruction):
    
    sample = read_list(input_list_name) #sample package
    sep_dict = input_sep_dict #the result of sperator, a dictionary of each cpt, with row and index
    #A = cpt_seperator_demo("Gal_bundle")

    
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
                
            elif len(Bar_can) == 1:
                row = Bar_can_row[0]
                index = Bar_can_index[0]
            
                sample[row][index]= "PrimBar"            
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
    #select base on mass bin 
    #morph Bin
    #specific feature bin
def prop_seperation():
    return None
# input morph info, 
# select feature, scan, create a dictionary of string, loop through it and create a set of lists correpsonding to the feature