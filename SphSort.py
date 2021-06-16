#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:46:49 2020

@author: dexter
"""

import numpy as np
import pickle

import matplotlib.pyplot as plt
from numpy.linalg import inv
import os


from SphRead import *


__all__ = ["remove_value","trim_value","replace_value",
           "str_zipping","remove_item", "add_item", "outliers",
           "selection_generic","cherry_pick" ,"cpt_seperator_demo",
           "cpt_classifier_demo", "plus_minus_seperator", "vdis_match", 
           "LTG_ETG_seperator", "prop_seperation"]

__author__="Dexter S.H. Hon"

#%% tested
def remove_value(input_array,value_list):
    """
    Remove the sepcfied value in an numpy array.

    Parameters
    ----------
    input_array : 1D numpy array
        The input.
    value_list : list
        The indices of value to be removed .

    Returns
    -------
    An resized numpy array without the value inside.


    """
    value = []
    input_array = list(input_array)
    
    for j in range(len(value_list)):
        value.append(input_array[value_list[j]])
        
    for i in range(len(value)):  
        input_array.remove(value[i])

    temp=np.array(input_array)        
    
    return temp

#%% tested
def trim_value(input_array,value):
    """
    Trim down an numpy array by removing a specific value.

    Parameters
    ----------
    input_array : 1D numpy array
        The input.
    value : str, float
        The value to be removed.

    Returns
    -------
    An resized numpy array without the value inside.

    """
    temp = []
    
    for i in range(np.size(input_array)):
        
        if input_array[i] == value:
            pass
        else:
            temp.append(input_array[i])
    temp=np.array(temp)        
    
    return temp

#%% tested
def replace_value(input_array,value1,value2):
    """
    Replace the sepcfied value in an numpy array.

    Parameters
    ----------
    input_array : 1D numpy array
        The input.
    value1 : str, float
        The value to be replaced.
    value2 : str, float
        The value to substitue.

    Returns
    -------
    An new numpy array with the replaced value inside.

    """
    temp = []
    
    for i in range(np.size(input_array)):
        
        if input_array[i] == value1:
            temp.append(value2)
        else:
            temp.append(input_array[i])
    temp=np.array(temp)        
    
    return temp

#%% tested
def str_zipping(list1,list2,zip_symbol=""):
    """
    Convert the content of two lists into str and link them together.
    e.g.: str_zipping([1,2,4],[4,5,3],zip_symbol = '+')
        >>> ['1+4','2+5','4+3']

    Parameters
    ----------
    list1, list2 : list
        The lists to be conjoin. They must have the same dimensions

    zip_symbol : str
        The symbol in which you want to add in between the elements.
        The default is "".
        
    Returns
    -------
    A list of conjoined elements.

    """
    
    list1, list2 = list(list1), list(list2)
    
    if len(list1) == len(list2): # matching list dimension
        pass
    else:
        ValueError
    
    list_product = []
    
    for i in range(len(list1)):
        list_product.append(str(list1[i])+zip_symbol+str(list2[i]))
        
    return list_product

#%%tested
def str_zipping_generic(*args):
    """
    Convert the content of multiple lists into str and link them together.
    e.g.: str_zipping('(',[1,2,4], '+', [4,5,3],')')
        >>> ['(1+4)','(2+5)','(4+3)']

    Parameters
    ----------
    *args: list or str or int
        The objects to be conjoin. The lists must have the same dimensions. 
        Otherwise, any str or int should be one per entry.


    Returns
    -------
    A list that looks like:
        ['(1+4)','(2+5)','(4+3)']
    """
    len_doc, master_list = [], []
    list_product = []

    for i in args:
        if type(i) == list:
            len_doc.append(len(i))
    #check input length consistency           
    for i in len_doc:
        if len_doc == i:
            pass
        else:
            ValueError
   #multiply singular element into list
    for i in args:
        if type(i) == list:
            master_list.append(i)
        elif type(i) == str: 
            A = []
            for j in range(len_doc[0]):
                A.append(i)
            master_list.append(A)
    #zipping    
    for i in range(len_doc[0]):
        
        string = ''
        for j in range(len(master_list)):
            string = string + str(master_list[j][i])            
            
        list_product.append(string)
        
    return list_product


#%% tested
def remove_item(input_list_name,keyword_list):
    """
    A method to remove lists base on specified keywords, 
    from the galaxy bundle.
    ----------                        
    input_list_name: str
        The galaxy bundle
        
    keyword_list: list
        A list containing the names of the galaxy
        
    Return
    ------
    Galaxy Bundle
    """    
    sample = read_list(input_list_name)
    row_list=[]
    row_example=[]

    for item in keyword_list:
        keyword = item

        for row in range(len(sample)):

            if sample[row][0] == keyword:
                print("Bingo",row,keyword)
                row_list.append(row)
                row_example.append(sample[row])
                pass
            else:
                pass
            
    sample2 = sample  
    
    for i in range(len(row_example)):
        sample2.remove(row_example[i])
    
    #print(row_list)
    #for i in row_list:
    ##    print(sample[i])
    #    sample2.remove(sample[i])
    ##    #(sample[i])
    #    print(i)
    #    print(len(sample2),len(sample))

            
    with open(input_list_name, 'wb') as f:
        pickle.dump(sample2, f)
        
    return sample2

#%% tested
def add_item(parent_list_name, target_list_name, keyword_list):
    """
    A method to add lists into a bundle based on specified keywords.
    ----------                        
    parent_list_name: str
        The parent galaxy bundle.
        
    target_list_name: str
        The target bindle.
    
    keyword_list: list
        A list containing the names of the galaxy.
        
    Return
    ------
    New target list
    """    
    parent_sample = read_list(parent_list_name)
    target_sample = read_list(target_list_name)
    
    for item in keyword_list:
       keyword = item
       for row in range(len(parent_sample)):
            if parent_sample[row][0] == keyword:
                target_sample.append(parent_sample[row])
            else:
                pass
            
    with open(target_list_name, 'wb') as f:
        pickle.dump(target_sample, f)
    return target_sample

#%%
def outliers(name,input_array,limit, direction="large"):
    """
    A generic function to select outliers larger than the limit
    ----------                        
    name : str
        The name of the outlier
            
    input_array : 1D numpy array         
        The array in question
            
    limit : float
        The outlier cut. 
        
    direction: str
        Define which direction of outlier to select.
        The options are "large" or "small".
        The default is "large".

    Return
    ------
    Dictionary containing the selected sample name and values
        
    """
    
    #the outlier name and attribute
    outliers_name, outliers =[], []
    if direction == "large":
        for row in range(len(input_array)):
            if (input_array[row] > limit):
                outliers_name.append(name[row])
                outliers.append(input_array[row])
            else: 
                pass
    elif direction == "small":
        for row in range(len(input_array)):
            if (input_array[row] < limit):
                outliers_name.append(name[row])
                outliers.append(input_array[row])
            else: 
                pass    
            
    return {"name":outliers_name, "Dist":outliers}

#%% tested
def selection_generic(input_list_x, input_list_y, func, 
                      direction="low", axis="y"):
    """
    A generic method (2D) to select a set of sample base on a cut.
        
    Parameters
    ----------
    input_list: 1D numpy array
        The quantity to compare the cut with.
        
    func: 1D numpy array
        An function with the same dimension as the input list.
    
    direction: str
        The direction of the selection
        The options are either "low" or "high".
        The default is "low".
        
    axis: str
        The axis of interest.
        The options are either "x" or "y"
        The default is "y"
    Returns
    -------
    Bag: Dict
    {index:[1,4,5],
             bag:[[],[]...]}

    """
    index_list,bag_list_y, bag_list_x=[],[],[]
    
    match_list_dim(input_list_x, func)
    match_list_dim(input_list_x, input_list_y)
    
    if axis == "y":
        if direction == "high":
            for i in range(len(input_list_y)):
                if input_list_y[i] > func[i]:
                
                    index_list.append(i)
                    bag_list_x.append(input_list_x[i])
                    bag_list_y.append(input_list_y[i])
                    
        elif direction == "low":
            for i in range(len(input_list_y)):
                if input_list_y[i] < func[i]:
                
                    index_list.append(i)
                    bag_list_x.append(input_list_x[i])
                    bag_list_y.append(input_list_y[i])  
                    
    elif axis == "x":
        if direction == "high":
            for i in range(len(input_list_x)):
                if input_list_x[i] > func[i]:
                
                    index_list.append(i)
                    bag_list_x.append(input_list_x[i])
                    bag_list_y.append(input_list_y[i])
                    
        elif direction == "down":
            for i in range(len(input_list_x)):
                if input_list_x[i] < func[i]:
                
                    index_list.append(i)
                    bag_list_x.append(input_list_x[i])
                    bag_list_y.append(input_list_y[i]) 
        
    Bag = {'index':index_list,
           'bag_x':bag_list_x,
           'bag_y':bag_list_y
               }
    return Bag 

#%%
# need to make zone_in_generic
# that does not loop every element
def zone_in_2D(input_matrix,low_bound=None,up_bound=None):
    """
    A function to select the sample within 
    low_bound and up_bound in 2D.
    
    This function loop through every elements.
    
    Parameters
    ----------
    input_matrix : list 
        The 2D input list.
    low_bound : list
        The list of lower bound.
        This need to be the same length as the 
        input_matrix.
        The default is [-1.0*np.inf,...]
    up_bound : TYPE
        The list of lower bound.
        This need to be the same length as the 
        input_matrix.
        The default is [-1.0*np.inf,...].

    Returns
    -------
    The dictionary that contain the indices, and the data matrix.

    """
    index_list, bag_list, bag_list_y, bag_list_x=[],[],[], []
    bag = {}
    
    dimension = len(input_matrix)
    
    # produce array of inf
    if low_bound == None:
        low_bound = np.repeat(-1.0*np.inf,dimension)

    if up_bound == None:
        up_bound = np.repeat(np.inf,dimension)

    # check if the dimension matches
    match_list_dim(low_bound, input_matrix)
    match_list_dim(up_bound, input_matrix)  
    
    for j in range(len(input_matrix[0])):
        if input_matrix[0][j] > low_bound[0] and \
            input_matrix[1][j] > low_bound[1] and \
            input_matrix[0][j] < up_bound[0] and \
                input_matrix[1][j] < up_bound[1]:
            
            index_list.append(j)
            bag_list_x.append(input_matrix[0][j])
            bag_list_y.append(input_matrix[1][j])
            
    bag_list.append(bag_list_x)
    bag_list.append(bag_list_y)
    
    bag["index"] = index_list
    bag["data"] = bag_list
    return bag

#%% tested
def seperator_label_generic(bundle, 
                            compartment_label,label_entry=None):
    """
    A generic method to compartmentalise a bundle.
    This function convert a parent list bundle into a list that contain 
    a number of dictionaries, based on labels.
    
    Example:
        [[name1,label1,value1,label2,value2,...],[...],...]
        -->input compartment_label = ["A","B","C",..]
        -->[{"A_index": [...],
             "A_sample":[[...],...]},
            {"B_index":[...].
             "B_sample":[[...],...] }...]

    Parameters
    ----------
    bundle : list
        A list bundle.
        
    compartment_label : list
        A list of label for matches.
        
    label_entry : float or int or list
        The location of the label entries.
        The default is None. We loop through the whole sample entry to FIND 
        the desired label.
        
        If the input is float, we assume all label_entries are in the 
        same position. e.g. Second column, input:1
        If the input is a list, we loop through the list and extract 
        information based on the location provided by a list. 
        e.g. [1,3,6], first loop: 2 nd column, second loop: 4th column, and 
        third loop, 6th column.
        

    Returns
    -------
    new_bundle.
        A new bundle that contain the compartmentalised data.
    """

    
    # check the nature of the label_entry
    
    new_bundle = []
    
    # loop through the list of labels
    for i in range(len(compartment_label)):
        Dict ={}
        
        index, data = [], []
        
        # Loop through everything if label_entry is None
        if label_entry == None:
            for j in range(len(bundle)):
                for k in range(len(bundle[j])):
                    if bundle[j][k] == compartment_label[i]:
                        index.append(j)
                        data.append(bundle[j])
                        pass
                    
        # Look up a fixed column if label_entry is int or float
        elif type(label_entry) == float or type(label_entry) == int:
            for j in range(len(bundle)):
                if bundle[j][int(label_entry)] == compartment_label[i]:
                    index.append(j)
                    data.append(bundle[j])
                    pass      
                
        # Look up specific column according to label_entry list
        elif type(label_entry) == list:
            #check if the dimension of label_entry and compartment_label matches
            match_list_dim(compartment_label, label_entry)
            
            for j in range(len(bundle)):
                if bundle[j][i] == compartment_label[i]:
                    index.append(j)
                    data.append(bundle[j])
                    pass       
                
        Dict[compartment_label[i]+"_index"] = index
        Dict[compartment_label[i]] = data
        
        new_bundle.append(Dict)
        
    return new_bundle

#%% tested
def cherry_pick(index_list,parent_list):
    """
    A generic method to cherry pick a set of subsample, given a list of indices 
    to select from.

    Parameters
    ----------
    index_list : 1D list
        The indices of the subsample.
    parent_list : 1D list
        The material of which one select from.

    Returns
    -------
    subsample : 1D list
        The result subsample.

    """
    subsample = []
    i=0
    for i in range(len(index_list)):
        subsample.append(parent_list[index_list[i]])
    return subsample

#%% tested
def morph_str_selection(index_list,morph_list):
    """
    Seperate an index list by the morphology type E, S0, or S string indicator.

    Parameters
    ----------
    index_list : 1d list
        The target index list.
    morph_list : 1D list
        The corresponding morphology list to the index list.
        
    Returns
    -------
    morph_dict: dict
        example
        {"E": [index_list[1],index_list[6],index_list[9],...],
         "S0": [index_list[2],index_list[4],index_list[8],...],
         "S": [index_list[3],index_list[6],index_list[10],...] }

    """
    morph_dict = {}
    E_list,S0_list,S_list = [], [], []
    for i in range(len(index_list)):
        if "E" in morph_list[i]:
            E_list.append(index_list[i])
        elif "0" in morph_list[i]:
            S0_list.append(index_list[i])
        else:
            S_list.append(index_list[i])

    morph_dict["E"] = E_list
    morph_dict["S0"] = S0_list
    morph_dict["S"] = S_list

    return morph_dict

#%% tested
def cpt_seperator_demo(input_list_name): 
    """
    A demo function for seperating galaxy components to different bins.
    It create the candidate list for component assignment.
    ...

    Parameters
    ----------
    input_list_name : str
        The file name of the galaxy bundle.

    Return
    -------

    
    """
    sample = read_list(input_list_name) #sample package
    
    master_Bulge_can, master_Bulge_can_index = [], []
    master_CoreBulge_can, master_CoreBulge_can_index = [], []

    master_Disk_can1, master_Disk_can2, master_Disk_can3 = [],[],[]
    master_Disk_can1_index, master_Disk_can2_index, master_Disk_can3_index = [],[],[]
    master_Bar_can, master_Bar_can_index = [], []
    master_Ring_can, master_Ring_can_index = [],[]
    master_Total_mag_can, master_Total_mag_can_index = [], []
    
    master_Bulge_can_row = []
    master_CoreBulge_can_row = []
    master_Disk_can1_row = []
    master_Disk_can2_row = []
    master_Disk_can3_row = []
    master_Bar_can_row = []
    master_Ring_can_row = []
    master_Total_mag_can_row = []
    
    for number in range(len(sample)):
        
        #print(sample[number])
        sample_indi = sample[number]    #indivdual sample
        
        Bulge_can,Bulge_can_index = [], []
        CoreBulge_can,CoreBulge_can_index = [],[]
        Disk_can1, Disk_can2, Disk_can3, Disk_can1_index, Disk_can2_index, Disk_can3_index = [],[],[],[],[],[]
        Bar_can, Bar_can_index = [], []
        Ring_can, Ring_can_index = [],[]
        Total_mag_can, Total_mag_can_index = [],[]
                  
                
        Bulge_can_row = []
        CoreBulge_can_row = []

        Disk_can1_row = []
        Disk_can2_row = []
        Disk_can3_row = []
        Bar_can_row = []
        Ring_can_row = []
        Total_mag_can_row = []
        
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
            # Anse, Rings
            elif sample_indi[index] == "Total_mag":
                Total_mag_can.append(sample_indi[index+1])
                Total_mag_can_index.append(index)  
                Total_mag_can_row.append(number)
                #print(sample_indi[0],"Gauss")
            ################################
        master_Bulge_can.append(Bulge_can), 
        master_Bulge_can_index.append(Bulge_can_index) 
        master_CoreBulge_can.append(CoreBulge_can), 
        master_CoreBulge_can_index.append(CoreBulge_can_index) 
        master_Disk_can1.append(Disk_can1), 
        master_Disk_can1_index.append(Disk_can1_index)
        master_Disk_can2.append(Disk_can2), 
        master_Disk_can2_index.append(Disk_can2_index)
        master_Disk_can3.append(Disk_can3), 
        master_Disk_can3_index.append(Disk_can3_index)
        master_Bar_can.append(Bar_can), 
        master_Bar_can_index.append(Bar_can_index)
        master_Ring_can.append(Ring_can), 
        master_Ring_can_index.append(Ring_can_index)
        
        master_Total_mag_can.append(Total_mag_can), 
        master_Total_mag_can_index.append(Total_mag_can_index)
        
        master_Bulge_can_row.append(Bulge_can_row)
        master_CoreBulge_can_row.append(CoreBulge_can_row)
        master_Disk_can1_row.append(Disk_can1_row)
        master_Disk_can2_row.append(Disk_can2_row)
        master_Disk_can3_row.append(Disk_can3_row)
        master_Bar_can_row.append(Bar_can_row)
        master_Ring_can_row.append(Ring_can_row)
        master_Total_mag_can_row.append(Total_mag_can_row)

    return{"master_Bulge_can":master_Bulge_can, 
           "master_Bulge_can_index":master_Bulge_can_index, 
           "master_Bulge_can_row":master_Bulge_can_row,
           "master_CoreBulge_can":master_CoreBulge_can, 
           "master_CoreBulge_can_index":master_CoreBulge_can_index, 
           "master_CoreBulge_can_row":master_CoreBulge_can_row,
           "master_Disk_can1":master_Disk_can1, 
           "master_Disk_can1_index":master_Disk_can1_index, 
           "master_Disk_can1_row":master_Disk_can1_row,
           "master_Disk_can2":master_Disk_can2, 
           "master_Disk_can2_index":master_Disk_can2_index,
           "master_Disk_can2_row":master_Disk_can2_row,
           "master_Disk_can3":master_Disk_can3, 
           "master_Disk_can3_index":master_Disk_can3_index,
           "master_Disk_can3_row":master_Disk_can3_row,
           "master_Bar_can":master_Bar_can, 
           "master_Bar_can_index":master_Bar_can_index,
           "master_Bar_can_row":master_Bar_can_row,
           "master_Ring_can":master_Ring_can, 
           "master_Ring_can_index":master_Ring_can_index, 
           "master_Ring_can_row":master_Ring_can_row,
           "master_Total_mag_can":master_Total_mag_can, 
           "master_Total_mag_can_index":master_Total_mag_can_index, 
           "master_Total_mag_can_row":master_Total_mag_can_row}            
                    
#%% tested
def cpt_classifier_demo(input_list_name, input_sep_dict ,output_list_name, 
                        override_instruction=[]):
    """
    A demo function for interpret the nature of the analytic function 
    describing a galaxy.
    
    ...
    The scheme of interpretation is as following:
        
        "Bulge": 
            The galactic spheroid.
            Bulges are described by Sersic function.
            In the presecne of more than one Sersic function in a galaxy,
            The one with the largest sersic index is considered the bulge.
            
        "Lens":
            Lens are described by Sersic function (usually n<1).
            In the presecne of more than one Sersic function in a galaxy,
            The one with the smallest sersic index is considered the lens.
            
        "CoreBulge":
            The core-depleted galactic spheroid.
            CoreBulge are desribed by core-Sersic function
            
        "Disk": 
            The extended Disk.
            The disk is described by either an expenetial function, a broken 
            exponential function, or a inclined disk model.
            
        
        "nucDisk": 
            The nuclear Disk.
            The nuclear disk is desribed by an exponential function with very
            steep slope.
        
        "PrimBar":
            The primary bar
        
        "SecBar":
            The secondary bar
        
        "Ring":
            The ring, anse and spiral arms
        
        
    
    
    This function replace the analytic functions label with the component name.
    
    The interpretation can be overrided in special case. 
    
    ...

    Parameters
    ----------
    input_list_name : str
        The file name of the galaxy bundle, from cpt_seperator_demo
        
    input_sep_dict: python dict
        A dictionary seperation 
    
    output_list_name: str
        The name of the output file. 
    
    Optional
    --------
    override_instruction: list
        A list for user 
    
    Return
    -------

    """
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
        Disk_can2_index = sep_dict["master_Disk_can2_index"][number]
        Disk_can2_row = sep_dict["master_Disk_can2_row"][number] 
                
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
        Disk_can3_index = sep_dict["master_Disk_can3_index"][number]
        Disk_can3_row = sep_dict["master_Disk_can3_row"][number]
        
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
        Bar_can_index = sep_dict["master_Bar_can_index"][number]
        Bar_can_row = sep_dict["master_Bar_can_row"][number]
        
        #print("Gal",sample[number][0], number, "Ferrer")
        if Bar_can == []:
            pass
        else:
            if len(Bar_can) > 1:
                r_out_list = [Bar_can[i][2] for i in range(len(Bar_can))]
                
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
        Ring_can_index = sep_dict["master_Ring_can_index"][number]
        Ring_can_row = sep_dict["master_Ring_can_row"][number]
        
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
                    print("replaced!", sample[i][cpt_index], 
                          override_instruction[number+2])
                else:
                    pass
    #############
    with open(output_list_name, 'wb') as f:
        pickle.dump(sample, f)
        
    return sample

#%%
def bundle_creation(obj_list,equvi=False,override_list=[]):
    """
    A method to create galaxy bundle, in an single execution.
    If you want to save the in
    ----------                        
    obj_list : str
        The name of ASCII file containing the file names of the profiler logs.
        
        e.g. 
        ./my_directory/galaxy1/galaxy1_major.txt
        ./my_directory/galaxy2/galaxy2_major.txt
        ./my_directory/galaxy3/galaxy3_major.txt

    equvi : bool         
        The value indicating if the logs are in equvilant axis
            
    override_list : list
        A list of overriding instruction.

        e.g.
        [name, index of element, replacement item,...]
        ["NGC2862",2,"Disk","NGC2872",5,"Point Source" ]
        
    Return
    ------
    bundle: list
        A galaxy bundle with the information of each structural components.       
    """
    run_list(obj_list,"gal_temp",equvi)

    C= cpt_seperator_demo('gal_temp')
    bundle = cpt_classifier_demo('Gal_bundle_equvi_bin2',C,
                        'Gal_bundle_equvi_bin2_cpt',override_list)
    os.remove("gal_temp")
    return bundle    

#%% tested
def plus_minus_seperator(name,input_array,limit):
    """
    A method to select seprate positive and negative number
    ----------                        
    name : str
        The name of the galaxy.
            
    input_array : 1D numpy array         
        The array in question
            
    limit : float
        The outlier cut. 
        
    Return
    ------
    A dictionary containing the selected sample name and values
        
    """
    
    outliers_name_positive, outliers_name_negative =[], [] 
    outliers_positive, outliers_negative = [], []
    
    for row in range(len(input_array)):
        if (input_array[row] < 0):
            outliers_name_negative.append(name[row])
            outliers_negative.append(input_array[row])
        elif (input_array[row] > 0): 
            outliers_name_positive.append(name[row])
            outliers_positive.append(input_array[row])
            
    return {"name_positive":outliers_name_positive, 
            "Dist_positive":outliers_positive,
            "name_negative":outliers_name_negative, 
            "Dist_negative":outliers_negative}

#%% tested
def LTG_ETG_seperator(input_list_name,output_list_name_LTG,
                      output_list_name_ETG):    
    ## correction, LTG and ETG are selected by if there's a spiral arm
    """
    A method to seprate Early-Type and Late-Type Galaxy.
    The seperation condition is whether the galaxy contain a extended disk.
    ----------                        
    input_lis_name : str
        The galaxy component bundle 
            
    output_list_name1,2 : str      
        The name for the output dictionary

        
    Return
    ------
    Dictionary containing two list: ETG and LTG.
    """
    input_list = read_list(input_list_name)
    ETG, LTG = [],[]
    
    for row in range(len(input_list)):
        #for index in range(len(input_list[row])):
        if "Ring" in input_list[row]:
            #print("LTG",input_list[row][0])
            LTG.append(input_list[row])
        elif "Ring" not in input_list[row]:
            #print("ETG",input_list[row][0])
            ETG.append(input_list[row])
            
    output = {"ETG": ETG, "LTG": LTG}
    
    with open(output_list_name_LTG, 'wb') as f:
        pickle.dump(output["LTG"], f)
    with open(output_list_name_ETG, 'wb') as f:
        pickle.dump(output["ETG"], f)
    return output    
#%% tested
def vdis_match(input_list_name, vdis_list, Dist_list, output_list_name):
    """
    A function to match the existing velocity dispersion with galaxy property 
    ----------                        
    input_list_name:
        Galaxy bundle with all the components
    vdis_list:
        
    Dist: py list, 2 name
        e.g. ["./distance/dir_dist_list.txt","./distance/dist_list.txt"]
    output_list_name:
        
    Return
    ------
    A dictionary with galaxy essential inforamtion 
    """
    input_list = read_list(input_list_name)
    
    sph_mag_p = grab_mag(input_list_name, ["Bulge","CoreBulge"])
    sph_Re_p = grab_parameter(input_list_name, ["Bulge","CoreBulge"], 1)
    profiler_total_mag_p = grab_total_mag(input_list_name)
    
    V1= np.genfromtxt(vdis_list, dtype='str') 
    V2= np.genfromtxt(vdis_list, dtype='float') 
    
    vdis_name_p, vdis_p, vdis_err_p = V1[:,0], V2[:,3], V2[:,4] #reading vdis file, specfic col numbers
    mag_g_p, mag_i_p = V2[:,10], V2[:,9]
    
    Dist_p = grab_dist(Dist_list[0],Dist_list[1])["Dist"]
    Dist_name_p = grab_dist(Dist_list[0],Dist_list[1])["Gal_name"]

    #Gal_bundle_BD_equvi

    sph_mag, sph_Re = [],[]
    name, vdis, vdis_err, Dist = [], [] ,[], []
    mag_g, mag_i = [], []
    profiler_total_mag =[]

    for row in range(len(input_list)):
        name_0 = input_list[row][0]
        
        for i in range(len(vdis_name_p)):
            if vdis_name_p[i] == name_0:
                name.append(name_0)
                vdis.append(vdis_p[i])
                vdis_err.append(vdis_err_p[i])
                
                mag_g.append(mag_g_p[i])
                mag_i.append(mag_i_p[i])
                
                profiler_total_mag.append(profiler_total_mag_p[row])
                
                sph_mag.append(sph_mag_p[row])
                sph_Re.append(sph_Re_p[row])                
            else:
                pass
            
        for j in range(len(Dist_name_p)):
            if  Dist_name_p[j] == name_0:
                Dist.append(Dist_p[j])
            else:
                pass
            
    vdis, vdis_err, Dist = np.array(vdis), np.array(vdis_err), np.array(Dist)
    mag_g, mag_i = np.array(mag_g), np.array(mag_i)
    profiler_total_mag = np.array(profiler_total_mag)
    
    sph_mag = np.array(sph_mag)
    sph_Re = np.array(sph_Re)
    
    output = {"name": name, 
              "vdis": vdis,
              "vdis_err": vdis_err,
              "Dist": Dist,
              "mag_g": mag_g,
              "mag_i": mag_i,
              "total_mag": profiler_total_mag,
              "sph_mag": sph_mag,
              "sph_Re": sph_Re
              }
    
    with open(output_list_name, 'wb') as f:
        pickle.dump(output, f)
    return output
    
#%% Under Construction no one allow in
    
    #select base on mass bin 
    #specific feature bin
def prop_seperation():
    return None
# input morph info, 
# select feature, scan, create a dictionary of string, loop through it and create a set of lists correpsonding to the feature