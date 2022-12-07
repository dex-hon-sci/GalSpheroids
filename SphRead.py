#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:42:10 2020

@author: dexter
"""

#import packages

import numpy as np
import pandas as pd
import pickle

import astropy as astropy 
#from astropy.table import Table
#from astropy.io import ascii

__all__ = ["read_list", "read_table", "match_list_dim", "convert_dict_ascii", 
           "convert_list_textable","lookup_bundle", "grab_name",
           "grab_parameter","grab_mag","grab_total_mag", "grab_dist", 
            "extract_match",
           "pd_read","run_list"]

__author__="Dexter S.H. Hon"

###########initial parameters###############
c, H0=    299792.458, 68.0        #speed of light in km/s, Hubble constant in(km/s)/Mpc
############################################
#%% tested
def read_list(input_file):
    """
    Just a function using pickle to load objects
    """
    with open(input_file, 'rb') as f:
        mylist = pickle.load(f)
    return mylist

#%% tested
def read_table(table,dtype='float'):
    """
    Just a function using genfromtext to load objects
    Assuming the data is float.
    """
    D=np.genfromtxt(table, dtype=dtype)
    return D
#%%
def pickle_save(input_list,output_name):
    """
    Simple save function via pickle.

    Parameters
    ----------
    input_list: list
        A python list
        
    output_name: str
        The filename for the dictionary to be saved via pickle
    
    Returns
    -------
    None.

    """
    with open(output_name, 'wb') as f: # saving the list 
        pickle.dump(input_list, f)
        
    return None
#%% tested # under construction
def match_list_dim(input1,input2):
    """
    A generic method matching the dimension of two input
    
    Parameters
    ----------
    input1, input2: list

    """
    if len(input1) == len(input2):
        pass
    else:
        raise Exception("Dimension not matching, \
                        size of %s =\= %s" %((input1),(input2))) 
#%%tested
def convert_list_ascii(input_name,output_name):
    """
    Convert python dictionary or list to a ASCII file.
    
    Parameters
    ----------
    Input_name: list
        The name of the list created by pickle.
        
    output_name: str
        The name of the output ASCII file.
    
    """    
    #with open(input_file, 'rb') as f:
    #    mylist = pickle.load(f)

    #print(value)
    #print(key)
    my_list = input_name
    with open(output_name, 'w') as f:
        for item in my_list:
            f.write("%s\n" % item)
    

#%%
def convert_dict_ascii(input_name,output_name):
    """
    Convert python dictionary or list to a ASCII file.
    
    Parameters
    ----------
    Input_name: str
        The name of the dict created by pickle.
        
    output_name: str
        The name of the output ASCII file.
    
    """
    input_file = read_list(input_name)
    
    #with open(input_file, 'rb') as f:
    #    mylist = pickle.load(f)

    value = list(input_file.values())
    key = list(input_file.keys())

    #print(value)
    #print(key)

    data = astropy.table.Table(value, names=key)
    astropy.io.ascii.write(data, output_name ,overwrite=True)
    
#%%WIP
def convert_list_textable(input_file, output_name):
    """
    Export tex format table from a python list to an ascii file

    Returns
    -------
    None.

    """
    
    with open(input_file, 'wb') as f:
        mylist = pickle.load(f)
        
    value = mylist.values()
    key =mylist.key()
    
    data = astropy.table.Table(value, names=key)

    astropy.io.ascii.write(data, output_name ,overwrite=True)

        
    return None

#%% tested
def lookup_bundle(gal_bundle,gal_name):
    
    """
    Lookup individual galaxy by name in a galaxy bundle.
    
    Parameters
    ----------
    gal_bundle: list
        The galaxy bundle.
        
    gal_name: str
        The name of the query 
        
    Return
    ----------
    Gal: list
    """
    GN= read_list(gal_bundle)
    for row in range(len(GN)):
        if GN[row][0]==gal_name:
            Gal = GN[row]
            print(GN[row])
        else:
            pass
    return Gal

#%% tested
def grab_name(filename):
    """
    A function for grab the name of the galaxies from a galaxy bundle.

    ...

    Parameters
    ----------
    filename : str, list
        The file name or the list of the galaxy bundle.


    Return
    -------
    storage
        A 1D numpy array of the apparant magnitude of set component.

    """
    if type(filename) == str:
        table = read_list(filename)
    elif type(filename) == list:
        table = filename
        
    storage = []
    for row in range(len(table)):
                storage.append(table[row][0])
                
    storage = np.array(storage)       
    return storage
#%%
def grab_name_specific(filename,keyword):
    """
    A function to grab the name of the galaxies from a galaxy bundle if it has 
    a specific keyword component

    Parameters
    ----------
    filename : str
        The file name of the galaxy bundle.
    keyword : str, list
        The name of the componets, such as: "Bulge", "Disk", "PrimBar", etc.

    Returns
    -------
    storage : list
        A 1D numpy array of the apparant magnitude of set component.


    """
    storage = []
    if type(filename) == str:
        table = read_list(filename)
    elif type(filename) == list:
        table = filename
        
    storage = []
    #generate a list of singular parameter
    for row in range(len(table)):
        for item in range(len(table[row])):
            if table[row][item] in (keyword):
                storage.append(table[row][0])
                
    storage = np.array(storage) 
    
    return storage

#%% tested
def grab_parameter(filename, keyword, number):
    """
    A function to grab a parameter of a component from a galaxy bundle.

    ...

    Parameters
    ----------
    filename : str
        The file name of the galaxy bundle.
    keyword: str, list
        The name of the componets, such as: "Bulge", "Disk", "PrimBar", etc.
    number: float
        The index of the function parameters

    Return
    -------
    storage
        A 1D numpy array of the apparant magnitude of set component.
    """
    if type(filename) == str:
        table = read_list(filename)
    elif type(filename) == list:
        table = filename
        
    storage = []
    #generate a list of singular parameter
    for row in range(len(table)):
        for item in range(len(table[row])):
            if table[row][item] in (keyword):
                storage.append(table[row][item+1][number])
                
    storage = np.array(storage)       
    return storage
#%%
def grab_parameter_whole(filename, keyword):
    """
    A function to grab the parameter set as a whole
    of a component from a galaxy bundle.

    ...

    Parameters
    ----------
    filename : str
        The file name of the galaxy bundle.
    keyword: str, list
        The name of the componets, such as: "Bulge", "Disk", "PrimBar", etc.
    number: float
        The index of the function parameters

    Return
    -------
    storage: 1D numpy array
        A 1D numpy array of the apparant magnitude of set component.
    """
    if type(filename) == str:
        table = read_list(filename)
    elif type(filename) == list:
        table = filename
        
    storage = []
    #generate a list of singular parameter
    for row in range(len(table)):
        for item in range(len(table[row])):
            if table[row][item] in (keyword):
                #print(table[row][item+1])
                storage.append(table[row][item+1])
                
    storage = np.array(storage)       
    return storage

#%% tested
def grab_mag(filename, keyword):
    """
    A function for grab the magnitude of a component from a galaxy bundle.

    ...

    Parameters
    ----------
    filename : str
        The file name of the galaxy bundle.
    keyword: str, list
        The name of the componets, such as: "Bulge", "Disk", "PrimBar", etc.

    Return
    -------
    storage
        A 1D numpy array of the apparant magnitude of set component.

    """
    if type(filename) == str:
        table = read_list(filename)
    elif type(filename) == list:
        table = filename
    
    storage = []
    for row in range(len(table)):
        for item in range(len(table[row])):
            if table[row][item] in (keyword):
                storage.append(table[row][item+2])
                
    storage = np.array(storage)       
    return storage

#%% tested
def grab_total_mag(filename):
    """
    A function for grabing the total magnitude of a 
    component from a galaxy bundle.

    ...

    Parameters
    ----------
    filename : str
        The file name of the galaxy bundle.

    Return
    -------
    storage
        A 1D numpy array of the total apparant magnitude.

    """
    if type(filename) == str:
        table = read_list(filename)
    elif type(filename) == list:
        table = filename
        
    storage = []
    for row in range(len(table)):
        for item in range(len(table[row])):
            if table[row][item] == 'Total_mag':
                storage.append(table[row][item+1][0])
                
    storage = np.array(storage)       
    return storage

#%% tested
def grab_dist(dir_dist_list,dist_list,
              name_index = 0, vel_index =3, vel_err_index = 4, scale_index = 5,
              dir_name_index = 0, dir_dist_index = 1, dir_dist_err_index =2,
              dir_method_index = 5, dir_method_flag_index = 7):
    """
    A function for grabing the distance of galaxy base on two lists.
    
    ...

    Parameters
    ----------
    dist_list: str
        The ASCII file name of a list with estimated velocity. 
        It is assume as default this list contains all the galaxy of interested 

    Optional
    --------
    dir_dist_list : str
        The ASCII file name of the redshift independent distance 
        measurement list. It override the dist_list distance estimation.
        This list is allowed to have different dimension as dist_list.
         
        
    Return
    -------
    
    A dictionary containing: 
        "Gal_name": the name of the galaxy, 
        "Dist": the distance in Mpc, 
        "Scale": the angular scale of kpc/arcsec
    
    """
    
    dir_dist_list1 ,dir_dist_list2 = np.genfromtxt(dir_dist_list, dtype='str'),\
    np.genfromtxt(dir_dist_list, dtype='float')
    
    dist_list1,dist_list2 = np.genfromtxt(dist_list, dtype='str'), \
    np.genfromtxt(dist_list, dtype='float')
    
    ## calculate v/H = dist
    ## replace dir_dist and calculate scale
    ## store in numpy array
    
    name = dist_list1[:,name_index]
    vel , vel_err= dist_list2[:,vel_index], dist_list2[:,vel_err_index]
    scale = dist_list2[:,scale_index]
    
    dist = vel / H0 
    dist_err = vel_err / H0
    
    #####
    dir_name = dir_dist_list1[:,dir_name_index]
    dir_dist =  dir_dist_list2[:,dir_dist_index]
    dir_dist_err =  dir_dist_list2[:, dir_dist_err_index]
    dir_method = dir_dist_list1[:,dir_method_index]
    dir_method_flag = dir_dist_list2[:,dir_method_flag_index]
    
    #####
    dir_scale = dir_dist* ((1e3) / 206264.806) #turn Mpc/rad to kpc/arcsec 
    
    ####replace####
    for row in range(len(name)):
        for row2 in range(len(dir_name)):
            if name[row] == dir_name[row2]:
                dist[row] = dir_dist[row2]
                dist_err[row] = dir_dist_err[row2]
                scale[row] = dir_scale[row2]
            else:
                pass
                
    return {"Gal_name": name, "Dist": dist,"Dist_err": dist_err, "Scale": scale}

#%% tested
def extract_match(keyword_list, match_list ,value_can):
    """
    Selecct a list of "values" if the keyword matches the element in the match 
    list.

    Parameters
    ----------
    keyword_list : TYPE
        DESCRIPTION.
    match_list : TYPE
        DESCRIPTION.
    value_can : TYPE
        DESCRIPTION.

    Returns
    -------
    value.

    """
    value = []
    for i in range(len(keyword_list)):
        #print(keyword_list[i])
        for j in range(len(match_list)):
            if keyword_list[i] == match_list[j]:
                #print(keyword_list[i],match_list[j])
                value.append(value_can[j])
    return value
#%% tested
def pd_read(filename,check_equvi):
    """
    A function for reading Profiler(by Bogdan) output log file, using pandas.
    It generate a machine readable list for further analysis.
    ...

    Parameters
    ----------
    filename : str
        The ASCII file name containing inoformation of one decomposition
    check_equvi: bool
        The indication if the input file is on equvalent axis. 
        check_equvi = TRUE if so.

    Return
    -------
    Gal_list
        A list containing the information of each components.
        e.g.: [NGC1234, 1.02,"Sersic", [3,12,1.5], 12, "Exp", [15,5],15, 
        "Total_mag",11] 
    """

    A = pd.read_table(filename,sep='\n')  #seperating the string by \n
    #print(filename)
    Gal_list=[]
    Gal_list.append(filename[20:27])  #recording the name of the galaxy
    #equvi = str(check_equvi)
    i=0 
    for i in range(A.shape[0]):
        
        # record the rms for the residual as the second element of the list
        if A.iloc[i][0] == "-----------------------------Detailed fit report:-------------------------------- ":
                Gal_list.append(float(A.iloc[i-1][0].split()[3]))
        
        # record the fitting parameters for the Sersic function 
        # ['Sersic', np.array([mu_e. r_e, n])], np.array([mag]])
        elif A.iloc[i][0][0:17] == "Sersic component:":
                
                mu_e = float(A.iloc[i+1][0].split()[2])
                r_e = float(A.iloc[i+2][0].split()[2])
                n = float(A.iloc[i+3][0].split()[2])
                Sersic_mag = float(A.iloc[i+4][0].split()[4]) if check_equvi else 0
                
                Gal_list.append('Sersic')
                Gal_list.append(np.array([mu_e,r_e,n]))
                Gal_list.append(Sersic_mag)
                
        # record the fitting parameters for the Core Sersic function 
        # ['CoreSersic', np.array([mu_p. r_e, r_b, n, alpha, gamma])], np.array([mag]])
        elif A.iloc[i][0][0:22] == "Core-Sersic component:":
            
                mu_p = float(A.iloc[i+1][0].split()[2])
                r_e = float(A.iloc[i+2][0].split()[2])
                n = float(A.iloc[i+4][0].split()[2])

                r_b = float(A.iloc[i+3][0].split()[2])
                alpha = float(A.iloc[i+5][0].split()[2])
                gamma = float(A.iloc[i+6][0].split()[2])
                
                CoreSersic_mag = float(A.iloc[i+7][0].split()[4]) if check_equvi else 0
                Gal_list.append('CoreSersic')
                Gal_list.append(np.array([mu_p,r_e,n,r_b,alpha,gamma]))
                Gal_list.append(CoreSersic_mag)
                
        # record the fitting parameters for the Exponential function 
        # ['Exp', np.array([mu_0. h])], np.array([mag]])
        elif A.iloc[i][0][0:22] == "Exponential component:":

                mu_0 = float(A.iloc[i+1][0].split()[2])
                h = float(A.iloc[i+2][0].split()[2])
                Exp_mag = float(A.iloc[i+3][0].split()[4]) if check_equvi else 0
            
                Gal_list.append('Exp')                
                Gal_list.append(np.array([mu_0,h]))
                Gal_list.append(Exp_mag)

        # record the fitting parameters for the Borken Exponential function 
        # ['BrokenExp', np.array([mu_0td. r_b , h1, h2])], np.array([BrokenExp_mag])
        elif A.iloc[i][0][0:25] == "Truncated disc component:":

                
                mu_0td = float(A.iloc[i+1][0].split()[2])
                r_b = float(A.iloc[i+2][0].split()[2])
                h1 = float(A.iloc[i+3][0].split()[2])
                h2 = float(A.iloc[i+4][0].split()[2])
                BrokenExp_mag = float(A.iloc[i+5][0].split()[4]) if check_equvi else 0

                Gal_list.append('BrokenExp')                                
                Gal_list.append(np.array([mu_0td,r_b,h1,h2]))
                Gal_list.append(BrokenExp_mag)
                
        # record the fitting parameters for the Ferrer function 
        # ['Ferrer', np.array([mu_0,r_0,alpha,beta])], np.array([Ferrer_mag]])
        elif A.iloc[i][0][0:17] == "Ferrer component:":
            
                mu_0 = float(A.iloc[i+1][0].split()[2])
                r_0 = float(A.iloc[i+2][0].split()[2])
                alpha = float(A.iloc[i+3][0].split()[2])
                beta = float(A.iloc[i+4][0].split()[2])
                Ferrer_mag = float(A.iloc[i+6][0].split()[4]) if check_equvi else 0
                
                Gal_list.append('Ferrer')                                
                Gal_list.append(np.array([mu_0,r_0,alpha,beta]))
                Gal_list.append(Ferrer_mag)
                
        # record the fitting parameters for the Gaussian function 
        # ['InclExp', np.array([mu_0z. z0])], np.array([InclExp_mag]])                
        elif A.iloc[i][0][0:24] == "Inclined disc component:":

                mu_0z = float(A.iloc[i+1][0].split()[2])
                z_0 = float(A.iloc[i+2][0].split()[2])
                InclExp_mag = float(A.iloc[i+3][0].split()[4]) if check_equvi else 0
                
                Gal_list.append('InclExp')                                
                Gal_list.append(np.array([mu_0z,z_0]))
                Gal_list.append(InclExp_mag)
                
                print(filename[2:9])
                
        # record the fitting parameters for the Gaussian function 
        # ['Gauss', np.array([mu_0. r_0, FWHM])], np.array([Gauss_mag]])
        elif A.iloc[i][0][0:19] == "Gaussian component:":

                mu_0 = float(A.iloc[i+1][0].split()[2])
                r_0 = float(A.iloc[i+2][0].split()[2])
                FWHM = float(A.iloc[i+3][0].split()[2])          
                Gauss_mag = float(A.iloc[i+4][0].split()[4]) if check_equvi else 0
                
                Gal_list.append('Gauss')                                
                Gal_list.append(np.array([mu_0,r_0,FWHM]))
                Gal_list.append(Gauss_mag)
                
        # record the fitting parameters for the Gaussian function 
        # ['PSF', np.array([mu_0])], np.array([PSF_mag]])              
        elif A.iloc[i][0][0:24] == "PSF [nuclear] component:":

                mu_0 = float(A.iloc[i+1][0].split()[2])
                PSF_mag = float(A.iloc[i+2][0].split()[4]) if check_equvi else 0
                
                Gal_list.append('PSF')                                                
                Gal_list.append(np.array([mu_0]))
                Gal_list.append(PSF_mag)
                
        # record the total magnitude of the galaxy 
        # ['Total_mag',  total_mag]
        elif A.iloc[i][0][0:23] == "Total galaxy magnitude:":

                total_mag = float(A.iloc[i][0].split()[3])
                
                Gal_list.append('Total_mag')                                                
                Gal_list.append([total_mag])            
    return Gal_list

#%% tested
def run_list(input_list,output_name,check_equvi):
    """
    A function for reading Profiler(by Bogdan) output log file, using pandas.
    It generate a machine readable list for further analysis.
    ...

    Parameters
    ----------
    input_list : str
        The ASCII file with each row pointing to one galaxy log file location.
        e.g.: ./my_gal/NGC0001/NGC0001_major.txt
    output_name: str
        The name of the output pickle file
    check_equvi: bool
        The indication if the input file is on equvalent axis. 
        check_equvi = TRUE if so.


    Return
    -------
    Gal_bundle
        A list of lists containing the information of each components.
        e.g.: [[NGC0001,...],[NGC0002,...],[NGC0003,...]]
    """

    input_gal = pd.read_table(input_list,sep='\n')
    Gal_bundle = []
    j=0 
    for j in range(input_gal.shape[0]): #read individual galaxy
    
        print(input_gal.iloc[j][0])
        Gal_list = pd_read(input_gal.iloc[j][0],check_equvi) 
        Gal_bundle.append(Gal_list)
        #print(Gal_list)
        
    with open(output_name, 'wb') as f: # saving the list 
        pickle.dump(Gal_bundle, f)
        
    return Gal_bundle
#%%
