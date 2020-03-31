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
import matplotlib.pyplot as plt
from numpy.linalg import inv

__all__ = ["read_list", "read_table", "match_list_dim", "lookup_bundle",
           "grab_parameter","grab_mag","grab_total_mag", "grab_dist", 
           "grab_info_mag","pd_read","run_list"]


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
def read_table(table):
    """
    Just a function using genfromtext to load objects
    assuming the data is float
    """
    D=np.genfromtxt(table, dtype='float')
    return D
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

#%% tested
def lookup_bundle(gal_bundle,gal_name):
    
    """
    Lookup individual galaxy by name in a galaxy bundle
    
    Parameters
    ----------
    gal_bundle: 
        
    gal_name:
        
    Return
    ----------
    
    Generates..functiona blah blah
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
def grab_parameter(filename, keyword, number):
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
    
    table = read_list(filename)
    storage = []
    #generate a list of singular parameter
    for row in range(len(table)):
        for item in range(len(table[row])):
            if table[row][item] in (keyword):
                storage.append(table[row][item+1][number])
                
    storage = np.array(storage)       
    return storage
#%%
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
    table = read_list(filename)
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
    table = read_list(filename)
    storage = []
    for row in range(len(table)):
        for item in range(len(table[row])):
            if table[row][item] == 'Total_mag':
                storage.append(table[row][item+1][0])
                
    storage = np.array(storage)       
    return storage
#%%
def grab_dist(dir_dist_list=None,dist_list):
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

def grab_info_mag():
    return None
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
    Gal_list.append(filename[2:9])  #recording the name of the galaxy
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
                r_b = float(A.iloc[i+3][0].split()[2])
                n = float(A.iloc[i+4][0].split()[2])
                alpha = float(A.iloc[i+5][0].split()[2])
                gamma = float(A.iloc[i+6][0].split()[2])
                CoreSersic_mag = float(A.iloc[i+7][0].split()[4]) if check_equvi else 0

                Gal_list.append('CoreSersic')
                Gal_list.append(np.array([mu_p,r_e,r_b,n,alpha,gamma]))
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
        Gal_list = pd_read(input_gal.iloc[j][0],check_equvi) 
        Gal_bundle.append(Gal_list)
        #print(Gal_list)
        
    with open(output_name, 'wb') as f: # saving the list 
        pickle.dump(Gal_bundle, f)
        
    return Gal_bundle
#%%

#%%
#%% require heavy modification or deletion, haven't decided yet
def create_table(Dtable,Ctable,M_sun,N):
    D1 = np.genfromtxt(Dtable , dtype='str')
    D2 = np.genfromtxt(Dtable , dtype='float')

    C1 = np.genfromtxt(Ctable , dtype='str')
    C2 = np.genfromtxt(Ctable , dtype='float')

    name = D1[:,0]
    n_maj = D2[:,1]
    Re_maj = D2[:,2]
    h_maj = D2[:,3]
    n_equ = D2[:,4]
    Re_equ = D2[:,5]
    h_equ = D2[:,6]
    mag_bulge = D2[:,7]
    mag_total = D2[:,8]

    zdist = C2[:,1]
    zdist_err = C2[:,2]
    AbsMag_F = C2[:,3]
    AbsMag_N = C2[:,4]
    AbsMag_u = C2[:,5]
    AbsMag_g = C2[:,6]
    AbsMag_r = C2[:,7]
    AbsMag_i = C2[:,8]
    AbsMag_z = C2[:,9]
    morph = C1[:,10]
    
    print(C2[:,12])
    maj_max = C2[:,12]
    equ_max = C2[:,13]

    scale, Re_kpc_maj, Re_kpc_equ, Dist = [], [], [], []

    for i in range(N):
        #print Re_maj[i]
        NWS=ned_wright_CosCal(zdist[i],0.3,0.7)
        scale = NWS.kpc_DA
        #print(scale)
        #print scale*Re_maj[i]
        #print Re_kpc_maj
        Re_kpc_maj.append(scale* Re_maj[i])    #effective radius
        Re_kpc_equ.append(scale* Re_equ[i])
        Dist.append(NWS.DCMR_Mpc)

    scale, Re_kpc_maj, Re_kpc_equ, Dist = np.array(scale), np.array(Re_kpc_maj), np.array(Re_kpc_equ), np.array(Dist)
    #print(name)
    #print(Dist)

    Mass_bulge_Into = cal_Mass(mag_bulge,Dist,ML_relation_Iband(AbsMag_g,AbsMag_i).Into13_MassRatio,M_sun)
    Mass_total_Into = cal_Mass(mag_total,Dist,ML_relation_Iband(AbsMag_g,AbsMag_i).Into13_MassRatio,M_sun)

    Mass_bulge_Taylor = cal_Mass(mag_bulge,Dist,ML_relation_Iband(AbsMag_g,AbsMag_i).Taylor11_MassRatio,M_sun)
    Mass_total_Taylor = cal_Mass(mag_total,Dist,ML_relation_Iband(AbsMag_g,AbsMag_i).Taylor11_MassRatio,M_sun)
    return tuple([Re_kpc_maj, Re_kpc_equ, Mass_bulge_Into, Mass_total_Into, Mass_bulge_Taylor, Mass_total_Taylor, name, n_equ, zdist, maj_max, equ_max,Dist])

