#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 01:11:25 2020

@author: dexter
"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import scipy as scipy
from scipy import stats

from astropy.table import Table, Column, MaskedColumn
import matplotlib.patches as mpatches


import SphRead as SRead
import SphPlot as SPlot

class AnalyticFunctions(object):
    """
    """
    def __init__(self,r=None,theta = None):
        self._r = r
        self.theta = theta
        
    @property
    def r(self):
        return(self._r)
        
    def sersic(mu_e,r_e,n):
        return None
    
    def exp():
        return None
    
    def broken_exp():
        return None


    def size_mass_powerlaw(M,a,b):
        """
        R_e = a*(M_*/(1e10*M_{\odot}))**b

        Returns
        -------
        Effective radius.

        """
        return a*(M/1e10)**b
    
    def size_mass_twopowerlaw(M,a,b,c,Mo,M_sun):
        """
        R_e = c*(M_*/(M_{\odot})*(1+(M/Mo))**(b-a)

        Returns
        -------
        Effective radius.

        """
        return c*(M/M_sun)*(1+(M/Mo))**(b-a)
    
    
#%%
from scipy.special import gamma
from scipy.special import gammainc

def b_value(n,x):
   bn=[]
   diff=[]
   b= np.linspace(0.1, 50.0, 50000) # works for most z
#   x=1/z=10 # 10, 2, 10/9    
#   b= np.linspace(0.1, 25.0, 2500)  #Re  2
#   b= np.linspace(0.1, 50.0, 25000)  #R10  10
#   b= np.linspace(0.1, 100.0, 50000)  #R90  10/9
   for j in b:        
        g2= x*gammainc(2*n,j)
        g22=np.round(g2,2)
        if g22==1.00: ## gammainc(2*n,j) is normalized by dividing with gamma(2*n)
            k=np.round(j,4)
            bn.append(k)
            dif=abs(1-g2)
            diff.append(dif)
   
   diff_min= min(diff)
   for s, d in zip(bn, diff):
       if d==diff_min:
           b_final=s
       else:
           continue     
   
   return b_final    
#%%    
class PlotProfile(AnalyticFunctions):
    
    def __init__(self,*args, **kwargs):
        super().__init__(*args, **kwargs)
#%%
class Fitting():
    def __init__(self,*args, **kwargs):
        super().__init__(*args, **kwargs)
#class 2Dplot(object)

#%%
def isofit_isophote():
    return None

def radius_by_percentage(input_file,max_radius=None,percentage=None):
    """
    A function to calculate the radius of the galaxy at some percentage of 
    light. 
    If one 
    
    Parameters
    ----------
    input_file: str
        The output file of ISOFIT
    max_radius: float
        The maximum radius of the galaxy, given in arcsec.
    percentage: float
        The percentage of light, from the centre extending outward. 

    Returns
    -------
    perc_radius: float
        The radius of the galaxy in which the amount of light reaches the 
        assigned percentage parameter.

    """
    #Read inputfile, b/a e and so on
    
    #calculate the light within each isophote
    
    #calculate the total light covered within max_radius
    
    #loop from the centre, and numerically determine the radius in which 
    # x percentage of light are included (stopping condition)
    
    
    
    return None