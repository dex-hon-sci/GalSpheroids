#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 20:50:25 2020

@author: dexter
"""
#%%
"""
import packages
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


__all__= []

__author__="Dexter S.H. Hon"


#%%
class Image():
    """
    A class containing the basic functions to manipulate the image
    
    Attributes
    ----------
        

    Methods
    -------
    
    ds9 example: wcs, marking
    
    Preprocessing, 
    
    
    
    detect seeing, 
    find sky value, -> plot 
    factorization
    
    image rotation
    image math, subtraction and division
    image mirror flip
    
    masking, and polygon
    
    scaling class or unsharp mask
    
    Dithering
    
    
    Source dectection (preselection process)
    
    GUI?
    
    """
    
    
    
