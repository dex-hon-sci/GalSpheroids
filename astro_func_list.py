#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 31 20:30:19 2020

@author: dexter
"""


import numpy as np

def Schechter_func(m, alpha, phi_0): # mass funciton formulation
    """
    
    """
    return np.log(10)*phi_0*10**(np.log10(m)*(alpha+1))*np.exp(-10**np.log10(m))