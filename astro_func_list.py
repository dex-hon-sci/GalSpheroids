#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 31 20:30:19 2020

@author: dexter
"""


import numpy as np

def Schechter_func(M, alpha, M_star, phi_0): # mass funciton formulation
    """
    The Schechter function that describe a mass function.
    """
    return np.log(10)*phi_0*10**(np.log10(M/M_star)*(alpha+1))*np.exp(
        -10**np.log10(M/M_star))

