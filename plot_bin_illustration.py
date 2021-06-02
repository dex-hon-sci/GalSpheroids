#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 15:29:42 2021

@author: dexter
"""
import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.style
import matplotlib as mpl

import matplotlib.gridspec as gridspec

plt.style.use('classic')

mpl.rcParams['grid.linewidth'] = 1.0

def plot_bin_illustration():
    
    plt.hlines(3, 2, 3)
    plt.hlines(2, 1, 2)
    plt.hlines(1, 0, 1)
 
    plt.vlines(3, 3, 6)
    plt.vlines(2, 2, 6)
    plt.vlines(1, 1, 6)    

    plt.hlines(1, 1, 3, linestyles='dashed')
    plt.hlines(2, 0, 2, linestyles='dashed',color='k')
    plt.hlines(3, 0, 3, linestyles='dashed',color='k')
    
    plt.vlines(3, 1, 6, linestyles='dashed')
    plt.vlines(2, 1, 6, linestyles='dashed')      
    
    plt.fill([0,1,1,0],[1,1,6,6], color='g', alpha =0.3)
    plt.fill([1,2,2,1],[2,2,6,6], color='g', alpha =0.2)
    plt.fill([3,2,2,3],[3,3,6,6], color='g', alpha =0.1)
    
    plt.text(0.4, 3.5, r"$\rm N_3$", fontsize = 20)
    plt.text(1.4, 3.5, r"$\rm N_2$", fontsize = 20)
    plt.text(2.4, 3.5, r"$\rm N_1$", fontsize = 20)

    plt.text(1.4, 1.5, r"$\rm N_2'$", fontsize = 20)
    plt.text(2.4, 1.5, r"$\rm N_1'$", fontsize = 20)
   
    plt.yticks([1,2,3],[r"$\rm M_3$", r"$\rm M_2$", r"$\rm M_1$"],fontsize=16) 
    plt.xticks([1,2,3],[r"$\rm V_3$", r"$\rm V_2$", r"$\rm V_1$"],fontsize=16) 

    plt.xlim(0,3.2)
    plt.ylim(0,5)
    
    plt.ylabel(r"$ M_*/\rm M_\odot$", fontsize=22)
    plt.xlabel(r"$ Volume (V)$", fontsize=22)

    return None

plot_bin_illustration()