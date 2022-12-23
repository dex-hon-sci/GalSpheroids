#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 16 20:19:27 2021

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

file_list=""

file = SRead.read_table(file_list,dtype="str")
for i in range(len(file)):
    
    SPlot.Plot2D.find_mode_sky(file_list[i])