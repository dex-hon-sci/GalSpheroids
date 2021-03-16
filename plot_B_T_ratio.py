#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 21:10:56 2021

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



bundle_name = "/home/dexter/SphProject/F_Gal_bundle_equvi_V_cpt"

morph_file = "/home/dexter/result/stat/completeness/morph_list.txt"

#morph_source = SRead.read_table(morph_file)
morph_source_t = SRead.read_table(morph_file,dtype="str")

morph_name, morph = morph_source_t[:,0], morph_source_t[:,1]

fig = plt.figure()

gs = gridspec.GridSpec(1,1) 

ax = plt.subplot(gs[0])

S = SPlot.ShowcaseIndi.cpt_to_total_by_type_plot(bundle_name, morph_name, 
                                                 morph, AX =ax)

