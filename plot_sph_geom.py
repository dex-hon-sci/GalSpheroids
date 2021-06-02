#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 12:14:07 2021

@author: dexter

This is a script dedicated to
visualise the spheroid parameters.

This script produce the following plots:
1) mu_0 - n plots
2) mu_0 - Re plots
3) size-mass plot with curved fit

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

## Read in the data

# Read in my data and the others.


# Use the Sersic function to calculate mu_0



