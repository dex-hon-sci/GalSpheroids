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

# seperate the data into two subset
n = np.array([[1,2,3],[4,5,6]])
mu0 = np.array([[1,2,3], [23,6,1]])
Mag = np.array([[1,2,3],[21,5,8]])


def plot_n_mu0_Mag_2plot(n,mu0,Mag,label=[],fit_instruc=0):
    """
    Plot a similar plot as in Graham 2019 figure 1 with my sample

    Parameters
    ----------
    n : numpy array
        Sersic indices
    mu0 : numpy array
        The central surface brightness.
    Mag : numpy array
        The absoulte magnitude.
    label : list, optional
        optional list of labels. The default is [].
    fit_instruc: float
        Fitting instruction. The default is the 0 element of the input array.
    """
    #Fitting the relevant data
    
    #Plotting
    fig = plt.figure()
    gs = gridspec.GridSpec(ncols=2, nrows=1,
                               hspace=0, wspace=0.0) 

    # The Mag vs Sersic n plot
    axt0 = plt.subplot(gs[0])
    
    # Check the dimension of the input array 
    if n.ndim == 1:
        axt0.plot(n,Mag,'o',ms=10)
    elif n.ndim > 1:
        print('yes')
        for i in range(len(n)):
            axt0.plot(n[i],Mag[i],'o',ms=10,label=label[i])
    
    axt0.legend()
    axt0.set_ylabel(r"$\mathfrak{M}_{i}$", fontsize=22)
    axt0.set_xlabel(r"$\mathrm{Sersic}~n $", fontsize=22)
    
    # The Mag vs mu_0 plot
    axt1 = plt.subplot(gs[1],sharey=axt0)
    
    # Check the dimension of the input array  
    if mu0.ndim == 1:
        axt1.plot(mu0,Mag,'o',ms=10)
    elif mu0.ndim > 1:
        print('yes')
        for i in range(len(mu0)):
            axt1.plot(mu0[i],Mag[i],'o',ms=10,label=label[i])    
    
    axt1.legend()
    axt1.set_xlabel(r"$ \mu_\mathrm{0,i}$", fontsize=22)
    plt.setp(axt1.get_yticklabels(), visible=False)
    plt.gca().invert_yaxis()
    plt.show()
    return fig

# plot the Mag vs n and mu0 plot
plot_n_mu0_Mag_2plot(n,mu0,Mag,label=[r"$type~1$",r"$type~2$"])
