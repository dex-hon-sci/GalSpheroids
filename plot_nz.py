#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 22:05:52 2020

@author: dexter
"""

from astropy.cosmology import FlatLambdaCDM
from matplotlib.colors import LogNorm
import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import matplotlib.pyplot as plt
import numpy as np


import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')

mpl.rcParams['grid.linewidth'] = 1.0

#####

D = np.array([45,75,110])

voll = ((D**3)/3)*((214-139)/180)*np.pi*(np.cos((np.pi/2)-(55*np.pi/180)-np.cos(np.pi/2)))

V1,V2,V3=voll[2],voll[1],voll[0]



Barro_data = SRead.read_table("/home/dexter/result/stat/completeness/nd_Barro.dat")
Barro_z = Barro_data[:,0]
Barro_nd = Barro_data[:,1]
Barro_nd_uerr = Barro_data[:,3]

vDokkum_data = SRead.read_table("/home/dexter/result/stat/completeness/nd_vDokkum.dat")
vDokkum_z = vDokkum_data[:,0]
vDokkum_nd = vDokkum_data[:,1]
vDokkum_nd_uerr =  vDokkum_data[:,2]

vdWel_data = SRead.read_table("/home/dexter/result/stat/completeness/nd_vdWel.dat")
vdWel_z = vdWel_data[:,0]
vdWel_nd = vdWel_data[:,3]
vdWel_nd_uerr = vdWel_data[:,4]

my_z = np.array([0.0127,0.0085,0.005])

###

#my_nd_Barro = np.array([6.75e-4,1.19e-4,2.94e-5])
#my_nd_Dam = np.array([7.68e-4,1.39e-4,2.94e-5])
#my_nd_vDokkum = np.array([2.14e-4, 1.99e-5, 6.31e-6])
#my_nd_vdWel = np.array([2.76e-4,2.65e-5,6.31e-6])
#
#my_nd_Graham = np.array([1.54e-4,1.33e-5,0])
#
####
my_nd_Barro = np.array([4.41e-5,1.19e-4,4.3e-4])
my_nd_Dam = np.array([4.83e-5,1.33e-4,4.3e-4])
my_nd_vDokkum = np.array([8.41e-6, 1.33e-5, 9.21e-5])
my_nd_vdWel = np.array([1.89e-5,2.65e-5,9.21e-5])

my_nd_Graham = np.array([8.41e-6,6.63e-6,0])


ms0 = 12


xlim = [-0.35,3.2]
ylim = [1e-6,1e-3]

def plot_nd_3bins(nd,marker,AX=plt):
    n_err_bin1 = np.sqrt(nd[0]*V1)/V1
    n_err_bin2 = np.sqrt(nd[1]*V2)/V2
    n_err_bin3 = np.sqrt(nd[2]*V3)/V3
    
    xedge1, yedge1 = [-5,5,5,-5],[nd[0]-n_err_bin1, 
                                  nd[0]-n_err_bin1, 
                                  nd[0]+n_err_bin1, nd[0]+n_err_bin1]
    xedge2, yedge2 = [-5,5,5,-5],[nd[1]-n_err_bin2, 
                                  nd[1]-n_err_bin2, 
                                  nd[1]+n_err_bin2, nd[1]+n_err_bin2]
    xedge3, yedge3 = [-5,5,5,-5],[nd[2]-n_err_bin3, 
                                  nd[2]-n_err_bin3, 
                                  nd[2]+n_err_bin3, nd[2]+n_err_bin3]
    
    AX.fill(xedge1, yedge1,alpha=0.3,color='#a5200b')
    AX.fill(xedge2, yedge2,alpha=0.3,color='#0b5786')
    AX.fill(xedge3, yedge3,alpha=0.3,color='#2a3236')

    AX.axhline(nd[0],xlim[0],xlim[1],linestyle="dashed",color='#a5200b')
    AX.axhline(nd[1],xlim[0],xlim[1],linestyle="dashed",color='#0b5786')
    AX.axhline(nd[2],xlim[0],xlim[1],linestyle="dashed",color='#2a3236')

    AX.axhline((nd[0]+nd[1]+nd[2]),xlim[0],xlim[1],linestyle="dashed", linewidth = 6,
               color='black')


    AX.errorbar(my_z[0],nd[0],yerr=n_err_bin1,ls='none',
                     linewidth=3, ecolor='#a5200b',zorder=20,mew=1,capsize=3)
    AX.errorbar(my_z[1],nd[1],yerr=n_err_bin2,ls='none',
                     linewidth=3, ecolor='#0b5786',zorder=20,mew=1,capsize=3)
    AX.errorbar(my_z[2],nd[2],yerr=n_err_bin3,ls='none',
                     linewidth=3, ecolor='#2a3236',zorder=20,mew=1,capsize=3) 
    

    AX.plot(my_z[0],nd[0], color ='#a5200b', label="Bin1",ms= ms0,marker=marker)
    AX.plot(my_z[1],nd[1], color ='#0b5786', label="Bin2",ms= ms0,marker=marker)
    AX.plot(my_z[2],nd[2], color ='#2a3236', label="Bin3",ms= ms0,marker=marker)



def plot_nd_Barro(AX):
    
    plot_nd_3bins(my_nd_Barro,'o',AX=AX)

    AX.plot(Barro_z,Barro_nd,'--o', ms= ms0, color='g',label="Barro et al. 2013")
    AX.errorbar(Barro_z,Barro_nd,yerr=Barro_nd_uerr-Barro_nd,ls='none',
                     linewidth=3, ecolor='g',zorder=20,mew=1,capsize=3) 


    xedge, yedge = list(Barro_z), list(Barro_nd)
    xedge.append(Barro_z[-1])
    yedge.append(0)
    xedge.append(Barro_z[0])
    yedge.append(0)

    #AX.fill(xedge, yedge, alpha=0.1, color='g')
    #AX.set_xlabel(r"$ z$",fontsize=18)
    AX.set_ylabel(r"$n (Mpc^{-3}$)",fontsize= 18)


    AX.set_ylim(ylim[0],ylim[1])
    AX.set_xlim(xlim[0],xlim[1])

    AX.set_yscale( 'log' )
    #AX.legend()
            

def plot_nd_vDokkum(AX):
    
    plot_nd_3bins(my_nd_vDokkum,'^',AX=AX)


    AX.plot(vDokkum_z,vDokkum_nd,'--^', ms= ms0,color='#e58b1a', label="van Dokkum et al. 2015")
    AX.errorbar(vDokkum_z,vDokkum_nd,yerr= vDokkum_nd_uerr - vDokkum_nd,ls='none',
                     linewidth=3, ecolor='#e58b1a',zorder=20,mew=1,capsize=3) 


    xedge, yedge = list(vDokkum_z), list(vDokkum_nd)
    xedge.append(vDokkum_z[-1])
    yedge.append(0)
    xedge.append(vDokkum_z[0])
    yedge.append(0)

    #AX.fill(xedge, yedge, alpha=0.1, color='y')
    #AX.set_xlabel(r"$ z$",fontsize=18)
    AX.set_ylabel(r"$ n (Mpc^{-3})$",fontsize= 18)
    AX.set_yscale( 'log' )
    #AX.legend()
        
    AX.set_ylim(ylim[0],ylim[1])
    AX.set_xlim(xlim[0],xlim[1])
    
def plot_nd_vdWel(AX):
    
    plot_nd_3bins(my_nd_vdWel,'s',AX=AX)


    AX.plot(vdWel_z,vdWel_nd,'--s', ms= ms0, color='b', label="van der Wel et al. 2014")
    AX.errorbar(vdWel_z,vdWel_nd,yerr= vdWel_nd_uerr - vdWel_nd,ls='none',
                     linewidth=3, ecolor='b',zorder=20,mew=1,capsize=3) 

    xedge, yedge = list(vdWel_z), list(vdWel_nd)
    xedge.append(vdWel_z[-1])
    yedge.append(0)
    xedge.append(vdWel_z[0])
    yedge.append(0)

    #AX.fill(xedge, yedge, alpha=0.1, color='b')
    #AX.set_xlabel(r"$ z$",fontsize=18)
    AX.set_ylabel(r"$ n (Mpc^{-3})$",fontsize= 18)
    AX.set_yscale( 'log' )
    #AX.legend()
    
    AX.set_ylim(ylim[0],ylim[1])
    AX.set_xlim(xlim[0],xlim[1])        
#############################
#
import matplotlib.gridspec as gridspec

def plot_nd_3plot():
    
    fig = plt.figure()
    gs = gridspec.GridSpec(ncols=1, nrows=3,
                               hspace=0.1, wspace=0.0) 

    #plot Panel (1)
    axs0 = plt.subplot(gs[0])  
    plot_nd_Barro(axs0)
    plt.setp(axs0.get_xticklabels(), visible=False)
    axs0.grid(True)

    
    #plot Panel (2)
    axs1 = plt.subplot(gs[1],sharex=axs0)  
    plot_nd_vDokkum(axs1)
    plt.setp(axs1.get_xticklabels(), visible=False)
    axs1.grid(True)

    
    #plot Panel (3)
    axs2 = plt.subplot(gs[2],sharex=axs0)      
    plot_nd_vdWel(axs2)
    axs2.set_xlabel(r"$ z$",fontsize=18)
    axs2.grid(True)
    plt.show()
    
plot_nd_3plot()