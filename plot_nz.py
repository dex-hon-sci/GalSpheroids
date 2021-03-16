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

#V1,V2,V3=voll[2],voll[1],voll[0]

V1_V = voll[2]-voll[1]
V2_V = voll[1]-voll[0]
V3_V = voll[0]

V1,V2,V3 = V1_V,V2_V,V3_V

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

Dam_data = SRead.read_table("/home/dexter/result/stat/completeness/nd_Damjanov15.dat")
Dam_z = Dam_data[:,0]
Dam_nd = 10**Dam_data[:,1]
Dam_nd_uerr = 10**Dam_data[:,2]


my_z = np.array([0.0127,0.0085,0.005])

###

#my_nd_Barro = np.array([6.75e-4,1.19e-4,2.94e-5])
#my_nd_Dam = np.array([7.68e-4,1.39e-4,2.94e-5])
#my_nd_vDokkum = np.array([2.14e-4, 1.99e-5, 6.31e-6])
#my_nd_vdWel = np.array([2.76e-4,2.65e-5,6.31e-6])
#
#my_nd_Graham = np.array([1.54e-4,1.33e-5,0])
#
#####
#Bin number density
#my_nd_Barro = np.array([4.41e-5,1.19e-4,4.3e-4])
#my_nd_Dam = np.array([4.83e-5,1.33e-4,4.3e-4])
#my_nd_vDokkum = np.array([8.41e-6, 1.33e-5, 9.21e-5])
#my_nd_vdWel = np.array([1.89e-5,2.65e-5,9.21e-5])
#
#my_nd_Graham = np.array([8.41e-6,6.63e-6,0])

#####
#Bin V number density
my_nd_Barro = np.array([4.92e-5,1.10e-4,7.37e-4])
my_nd_Dam = np.array([4.92e-5,1.27e-4,7.98e-4])
my_nd_vDokkum = np.array([9.23e-6, 1.69e-5, 1.23e-4])
my_nd_vdWel = np.array([6.15e-6,2.54e-5,1.23e-4])

my_nd_Graham = np.array([1.23e-6,8.46e-6,0])

#my_nd_E = np.array([6/V1, 3/V2, 2/V3])#8.97e-5]
my_nd_E = np.array([6/V1, 3/V2, 3/V3])#8.97e-5] BinV

my_nd_E_Bin1 = np.array([my_nd_E[0],0])
my_nd_E_Bin2 = np.array([my_nd_E[1],0])
my_nd_E_Bin3 = np.array([my_nd_E[2],0])


#my_nd_oldE = np.array([15/V1, 11/V2, 2/V3])
my_nd_oldE = np.array([13/V1, 9/V2, 6/V3]) #BinV

my_nd_oldE_Bin1 = np.array([my_nd_oldE[0],0])
my_nd_oldE_Bin2 = np.array([my_nd_oldE[1],0])
my_nd_oldE_Bin3 = np.array([my_nd_oldE[2],0])

ms0 = 12

print(np.sum(my_nd_Barro), np.sum(my_nd_vdWel), np.sum(my_nd_vDokkum))

print("cE to E ratio:",
      np.sum(my_nd_Barro)/np.sum(my_nd_E), 
      np.sum(my_nd_vdWel)/np.sum(my_nd_E) , 
      np.sum(my_nd_vDokkum)/np.sum(my_nd_E))

print("E-peak RN /peak RN:",
      (np.sum(my_nd_E)/max(Barro_nd)), 
      (np.sum(my_nd_E)/max(vdWel_nd)), 
      (np.sum(my_nd_E)/max(vDokkum_nd)))

print("csph-peak RN /peak RN:",
      (np.sum(my_nd_Barro)/max(Barro_nd)), 
      (np.sum(my_nd_vdWel)/max(vdWel_nd)), 
      (np.sum(my_nd_vDokkum)/max(vDokkum_nd)))

xlim = [-0.35,3.0]
ylim = [1e-6,1e-3]

def plot_compact_sum(nd0, linestyle = "dashed", colour="black", label = "", AX=plt):
    """
    plot the sum of compact sample and the Ellipitcals

    Returns
    -------
    None.

    """ 
 
    AX.axhline(sum(nd0),xlim[0],xlim[1],linestyle=linestyle, linewidth = 5,
               color=colour , label= label)
    

    

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
    
    #AX.fill(xedge1, yedge1,alpha=0.3,color='#a5200b')
    #AX.fill(xedge2, yedge2,alpha=0.3,color='#0b5786')
    #AX.fill(xedge3, yedge3,alpha=0.3,color='#2a3236')

    #AX.axhline(nd[0],xlim[0],xlim[1],linestyle="dashed",color='#a5200b')
    #AX.axhline(nd[1],xlim[0],xlim[1],linestyle="dashed",color='#0b5786')
    #AX.axhline(nd[2],xlim[0],xlim[1],linestyle="dashed",color='#2a3236')



    AX.errorbar(my_z[0],nd[0],yerr=n_err_bin1,ls='none',
                     linewidth=3, ecolor='#a5200b',mew=1,capsize=3)
    AX.errorbar(my_z[1],nd[1],yerr=n_err_bin2,ls='none',
                     linewidth=3,  ecolor='#0b5786',mew=1,capsize=3)
    AX.errorbar(my_z[2],nd[2],yerr=n_err_bin3,ls='none',
                     linewidth=3,  ecolor='#2a3236',mew=1,capsize=3) 
    

    AX.plot(my_z[0],nd[0], color ='#a5200b',ms= ms0,marker=marker)
    AX.plot(my_z[1],nd[1], color ='#0b5786',ms= ms0,marker=marker)
    AX.plot(my_z[2],nd[2], color ='#2a3236',ms= ms0,marker=marker)

def plot_nd_Dam(AX):
    AX.plot(Dam_z,Dam_nd,'--o', ms= ms0, color='r',label="Damjanov et al. 2015")
    AX.errorbar(Dam_z,Dam_nd,yerr=Dam_nd_uerr-Dam_nd,ls='none',
                     linewidth=3, ecolor='r',mew=1,capsize=3) 


    xedge, yedge = list(Dam_z), list(Dam_nd)
    xedge.append(Dam_z[-1])
    yedge.append(0)
    xedge.append(Dam_z[0])
    yedge.append(0)

    #AX.fill(xedge, yedge, alpha=0.1, color='g')
    #AX.set_xlabel(r"$ z$",fontsize=18)
    
    AX.set_ylabel(r"$n (Mpc^{-3}$)",fontsize= 18)


    AX.set_ylim(ylim[0],ylim[1])
    AX.set_xlim(xlim[0],xlim[1])

    AX.set_yscale( 'log' )
    #AX.legend() 

def plot_nd_Barro(AX):

    AX.plot(Barro_z,Barro_nd,'--o', ms= ms0, color='g',label="Barro et al. 2013")
    AX.errorbar(Barro_z,Barro_nd,yerr=Barro_nd_uerr-Barro_nd,ls='none',
                     linewidth=3, ecolor='g',mew=1,capsize=3) 


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

    AX.plot(vDokkum_z,vDokkum_nd,'--^', ms= ms0,color='#e58b1a', label="van Dokkum et al. 2015")
    AX.errorbar(vDokkum_z,vDokkum_nd,yerr= vDokkum_nd_uerr - vDokkum_nd,ls='none',
                     linewidth=3, ecolor='#e58b1a',mew=1,capsize=3) 


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


    AX.plot(vdWel_z,vdWel_nd,'--s', ms= ms0, color='b', label="van der Wel et al. 2014")
    AX.errorbar(vdWel_z,vdWel_nd,yerr= vdWel_nd_uerr - vdWel_nd,ls='none',
                     linewidth=3, ecolor='b',mew=1,capsize=3) 

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


def plot_nd_E_3plot():
        
    fig = plt.figure()
    gs = gridspec.GridSpec(ncols=1, nrows=4,
                               hspace=0.1, wspace=0.0) 

    #plot Panel (1)
    axs0 = plt.subplot(gs[0])  
    
    plot_compact_sum(my_nd_oldE_Bin1, linestyle="dashed", colour='blue',label = r"old $n_{E}$ in Bin 1", AX=axs0)
    plot_compact_sum(my_nd_E_Bin1, linestyle = "dashed", colour='red',label = r"new $n_{E}$ in Bin 1", AX=axs0)

    plot_compact_sum(my_nd_oldE_Bin2, linestyle="dotted", colour='blue',label = r"old $n_{E}$ in Bin 2", AX=axs0)
    plot_compact_sum(my_nd_E_Bin2, linestyle = "dotted", colour='red',label = r"new $n_{E}$ in Bin 2", AX=axs0)

    plot_compact_sum(my_nd_oldE_Bin3, linestyle="dashdot", colour='blue',label = r"old $n_{E}$ in Bin 3", AX=axs0)
    plot_compact_sum(my_nd_E_Bin3, linestyle = "dashdot", colour='red',label = r"new $n_{E}$ in Bin 3", AX=axs0)

    
    plot_nd_3bins(my_nd_Barro,'o',AX=axs0)

    
    plot_nd_Barro(axs0)
    plt.setp(axs0.get_xticklabels(), visible=False)
    #axs0.grid(True)
    axs0.legend(loc=4)


    
    #plot Panel (2)
    axs1 = plt.subplot(gs[1],sharex=axs0) 
    

    plot_compact_sum(my_nd_oldE_Bin1, linestyle="dashed", colour='blue',label = r"old $n_{E}$ in Bin 1", AX=axs1)
    plot_compact_sum(my_nd_E_Bin1, linestyle = "dashed", colour='red',label = r"new $n_{E}$ in Bin 1", AX=axs1)

    plot_compact_sum(my_nd_oldE_Bin2, linestyle="dotted", colour='blue',label = r"old $n_{E}$ in Bin 2", AX=axs1)
    plot_compact_sum(my_nd_E_Bin2, linestyle = "dotted", colour='red',label = r"new $n_{E}$ in Bin 2", AX=axs1)

    plot_compact_sum(my_nd_oldE_Bin3, linestyle="dashdot", colour='blue',label = r"old $n_{E}$ in Bin 3", AX=axs1)
    plot_compact_sum(my_nd_E_Bin3, linestyle = "dashdot", colour='red',label = r"new $n_{E}$ in Bin 3", AX=axs1)
    
    plot_nd_3bins(my_nd_vdWel,'s',AX=axs1)
    
    plot_nd_vdWel(axs1)
    plt.setp(axs1.get_xticklabels(), visible=False)
    #axs1.grid(True)
    axs1.legend(loc=4)


    
    #plot Panel (3)
    axs2 = plt.subplot(gs[2],sharex=axs0)  
    

    plot_compact_sum(my_nd_oldE_Bin1, linestyle="dashed", colour='blue',label = r"old $n_{E}$ in Bin 1", AX=axs2)
    plot_compact_sum(my_nd_E_Bin1, linestyle = "dashed", colour='red',label = r"new $n_{E}$ in Bin 1", AX=axs2)

    plot_compact_sum(my_nd_oldE_Bin2, linestyle="dotted", colour='blue',label = r"old $n_{E}$ in Bin 2", AX=axs2)
    plot_compact_sum(my_nd_E_Bin2, linestyle = "dotted", colour='red',label = r"new $n_{E}$ in Bin 2", AX=axs2)

    plot_compact_sum(my_nd_oldE_Bin3, linestyle="dashdot", colour='blue',label = r"old $n_{E}$ in Bin 3", AX=axs2)
    plot_compact_sum(my_nd_E_Bin3, linestyle = "dashdot", colour='red',label = r"new $n_{E}$ in Bin 3", AX=axs2)

    plot_nd_3bins(my_nd_vDokkum,'^',AX=axs2)
    
    plot_nd_vDokkum(axs2)

    axs2.set_xlabel(r"$ z$",fontsize=18)
    #axs2.grid(True)
    axs2.legend(loc=4)
    
    #plot Panel (4)
    axs3= plt.subplot(gs[3],sharex=axs0)  

    
    plot_compact_sum(my_nd_oldE_Bin1, linestyle="dashed", colour='blue',label = r"old $n_{E}$ in Bin 1", AX=axs3)
    plot_compact_sum(my_nd_E_Bin1, linestyle = "dashed", colour='red',label = r"new $n_{E}$ in Bin 1", AX=axs3)

    plot_compact_sum(my_nd_oldE_Bin2, linestyle="dotted", colour='blue',label = r"old $n_{E}$ in Bin 2$", AX=axs3)
    plot_compact_sum(my_nd_E_Bin2, linestyle = "dotted", colour='red',label = r"new $n_{E}$ in Bin 2", AX=axs3)

    plot_compact_sum(my_nd_oldE_Bin3, linestyle="dashdot", colour='blue',label = r"old $n_{E}$ in Bin 3", AX=axs3)
    plot_compact_sum(my_nd_E_Bin3, linestyle = "dashdot", colour='red',label = r"new $n_{E}$ in Bin 3", AX=axs3)
    plot_nd_3bins(my_nd_Dam,'^',AX=axs3)
    
    plot_nd_Dam(axs3)

    axs3.set_xlabel(r"$ z$",fontsize=18)
    #axs2.grid(True)
    axs3.legend(loc=4)
    
    plt.show()
    

def plot_nd_3plot():
    
    fig = plt.figure()
    gs = gridspec.GridSpec(ncols=1, nrows=4,
                               hspace=0.1, wspace=0.0) 

    #plot Panel (1)
    axs0 = plt.subplot(gs[0])  
    
    plot_compact_sum(my_nd_Barro, colour='black',label = r"$n_{sum}$", AX=axs0)
    plot_compact_sum(my_nd_E, colour='red',label = r"$n_{E}$", AX=axs0)    
    plot_nd_3bins(my_nd_Barro,'o',AX=axs0)

    
    plot_nd_Barro(axs0)
    plt.setp(axs0.get_xticklabels(), visible=False)
    #axs0.grid(True)
    axs0.legend(loc=4)


    
    #plot Panel (2)
    axs1 = plt.subplot(gs[1],sharex=axs0) 
    

    plot_compact_sum(my_nd_vdWel, colour='black',label = r"$n_{sum}$", AX=axs1)
    plot_compact_sum(my_nd_E, colour='red',label = r"$n_{E}$", AX=axs1)    
    plot_nd_3bins(my_nd_vdWel,'s',AX=axs1)
    
    plot_nd_vdWel(axs1)
    plt.setp(axs1.get_xticklabels(), visible=False)
    #axs1.grid(True)
    axs1.legend(loc=4)


    
    #plot Panel (3)
    axs2 = plt.subplot(gs[2],sharex=axs0)  

    

    plot_compact_sum(my_nd_vDokkum, colour='black',label = r"$n_{sum}$", AX=axs2)
    plot_compact_sum(my_nd_E, colour='red',label = r"$n_{E}$", AX=axs2)    
    plot_nd_3bins(my_nd_vDokkum,'^',AX=axs2)
    
    plot_nd_vDokkum(axs2)

    axs2.set_xlabel(r"$ z$",fontsize=18)
    #axs2.grid(True)
    axs2.legend(loc=4)
    
    #plot Panel (4)
    axs3= plt.subplot(gs[3],sharex=axs0)  


    plot_compact_sum(my_nd_Dam, colour='black',label = r"$n_{sum}$", AX=axs3)
    plot_compact_sum(my_nd_E, colour='red',label = r"$n_{E}$", AX=axs3)    
    plot_nd_3bins(my_nd_Dam,'^',AX=axs3)
    
    plot_nd_Dam(axs3)

    axs3.set_xlabel(r"$ z$",fontsize=18)
    #axs2.grid(True)
    axs3.legend(loc=4)
    
    plt.show()
    
plot_nd_3plot()
plot_nd_E_3plot()