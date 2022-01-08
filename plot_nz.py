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
mpl.rcParams["legend.numpoints"] = 1.0
mpl.rcParams["legend.scatterpoints"] = 1.0
#####

D = np.array([45,75,110])

voll = ((D**3)/3)*((214-139)/180)*np.pi*(np.cos((np.pi/2)-(55*np.pi/180)-np.cos(np.pi/2)))

print(voll)
#V1,V2,V3=voll[2],voll[1],voll[0]

V1_V = voll[2]-voll[1]
V2_V = voll[1]-voll[0]
V3_V = voll[0]

V1,V2,V3 = V1_V,V2_V,V3_V

V = np.array([V1_V, V2_V, V3_V])

print('-----volume-----')
print(voll[2], voll[1], voll[0],'sum')
print(V1_V,V2_V,V3_V,'sum',V1_V+V2_V+V3_V)

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

# Damjanov 2015 Barro cut
Dam_data = SRead.read_table("/home/dexter/result/stat/completeness/nd_Damjanov15.dat")
Dam_z = Dam_data[:,0]
Dam_nd = 10**Dam_data[:,1]
Dam_nd_uerr = 10**Dam_data[:,2]

Dam_data1 = SRead.read_table("/home/dexter/result/stat/completeness/nd_Damjanov15_COSMOS.dat")
Dam_z1 = Dam_data1[:,0]
Dam_nd1 = 10**Dam_data1[:,1]
Dam_nd_uerr1 = 10**Dam_data1[:,2]

Dam_data2 = SRead.read_table("/home/dexter/result/stat/completeness/nd_Damjanov14_BOSS.dat")
Dam_z2 = Dam_data2[:,0]
Dam_nd2 = 10**Dam_data2[:,1]
Dam_nd_uerr2 = 10**Dam_data2[:,2]

Dam_data3 = SRead.read_table("/home/dexter/result/stat/completeness/nd_Damjanov15_vdwel.dat")
Dam_z3 = Dam_data3[:,0]
Dam_nd3 = 10**Dam_data3[:,1]
Dam_nd_uerr3 = 10**Dam_data3[:,2]

print("max",max(Barro_nd),max(vDokkum_nd),max(vdWel_nd))

Poggianti_z = 0.06009445742113899
Poggianti_nd = 10**-3.8849206349206353
Poggianti_nd_u = 10**-3.5238095238095237

my_z = np.array([0.0127,0.0085,0.005])

my_z = np.array([0,0,0])
###
#my_nd_Barro = np.array([6.75e-4,1.19e-4,2.94e-5])
#my_nd_Dam = np.array([7.68e-4,1.39e-4,2.94e-5])
#my_nd_vDokkum = np.array([2.14e-4, 1.99e-5, 6.31e-6])
#my_nd_vdWel = np.array([2.76e-4,2.65e-5,6.31e-6])
#
#my_nd_Graham = np.array([1.54e-4,1.33e-5,0])
#
#####
#Bin H number density
#my_nd_Barro = np.array([4.41e-5,1.19e-4,4.3e-4])
#my_nd_Dam = np.array([4.83e-5,1.33e-4,4.3e-4])
#my_nd_vDokkum = np.array([8.41e-6, 1.33e-5, 9.21e-5])
#my_nd_vdWel = np.array([1.89e-5,2.65e-5,9.21e-5])
#
#my_nd_Graham = np.array([8.41e-6,6.63e-6,0])

#####
#Bin V number density RC15
my_nd_Barro = np.array([5.23e-5,1.07e-4,7.06e-4])
my_nd_Dam = np.array([5.23e-5,1.35e-4,7.98e-4])
my_nd_vDokkum = np.array([1.54e-5, 2.54e-5, 3.07e-5])
my_nd_vdWel = np.array([1.54e-5,2.54e-5,6.14e-5])
my_nd_Graham = np.array([1.54e-5,1.69e-5,0])

# The lower limit for the overall num dens
my_nd_sum_Barro_low = sum(my_nd_Barro*V)/voll[2]
my_nd_sum_Dam_low = sum(my_nd_Dam*V)/voll[2]
my_nd_sum_vDokkum_low =  sum(my_nd_vDokkum*V)/voll[2]
my_nd_sum_vdWel_low =  sum(my_nd_vdWel*V)/voll[2]
my_nd_sum_Graham_low =  sum(my_nd_Graham*V)/voll[2]

my_nd_sum_E_low = (7+ 3+ 3)/voll[2]

#Bin V number density T11
my_nd_Barro_T11 = np.array([3.69e-5,7.61e-5,4.61e-4])
my_nd_Dam_T11 = np.array([4.92e-5,8.46e-5,5.83e-4])
my_nd_vDokkum_T11 = np.array([0,0,0])
my_nd_vdWel_T11 = np.array([0,0,0])
my_nd_Graham_T11 = np.array([3.08e-6,0,0])

#Bin V numbder density Z09
my_nd_Barro_Z09 = np.array([4.62e-5,1.02e-4,6.14e-4])
my_nd_Dam_Z09 = np.array([5.23e-5,1.10e-4,7.37e-4])
my_nd_vDokkum_Z09 = np.array([6.15e-6,8.46e-6,0])
my_nd_vdWel_Z09 = np.array([3.08e-6,2.54e-5,0])
my_nd_Graham_Z09 = np.array([9.23e-6,8.46e-6,0])

#Bin V numer density IP13
my_nd_Barro_IP13 = np.array([6.15e-5,1.44e-4,1.04e-3])
my_nd_Dam_IP13 = np.array([6.46e-5,1.69e-4,1.17e-3])
my_nd_vDokkum_IP13 = np.array([4.62e-5,8.46e-5,2.76e-4])
my_nd_vdWel_IP13 = np.array([3.39e-5,6.77e-5,1.84e-4])
my_nd_Graham_IP13 = np.array([1.85e-5,4.23e-5,9.21e-5])
#my_nd_Barro = my_nd_Barro_T11
#my_nd_Dam = my_nd_Dam_T11
#my_nd_vDokkum = my_nd_vDokkum_T11
#my_nd_vdWel = my_nd_vdWel_T11
#my_nd_Graham = my_nd_Graham_T11

#my_nd_E = np.array([6/V1, 3/V2, 2/V3])#8.97e-5]
my_nd_E = np.array([7/V1, 3/V2, 3/V3])#8.97e-5] BinV

my_nd_E_Bin1 = np.array([my_nd_E[0],0])
my_nd_E_Bin2 = np.array([my_nd_E[1],0])
my_nd_E_Bin3 = np.array([my_nd_E[2],0])


#my_nd_oldE = np.array([15/V1, 11/V2, 2/V3])
my_nd_oldE = np.array([13/V1, 9/V2, 6/V3]) #BinV

my_nd_oldE_Bin1 = np.array([my_nd_oldE[0],0])
my_nd_oldE_Bin2 = np.array([my_nd_oldE[1],0])
my_nd_oldE_Bin3 = np.array([my_nd_oldE[2],0])

peak_RN_Barro = 2.413e-4
peak_RN_vdWel = 1.692e-4
peak_RN_vDokkum = 1.458e-4 
peak_RN_Dam = 5.08e-5


ms0 = 12

print(np.sum(my_nd_Barro), np.sum(my_nd_vdWel), np.sum(my_nd_vDokkum))


print("nd_E+ES",my_nd_E,sum(my_nd_E))

print("Csph", my_nd_Barro,sum(my_nd_Barro))

print("cSph to E ratio:",
      np.sum(my_nd_Barro)/np.sum(my_nd_E), 
      np.sum(my_nd_vdWel)/np.sum(my_nd_E) , 
      np.sum(my_nd_vDokkum)/np.sum(my_nd_E))

print("E-peak RN /peak RN:",
      (np.sum(my_nd_E)/max(Barro_nd)), 
      (np.sum(my_nd_E)/max(vdWel_nd)), 
      (np.sum(my_nd_E)/max(vDokkum_nd)))

print("csph-peak /peak RN:",
      (np.sum(my_nd_Barro)/max(Barro_nd)), 
      (np.sum(my_nd_vdWel)/max(vdWel_nd)), 
      (np.sum(my_nd_vDokkum)/max(vDokkum_nd)))

xlim = [-0.35,3.0]
ylim = [0.8e-6,1.2e-3]

xlim_n = [0.05,0.15]

def plot_compact_sum(nd0, linestyle = "solid", colour="black", label = "", 
                     xlim_n = xlim_n, 
                     AX=plt, nsum=True):
    """
    plot the sum of compact sample and the Ellipitcals

    Returns
    -------
    None.

    """ 
    #nd = sum(nd0)
    if nsum == True:
        nd = sum(nd0*V)/voll[2]

    elif nsum == False:
        nd = nd0

    AX.axhline(nd,xlim_n[0],xlim_n[1],linestyle=linestyle, linewidth = 3,
               color=colour , label= label)
    print("sum(nd0*V)",sum(nd0*V),"nd:",nd)
    

def plot_nd_3bins(nd,marker, my_z=my_z,AX=plt):
    n_err_bin1 = np.sqrt(nd[0]*V1)/V1
    n_err_bin2 = np.sqrt(nd[1]*V2)/V2
    n_err_bin3 = np.sqrt(nd[2]*V3)/V3
    
    xedge1, yedge1 = [-5,5,5,-5],[nd[0]-n_err_bin1, nd[0]-n_err_bin1,
                                  nd[0]+n_err_bin1, nd[0]+n_err_bin1]
    xedge2, yedge2 = [-5,5,5,-5],[nd[1]-n_err_bin2, nd[1]-n_err_bin2,
                                  nd[1]+n_err_bin2, nd[1]+n_err_bin2]
    xedge3, yedge3 = [-5,5,5,-5],[nd[2]-n_err_bin3, nd[2]-n_err_bin3,
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

    #AX.scatter(my_z[0],nd[0], s= 50, c ='#a5200b',marker=marker)
    #AX.scatter(my_z[1],nd[1], s= 50, c ='#0b5786',marker=marker)
    #AX.scatter(my_z[2],nd[2], s= 50, c ='#2a3236',marker=marker)

def plot_nd_Dam(AX):
    AX.plot(Dam_z,Dam_nd,'--d',lw =5,  ms= ms0, color='r',
            label=r"$\rm Damjanov~et~al.~(2015)~(Barro~criteria)$")
    AX.errorbar(Dam_z,Dam_nd,yerr=Dam_nd_uerr-Dam_nd,ls='none',
                     linewidth=3, ecolor='r',mew=1,capsize=3) 


    xedge, yedge = list(Dam_z), list(Dam_nd)
    xedge.append(Dam_z[-1])
    yedge.append(0)
    xedge.append(Dam_z[0])
    yedge.append(0)

    #AX.fill(xedge, yedge, alpha=0.1, color='g')
    #AX.set_xlabel(r"$ z$",fontsize=18)
    
    AX.set_ylabel(r"$n \rm (Mpc^{-3}$)",fontsize= 20)

    AX.set_ylim(ylim[0],ylim[1])
    AX.set_xlim(xlim[0],xlim[1])

    AX.set_yscale( 'log' )
    #AX.legend() 

def plot_nd_Dam1(AX):
    AX.plot(Dam_z1,Dam_nd1,'--d',lw =5,  ms= ms0, color='r',
            label=r"$\rm Damjanov~et~al.~(2015)$")

    AX.errorbar(Dam_z1,Dam_nd1,yerr=Dam_nd_uerr1-Dam_nd1,ls='none',
                     linewidth=3, ecolor='r',mew=1,capsize=3) 


    xedge, yedge = list(Dam_z1), list(Dam_nd1)
    xedge.append(Dam_z1[-1])
    yedge.append(0)
    xedge.append(Dam_z1[0])
    yedge.append(0)

    #AX.fill(xedge, yedge, alpha=0.1, color='g')
    #AX.set_xlabel(r"$ z$",fontsize=18)
    
    AX.set_ylabel(r"$n \rm (Mpc^{-3}$)",fontsize= 20)

    AX.set_ylim(ylim[0],ylim[1])
    AX.set_xlim(xlim[0],xlim[1])

    AX.set_yscale( 'log' )
    #AX.legend() 
    
def plot_nd_Dam2(AX):
    AX.plot(Dam_z2,Dam_nd2,'--d',lw =5,  ms= ms0, color='r',
            label=r"$\rm Damjanov~et~al.~(2015)$")
    AX.errorbar(Dam_z2,Dam_nd2,yerr=Dam_nd_uerr2-Dam_nd2,ls='none',
                     linewidth=3, ecolor='r',mew=1,capsize=3) 


    xedge, yedge = list(Dam_z2), list(Dam_nd2)
    xedge.append(Dam_z2[-1])
    yedge.append(0)
    xedge.append(Dam_z2[0])
    yedge.append(0)

    #AX.fill(xedge, yedge, alpha=0.1, color='g')
    #AX.set_xlabel(r"$ z$",fontsize=18)
    
    AX.set_ylabel(r"$n \rm (Mpc^{-3}$)",fontsize= 20)

    AX.set_ylim(ylim[0],ylim[1])
    AX.set_xlim(xlim[0],xlim[1])

    AX.set_yscale( 'log' )
    #AX.legend() 

def plot_nd_Dam3(AX):
    AX.plot(Dam_z3,Dam_nd3,'--d',lw =5,  ms= ms0, color='r',
            label=r"$\rm Damjanov~et~al.~(2015)~(van~der~Wel~criteria)$")
    AX.errorbar(Dam_z3,Dam_nd3,yerr=Dam_nd_uerr3-Dam_nd3,ls='none',
                     linewidth=3, ecolor='r',mew=1,capsize=3) 


    xedge, yedge = list(Dam_z3), list(Dam_nd3)
    xedge.append(Dam_z3[-1])
    yedge.append(0)
    xedge.append(Dam_z3[0])
    yedge.append(0)

    #AX.fill(xedge, yedge, alpha=0.1, color='g')
    #AX.set_xlabel(r"$ z$",fontsize=18)
    
    AX.set_ylabel(r"$n \rm (Mpc^{-3}$)",fontsize= 20)


    AX.set_ylim(ylim[0],ylim[1])
    AX.set_xlim(xlim[0],xlim[1])

    AX.set_yscale( 'log' )
    #AX.legend() 

def plot_nd_Barro(AX):

    AX.plot(Barro_z,Barro_nd,'--d', lw =5,  ms= ms0, color='g',
            label=r"$\rm Barro~et~al.~(2013)$")
    AX.errorbar(Barro_z,Barro_nd,yerr=Barro_nd_uerr-Barro_nd,ls='none',
                     linewidth=3, ecolor='g',mew=1,capsize=3) 


    xedge, yedge = list(Barro_z), list(Barro_nd)
    xedge.append(Barro_z[-1])
    yedge.append(0)
    xedge.append(Barro_z[0])
    yedge.append(0)

    #AX.fill(xedge, yedge, alpha=0.1, color='g')
    #AX.set_xlabel(r"$ z$",fontsize=18)
    
    AX.set_ylabel(r"$n \rm (Mpc^{-3}$)",fontsize= 20)


    AX.set_ylim(ylim[0],ylim[1])
    AX.set_xlim(xlim[0],xlim[1])

    AX.set_yscale( 'log' )
    #AX.legend()
            

def plot_nd_vDokkum(AX):

    AX.plot(vDokkum_z,vDokkum_nd,'--d', lw=5, ms= ms0,color='#e58b1a', 
            label=r"$\rm van~Dokkum~et~al.~(2015)$")
    AX.errorbar(vDokkum_z,vDokkum_nd,yerr= vDokkum_nd_uerr - vDokkum_nd,ls='none',
                     linewidth=3, ecolor='#e58b1a',mew=1,capsize=3) 


    xedge, yedge = list(vDokkum_z), list(vDokkum_nd)
    xedge.append(vDokkum_z[-1])
    yedge.append(0)
    xedge.append(vDokkum_z[0])
    yedge.append(0)

    #AX.fill(xedge, yedge, alpha=0.1, color='y')
    #AX.set_xlabel(r"$ z$",fontsize=18)
    AX.set_ylabel(r"$ n (Mpc^{-3})$",fontsize= 20)
    AX.set_yscale( 'log' )
    #AX.legend()
        
    AX.set_ylim(ylim[0],ylim[1])
    AX.set_xlim(xlim[0],xlim[1])
    
def plot_nd_vdWel(AX):


    AX.plot(vdWel_z,vdWel_nd,'--d', lw=5, ms= ms0, color='b', 
            label=r"$\rm van~der~Wel~et~al.~(2014)$")
    AX.errorbar(vdWel_z,vdWel_nd,yerr= vdWel_nd_uerr - vdWel_nd,ls='none',
                     linewidth=3, ecolor='b',mew=1,capsize=3) 

    xedge, yedge = list(vdWel_z), list(vdWel_nd)
    xedge.append(vdWel_z[-1])
    yedge.append(0)
    xedge.append(vdWel_z[0])
    yedge.append(0)

    #AX.fill(xedge, yedge, alpha=0.1, color='b')
    #AX.set_xlabel(r"$ z$",fontsize=18)
    AX.set_ylabel(r"$ n (Mpc^{-3})$",fontsize= 20)
    AX.set_yscale( 'log' )
    #AX.legend()
    
    AX.set_ylim(ylim[0],ylim[1])
    AX.set_xlim(xlim[0],xlim[1])    
    
#############################
#
import matplotlib.gridspec as gridspec


def plot_nd_E_3plot():
        
    fig = plt.figure(figsize=(6.4, 4.8))
    gs = gridspec.GridSpec(ncols=1, nrows=4,
                               hspace=0.0, wspace=0.0) 

    #plot Panel (1)
    axs0 = plt.subplot(gs[0])  
    
    plot_compact_sum(my_nd_oldE_Bin1, linestyle="dashed", colour='blue',label = r"old $n_{E+ES}$ in Bin 1", AX=axs0)
    plot_compact_sum(my_nd_E_Bin1, linestyle = "dashed", colour='red',label = r"new $n_{E+ES}$ in Bin 1", AX=axs0)

    plot_compact_sum(my_nd_oldE_Bin2, linestyle="dotted", colour='blue',label = r"old $n_{E+ES}$ in Bin 2", AX=axs0)
    plot_compact_sum(my_nd_E_Bin2, linestyle = "dotted", colour='red',label = r"new $n_{E+ES}$ in Bin 2", AX=axs0)

    plot_compact_sum(my_nd_oldE_Bin3, linestyle="dashdot", colour='blue',label = r"old $n_{E+ES}$ in Bin 3", AX=axs0)
    plot_compact_sum(my_nd_E_Bin3, linestyle = "dashdot", colour='red',label = r"new $n_{E+ES}$ in Bin 3", AX=axs0)

    
    plot_nd_3bins(my_nd_Barro,'o',AX=axs0)
    
    
    plot_nd_Barro(axs0)
    plt.setp(axs0.get_xticklabels(), visible=False)
    #axs0.grid(True)
    axs0.legend(loc=4,fontsize=10)


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
    axs1.legend(loc=4,fontsize=10)


    
    #plot Panel (3)
    axs2 = plt.subplot(gs[2],sharex=axs0)  

    plot_compact_sum(my_nd_oldE_Bin1, linestyle="dashed", colour='blue',label = r"old $n_{E+ES}$ in Bin 1", AX=axs2)
    plot_compact_sum(my_nd_E_Bin1, linestyle = "dashed", colour='red',label = r"new $n_{E+ES}$ in Bin 1", AX=axs2)

    plot_compact_sum(my_nd_oldE_Bin2, linestyle="dotted", colour='blue',label = r"old $n_{E+ES}$ in Bin 2", AX=axs2)
    plot_compact_sum(my_nd_E_Bin2, linestyle = "dotted", colour='red',label = r"new $n_{E+ES}$ in Bin 2", AX=axs2)

    plot_compact_sum(my_nd_oldE_Bin3, linestyle="dashdot", colour='blue',label = r"old $n_{E+ES}$ in Bin 3", AX=axs2)
    plot_compact_sum(my_nd_E_Bin3, linestyle = "dashdot", colour='red',label = r"new $n_{E+ES}$ in Bin 3", AX=axs2)

    plot_nd_3bins(my_nd_vDokkum,'o',AX=axs2)
    
    plot_nd_vDokkum(axs2)

    #axs2.grid(True)
    axs2.legend(loc=4,fontsize=10)
    
    #plot Panel (4)
    axs3= plt.subplot(gs[3],sharex=axs0)  

    
    plot_compact_sum(my_nd_oldE_Bin1, linestyle="dashed", colour='blue',label = r"old $n_{E+ES}$ in Bin 1", AX=axs3)
    plot_compact_sum(my_nd_E_Bin1, linestyle = "dashed", colour='red',label = r"new $n_{E+ES}$ in Bin 1", AX=axs3)

    plot_compact_sum(my_nd_oldE_Bin2, linestyle="dotted", colour='blue',label = r"old $n_{E+ES}$ in Bin 2$", AX=axs3)
    plot_compact_sum(my_nd_E_Bin2, linestyle = "dotted", colour='red',label = r"new $n_{E+ES}$ in Bin 2", AX=axs3)

    plot_compact_sum(my_nd_oldE_Bin3, linestyle="dashdot", colour='blue',label = r"old $n_{E+ES}$ in Bin 3", AX=axs3)
    plot_compact_sum(my_nd_E_Bin3, linestyle = "dashdot", colour='red',label = r"new $n_{E+ES}$ in Bin 3", AX=axs3)
    plot_nd_3bins(my_nd_Dam,'o',AX=axs3)
    
    plot_nd_Dam(axs3)

    axs3.set_xlabel(r"$ z$",fontsize=18)
    #axs2.grid(True)
    axs3.legend(loc=4,fontsize=10)
    plt.tight_layout()

    plt.show()
    

def plot_nd_3plot():
    
    fig = plt.figure(figsize=(6.4, 14))
    gs = gridspec.GridSpec(ncols=1, nrows=4,
                               hspace=0.0, wspace=0.0) 

    #plot Panel (1)
    axs0 = plt.subplot(gs[0])  
    

    axs0.plot(Poggianti_z,Poggianti_nd,'d', lw=5, ms= ms0, color='purple', 
              label=r"$\rm Poggianti~et~al.~(2013)$")
    #axs0.errorbar(Poggianti_z,Poggianti_nd,yerr= Poggianti_nd_u - Poggianti_nd
    #              ,ls='none',linewidth=3, ecolor='purple',mew=1,capsize=3) 

    plot_nd_3bins(my_nd_Barro,'o',AX=axs0)
    print('---n_c,Sph (Barro)-------')
    plot_compact_sum(my_nd_Barro, colour='black',label = r"", AX=axs0)    
    print('---n_E+ES--------')
    plot_compact_sum(my_nd_E, colour='red',label = r"", AX=axs0)    
    
    print('-----------------')


    plot_nd_Barro(axs0)
    plot_nd_Dam(axs0)

    plt.setp(axs0.get_xticklabels(), visible=False)
    #axs0.grid(True)
    axs0.text(1.6,4e-4, r'$\rm Barro~criteria$', weight = "bold", color= "g", fontsize=24)
    axs0.legend(loc=4,fontsize=11)

    twin0=axs0.twinx()
    
    #plot Panel (2)
    axs1 = plt.subplot(gs[1],sharex=axs0) 
       
    plot_nd_3bins(my_nd_vdWel,'o',AX=axs1)
    
    print('---n_c,Sph (vdWel)-------')
    plot_compact_sum(my_nd_vdWel, colour='black',label = r"", AX=axs1)
    print('---n_E+ES--------')
    plot_compact_sum(my_nd_E, colour='red',label = r"", AX=axs1) 
    print('-----------------')

    plot_nd_vdWel(axs1)
    plot_nd_Dam3(axs1)
    plt.setp(axs1.get_xticklabels(), visible=False)
    #axs1.grid(True)
    axs1.text(1.0,4e-4, r'$\rm van~der~Wel~criteria$', weight = "bold", color= "b", fontsize=24)
    axs1.legend(loc=4,fontsize=11)

    #plot Panel (3)
    axs2 = plt.subplot(gs[2],sharex=axs0)  

  
    plot_nd_3bins(my_nd_vDokkum,'o',AX=axs2)
    print('---n_c,Sph (vDokkum)-------')
    plot_compact_sum(my_nd_vDokkum, colour='black',label = r"", AX=axs2)
    print('---n_E+ES--------')
    plot_compact_sum(my_nd_E, colour='red',label = r"", AX=axs2)  
    print('-----------------')

    plot_nd_vDokkum(axs2)
    plt.setp(axs2.get_xticklabels(), visible=False)
    
    #axs2.grid(True)
    axs2.text(0.9,4e-4, r'$\rm van~Dokkum~criteria$', weight = "bold", color= "#e58b1a", fontsize=24)
    axs2.legend(loc=4,fontsize=11)
    
    #plot Panel (4)
    axs3= plt.subplot(gs[3],sharex=axs0)  
    
    plot_nd_3bins(my_nd_Dam,'o',AX=axs3)

    print('---n_c,Sph (Dam)-------')
    plot_compact_sum(my_nd_Dam, colour='black',label = r"", AX=axs3)

    print('---n_E+ES (Dam)-------')

    plot_compact_sum(my_nd_E, colour='red',label = r"", AX=axs3)  
    print('-----------------')
   
    plot_nd_Dam1(axs3)
    #plot_nd_Dam2(axs3)

    axs3.set_xlabel(r"$ z$",fontsize=20)
    #axs2.grid(True)
    axs3.text(1.2,4e-4, r'$\rm Damjanov~criteria$', weight = "bold", color= "r", fontsize=24)
    axs3.legend(loc=4,fontsize=11)
    
    twin0.scatter([],[],label=r"$\rm c,Sph~in~Bin~1 $", color ='#a5200b', 
               marker ="o")
    twin0.scatter([],[],label=r"$\rm c,Sph~in~Bin~2 $", color ='#0b5786', 
               marker ="o")
    twin0.scatter([],[],label=r"$\rm c,Sph~in~Bin~3 $", color ='#2a3236', 
               marker ="o")


    plot_compact_sum([np.nan], colour='black',label = r"$n_\mathrm{c,Sph}$", AX=twin0)
    plot_compact_sum([np.nan], colour='red',label = r"$n_\mathrm{E+ES}$", AX=twin0) 

    twin0.legend(loc='upper center',fontsize=13, bbox_to_anchor=(0.5, 1.3),
          fancybox=True, shadow=False, ncol=3)
    twin0.set_yscale('log')
    plt.setp(twin0.get_yticklabels(), visible=False)
    plt.tight_layout()

    plt.show()
    
    
def plot_nd_all_mass():
    
    ylim_allmass = [0.1e-6,2e-3,2e-3,0.1e-6]
    xlim_allmass = [-0.2,-0.2,0.2,0.2]
    alpha_allmass = 0.25
    
    fig = plt.figure(figsize=(6.4, 12))
    gs = gridspec.GridSpec(ncols=4, nrows=4,
                               hspace=0.0, wspace=0.1) 

    # T11 Barro cut
    axs0 = plt.subplot(gs[0])  
    
    plot_compact_sum(my_nd_Barro_T11, colour='black', xlim_n = xlim,
                     label = r"$n_\mathrm{c,Sph}$", 
                     AX=axs0,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"$n_\mathrm{E+ES}$", 
                     AX=axs0,nsum = True) 
    axs0.axhline(peak_RN_Barro,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="g" )
    axs0.axhline(0,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="k",label=r"$\rm max(n_\mathrm{RN})$" )
    
    plot_nd_3bins(my_nd_Barro_T11,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs0)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs0)

    
    #plt.setp(axs0.get_xticklabels(), visible=False)
    axs0.set_ylabel(r"$ n (\rm M p c^{-3})$",fontsize=16)

    axs0.set_ylim(0.1e-6,2e-3)
    axs0.set_xlim(-0.2,0.2)
    axs0.set_yscale('log')    
    axs0.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='g')
    plt.xticks([-0.1, 0.1],[],fontsize=18)

    # Z09 Barro cut
    axs1 = plt.subplot(gs[1])  

    plot_compact_sum(my_nd_Barro_Z09, colour='black', xlim_n = xlim,
                     label = r"", AX=axs1, nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs1,nsum = True) 
    axs1.axhline(peak_RN_Barro,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="g" ) 
    
    plot_nd_3bins(my_nd_Barro_Z09,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs1)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs1)
     
    axs1.set_ylim(0.1e-6,2e-3)
    axs1.set_xlim(-0.2,0.2)
    axs1.set_yscale('log')
    axs1.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='g')
    plt.xticks([-0.1, 0.1],[],fontsize=18)

    
    #plt.setp(axs1.get_xticklabels(), visible=False)
    plt.setp(axs1.get_yticklabels(), visible=False)

    # RC15 Barro cut
    axs2 = plt.subplot(gs[2])  
    
    plot_compact_sum(my_nd_Barro, colour='black', xlim_n = xlim,
                     label = r"", AX=axs2,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs2,nsum = True) 
    axs2.axhline(peak_RN_Barro,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="g" )
    
    plot_nd_3bins(my_nd_Barro,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs2)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs2)
     
    axs2.set_ylim(0.1e-6,2e-3)
    axs2.set_xlim(-0.2,0.2)
    axs2.set_yscale('log')
    axs2.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='g')
    plt.xticks([-0.1, 0.1],[],fontsize=18)

    #plt.setp(axs2.get_xticklabels(), visible=False)
    plt.setp(axs2.get_yticklabels(), visible=False)
    
    # IP13 Barro cut
    axs3 = plt.subplot(gs[3])  
    
    plot_compact_sum(my_nd_Barro_IP13, colour='black', xlim_n = xlim,
                     label = r"", AX=axs3,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs3,nsum = True) 
    axs3.axhline(peak_RN_Barro,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="g" )
    
    plot_nd_3bins(my_nd_Barro_IP13,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs3)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs3)
     
    axs3.set_ylim(0.1e-6,2e-3)
    axs3.set_xlim(-0.2,0.2)
    axs3.set_yscale('log')
    
    twin3=axs3.twinx()
    twin3.set_ylabel(r'$\rm Barro~criteria$',fontsize=16)
    
    axs3.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='g')
    plt.xticks([-0.1, 0.1],[],fontsize=18)

    #plt.setp(axs3.get_xticklabels(), visible=False)
    plt.setp(axs3.get_yticklabels(), visible=False)
    #plt.setp(twin3.get_xticklabels(), visible=False)
    plt.setp(twin3.get_yticklabels(), visible=False)
    # T11 vdWel cut
    axs4 = plt.subplot(gs[4])  
    plot_compact_sum(my_nd_vdWel_T11, colour='black', xlim_n = xlim,
                     label = r"", AX=axs4,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs4,nsum = True) 
    axs4.axhline(peak_RN_vdWel,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="b" )
    
    plot_nd_3bins(my_nd_vdWel_T11,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs4)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs4)
     
    axs4.set_ylim(0.1e-6,2e-3)
    axs4.set_xlim(-0.2,0.2)
    axs4.set_yscale('log')

    #plt.setp(axs4.get_xticklabels(), visible=False)
    axs4.set_ylabel(r"$ n (\rm M p c^{-3})$",fontsize=16)
    
    axs4.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='b')
    plt.xticks([-0.1, 0.1],[],fontsize=18)

    # Z09 vdWel cut 
    axs5 = plt.subplot(gs[5])  
    plot_compact_sum(my_nd_vdWel_Z09, colour='black', xlim_n = xlim,
                     label = r"", AX=axs5,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs5,nsum = True) 
    axs5.axhline(peak_RN_vdWel,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="b" )
    
    plot_nd_3bins(my_nd_vdWel_Z09,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs5)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs5)
     
    axs5.set_ylim(0.1e-6,2e-3)
    axs5.set_xlim(-0.2,0.2)
    axs5.set_yscale('log')   
    
    plt.setp(axs5.get_xticklabels(), visible=False)
    plt.setp(axs5.get_yticklabels(), visible=False)
    
    axs5.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='b')
    plt.xticks([-0.1, 0.1],[],fontsize=18)

    # RC15 vdWel cut
    axs6 = plt.subplot(gs[6])  
    
    plot_compact_sum(my_nd_vdWel, colour='black', xlim_n = xlim,
                     label = r"", AX=axs6,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs6,nsum = True)
    axs6.axhline(peak_RN_vdWel,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="b" )
    
    plot_nd_3bins(my_nd_vdWel,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs6)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs6)
    
    axs6.set_ylim(0.1e-6,2e-3)
    axs6.set_xlim(-0.2,0.2)
    axs6.set_yscale('log')
    axs6.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='b')
    plt.xticks([-0.1, 0.1],[],fontsize=18)

    #plt.setp(axs6.get_xticklabels(), visible=False)
    plt.setp(axs6.get_yticklabels(), visible=False)
    
    # IP13 vdWel cut
    axs7 = plt.subplot(gs[7]) 
    plot_compact_sum(my_nd_vdWel_IP13, colour='black', xlim_n = xlim,
                     label = r"", AX=axs7,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs7,nsum = True) 
    axs7.axhline(peak_RN_vdWel,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="b" )
    
    plot_nd_3bins(my_nd_vdWel_IP13,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs7)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs7)
     
    axs7.set_ylim(0.1e-6,2e-3)
    axs7.set_xlim(-0.2,0.2)
    axs7.set_yscale('log')
        
    twin7 = axs7.twinx()
    #twin7.set_yscale('log')
    twin7.set_ylabel(r'$\rm van~der~Wel~criteria$',fontsize=16)
    #axs7.yaxis.tick_right()
    #axs7.yaxis.set_label_position("right")
    #axs7.set_ylabel(r'$\rm van~der~Wel~cut$',fontsize=16)
    
    #plt.setp(axs7.get_xticklabels(), visible=False)
    plt.setp(axs7.get_yticklabels(), visible=False)
    #plt.setp(twin7.get_xticklabels(), visible=False)
    plt.setp(twin7.get_yticklabels(), visible=False)   
    
    axs7.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='b')
    plt.xticks([-0.1, 0.1],[],fontsize=18)

    # T11 vDokkum cut
    axs8 = plt.subplot(gs[8]) 
    
    plot_compact_sum(my_nd_vDokkum_T11, colour='black', xlim_n = xlim,
                     label = r"", AX=axs8,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs8,nsum = True) 
    axs8.axhline(peak_RN_vDokkum,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="#b48f11" )
    
    plot_nd_3bins(my_nd_vDokkum_T11,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs8)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs8)
     
    axs8.set_ylim(0.1e-6,2e-3)
    axs8.set_xlim(-0.2,0.2)
    axs8.set_yscale('log')
    axs8.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='y')
    plt.xticks([-0.1, 0.1],[],fontsize=18)

    #plt.setp(axs8.get_xticklabels(), visible=False)
    axs8.set_ylabel(r"$ n (\rm M p c^{-3})$",fontsize=16)

    #Z09 vDokkum cut
    axs9 = plt.subplot(gs[9]) 
    
    plot_compact_sum(my_nd_vDokkum_Z09, colour='black', xlim_n = xlim,
                     label = r"", AX=axs9,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs9,nsum = True) 
    axs9.axhline(peak_RN_vDokkum,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="#b48f11" )
    
    plot_nd_3bins(my_nd_vDokkum_Z09,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs9)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs9)
     
    axs9.set_ylim(0.1e-6,2e-3)
    axs9.set_xlim(-0.2,0.2)
    axs9.set_yscale('log')
    
    #plt.setp(axs9.get_xticklabels(), visible=False)
    plt.setp(axs9.get_yticklabels(), visible=False)
    
    axs9.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='y')
    plt.xticks([-0.1, 0.1],[],fontsize=18)

    # RC15 vDokkum cut
    axs10 = plt.subplot(gs[10]) 
    plot_compact_sum(my_nd_vDokkum, colour='black', xlim_n = xlim,
                     label = r"", AX=axs10,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs10,nsum = True) 
    axs10.axhline(peak_RN_vdWel,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="#b48f11" )
    
    plot_nd_3bins(my_nd_vDokkum,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs10)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs10)
     
    axs10.set_ylim(0.1e-6,2e-3)
    axs10.set_xlim(-0.2,0.2)
    axs10.set_yscale('log')
    axs10.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='y')
    plt.xticks([-0.1, 0.1],[],fontsize=18)

    plt.setp(axs10.get_xticklabels(), visible=False)
    plt.setp(axs10.get_yticklabels(), visible=False)
    
    # IP13 vDokkum cut
    axs11 = plt.subplot(gs[11]) 
    plot_compact_sum(my_nd_vDokkum_IP13, colour='black', xlim_n = xlim,
                     label = r"", AX=axs11,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs11,nsum = True) 
    axs11.axhline(peak_RN_vdWel,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="#b48f11" )
    
    plot_nd_3bins(my_nd_vDokkum_IP13,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs11)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs11)
     
    axs11.set_ylim(0.1e-6,2e-3)
    axs11.set_xlim(-0.2,0.2)
    axs11.set_yscale('log')

    twin11 = axs11.twinx()
    twin11.set_ylabel(r'$\rm van~Dokkum~criteria$',fontsize=14)
    
    axs11.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='y')
    plt.xticks([-0.1, 0.1],[])

    #plt.setp(axs11.get_xticklabels(), visible=False)
    plt.setp(axs11.get_yticklabels(), visible=False)
    #plt.setp(twin11.get_xticklabels(), visible=False)
    plt.setp(twin11.get_yticklabels(), visible=False)    
    # T11 Dam cut
    axs12 = plt.subplot(gs[12]) 
    
    plot_compact_sum(my_nd_Dam_T11, colour='black', xlim_n = xlim,
                     label = r"", AX=axs12,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs12,nsum = True) 
    axs12.axhline(peak_RN_Dam,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="r" )
    
    plot_nd_3bins(my_nd_Dam_T11,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs12)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs12)    
    
    axs12.set_ylim(0.1e-6,2e-3)
    axs12.set_xlim(-0.2,0.2)
    axs12.set_yscale('log')
    axs12.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='r')
    plt.xticks([-0.1, 0.1],[r'$\rm c,Sph$', r'$\rm E+ES$'],fontsize=12)

    #plt.setp(axs12.get_xticklabels(), visible=False)

    axs12.set_xlabel(r"$\rm T11$",fontsize=16)
    axs12.set_ylabel(r"$ n (\rm M p c^{-3})$",fontsize=14)

    # Z09 Dam cut
    axs13 = plt.subplot(gs[13]) 
    plot_compact_sum(my_nd_Dam_Z09, colour='black', xlim_n = xlim,
                     label = r"", AX=axs13,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs13,nsum = True) 
    axs13.axhline(peak_RN_Dam,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="r" )
    
    plot_nd_3bins(my_nd_Dam_Z09,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs13)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs13)
     
    axs13.set_ylim(0.1e-6,2e-3)
    axs13.set_xlim(-0.2,0.2)
    axs13.set_yscale('log')
    plt.xticks([-0.1, 0.1],[r'$\rm c,Sph$', r'$\rm E+ES$'],fontsize=12)

    plt.setp(axs13.get_yticklabels(), visible=False)
    #plt.setp(axs13.get_xticklabels(), visible=False)
    
    axs13.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='r')

    axs13.set_xlabel(r"$\rm Z09$",fontsize=16)


    # RC15 Dam cut
    axs14 = plt.subplot(gs[14]) 
    plot_compact_sum(my_nd_Dam, colour='black', xlim_n = xlim,
                     label = r"", AX=axs14,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs14,nsum = True) 
    axs14.axhline(peak_RN_Dam,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="r" )
    
    plot_nd_3bins(my_nd_Dam,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs14)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs14)
     
    axs14.set_ylim(0.1e-6,2e-3)
    axs14.set_xlim(-0.2,0.2)
    axs14.set_yscale('log')
    axs14.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='r')
    plt.xticks([-0.1, 0.1],[r'$\rm c,Sph$', r'$\rm E+ES$'],fontsize=12)

    
    plt.setp(axs14.get_yticklabels(), visible=False)
    #plt.setp(axs14.get_xticklabels(), visible=False)
    axs14.set_xlabel(r"$\rm RC15$",fontsize=16)

    # IP13 Dam cut    plt.tight_layout()

    axs15 = plt.subplot(gs[15]) 
    plot_compact_sum(my_nd_Dam_IP13, colour='black', xlim_n = xlim,
                     label = r"", AX=axs15,nsum = True)
    plot_compact_sum(my_nd_E, colour='red', xlim_n = xlim,
                     label = r"", AX=axs15,nsum = True) 
    axs15.axhline(peak_RN_Dam,xlim[0],xlim[1],linestyle="dashed", linewidth = 3,
               color="r" )
    
    plot_nd_3bins(my_nd_Dam_IP13,'o',my_z=np.array([-0.1,-0.1,-0.1]) ,AX=axs15)
    plot_nd_3bins(my_nd_E, 's', my_z=np.array([0.1,0.1,0.1]),AX=axs15)
     
    axs15.set_ylim(0.1e-6,2e-3)
    axs15.set_xlim(-0.2,0.2)
    axs15.set_yscale('log')
    
    axs15.fill(xlim_allmass, ylim_allmass, alpha=alpha_allmass, color='r')
    plt.xticks([-0.1, 0.1],[r'$\rm c,Sph$', r'$\rm E+ES$'],fontsize=12)

    twin15 = axs15.twinx()
    #twin15.set_yscale('log')
    twin15.set_ylabel(r'$\rm Damjanov~criteria$',fontsize=14)
    #axs15.yaxis.tick_right()
    #axs15.set_yticks([])
    #axs15.yaxis.set_label_position("right")
    #axs15.set_ylabel(r'$\rm Damjanov~cut$',fontsize=16)


    axs15.set_xlabel(r"$\rm IP13$",fontsize=16)
    #twin15.set_xticks([-0.2, 0.2],['c,sph', 'E'])
    plt.setp(axs15.get_yticklabels(), visible=False)
    #plt.setp(axs15.get_xticklabels(), visible=False)  

    plt.setp(twin15.get_yticklabels(), visible=False)
    #plt.setp(twin15.get_xticklabels(), visible=False)  
    
    axs0.scatter([],[],label=r"$\rm c,~Sph~in~Bin~1 $", color ='#a5200b', marker="o")
    axs0.scatter([],[],label=r"$\rm c,~Sph~in~Bin~2 $", color ='#0b5786', marker ="o")
    axs0.scatter([],[],label=r"$\rm c,~Sph~in~Bin~3 $", color ='#2a3236', marker ="o")

    axs0.scatter([],[],label=r"$\rm E+ES~in~Bin~1 $", color ='#a5200b', marker ="s")
    axs0.scatter([],[],label=r"$\rm E+ES~in~Bin~2 $", color ='#0b5786', marker = "s")
    axs0.scatter([],[],label=r"$\rm E+ES~in~Bin~3 $", color ='#2a3236', marker = "s")
    axs0.legend(loc='upper center',fontsize=12, bbox_to_anchor=(2.2, 1.6),
          fancybox=True, shadow=False, ncol=3)
    plt.tight_layout()
    plt.show()
    return fig
    
plot_nd_3plot()
#plot_nd_E_3plot()

plot_nd_all_mass()