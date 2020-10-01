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
my_nd_Barro = np.array([4.62e-5,1.19e-4,4.3e-4])
my_nd_Dam = np.array([5.26e-5,1.39e-4,4.3e-4])
my_nd_vDokkum = np.array([1.47e-5, 1.99e-5, 9.21e-5])
my_nd_vdWel = np.array([1.89e-5,2.65e-5,9.21e-5])

my_nd_Graham = np.array([1.05e-5,1.33e-5,0])


ms0 = 12



def plot_nd_3bins(nd,marker):
    n_err_bin1 = np.sqrt(nd[0]*V1)/V1
    n_err_bin2 = np.sqrt(nd[1]*V2)/V2
    n_err_bin3 = np.sqrt(nd[2]*V3)/V3

    plt.errorbar(my_z[0],nd[0],yerr=n_err_bin1,ls='none',
                     linewidth=3, ecolor='#a5200b',zorder=20,mew=1,capsize=3)
    plt.errorbar(my_z[1],nd[1],yerr=n_err_bin2,ls='none',
                     linewidth=3, ecolor='#0b5786',zorder=20,mew=1,capsize=3)
    plt.errorbar(my_z[2],nd[2],yerr=n_err_bin3,ls='none',
                     linewidth=3, ecolor='#2a3236',zorder=20,mew=1,capsize=3) 

    plt.plot(my_z[0],nd[0], color ='#a5200b', label="Bin1",ms= ms0,marker=marker)
    plt.plot(my_z[1],nd[1], color ='#0b5786', label="Bin2",ms= ms0,marker=marker)
    plt.plot(my_z[2],nd[2], color ='#2a3236', label="Bin3",ms= ms0,marker=marker)


fig, ax = plt.subplots()
plot_nd_3bins(my_nd_Barro,'o')
plt.plot(Barro_z,Barro_nd,'--o', ms= ms0, color='g',label="Barro et al. 2013")
plt.errorbar(Barro_z,Barro_nd,yerr=Barro_nd_uerr-Barro_nd,ls='none',
                     linewidth=3, ecolor='g',zorder=20,mew=1,capsize=3) 


xedge, yedge = list(Barro_z), list(Barro_nd)
xedge.append(Barro_z[-1])
yedge.append(0)
xedge.append(Barro_z[0])
yedge.append(0)

plt.fill(xedge, yedge, alpha=0.1, color='g')
plt.xlabel("z",fontsize=18)
plt.ylabel("n ($Mpc^{-3}$)",fontsize= 18)

plt.yscale( 'log' )
plt.legend()
        
plt.show()

Barro_nd_uerr


fig, ax = plt.subplots()

plt.plot(vDokkum_z,vDokkum_nd,'--^', ms= ms0,color='y', label="van Dokkum et al. 2015")
plt.errorbar(vDokkum_z,vDokkum_nd,yerr= vDokkum_nd_uerr - vDokkum_nd,ls='none',
                     linewidth=3, ecolor='y',zorder=20,mew=1,capsize=3) 

plot_nd_3bins(my_nd_vDokkum,'^')

xedge, yedge = list(vDokkum_z), list(vDokkum_nd)
xedge.append(vDokkum_z[-1])
yedge.append(0)
xedge.append(vDokkum_z[0])
yedge.append(0)

plt.fill(xedge, yedge, alpha=0.1, color='y')
plt.xlabel("z",fontsize=18)
plt.ylabel("n ($Mpc^{-3}$)",fontsize= 18)
plt.yscale( 'log' )
plt.legend()
        
plt.show()


fig, ax = plt.subplots()

plt.plot(vdWel_z,vdWel_nd,'--s', ms= ms0, color='b', label="van der Wel et al. 2014")
plt.errorbar(vdWel_z,vdWel_nd,yerr= vdWel_nd_uerr - vdWel_nd,ls='none',
                     linewidth=3, ecolor='b',zorder=20,mew=1,capsize=3) 
plot_nd_3bins(my_nd_vdWel,'s')

xedge, yedge = list(vdWel_z), list(vdWel_nd)
xedge.append(vdWel_z[-1])
yedge.append(0)
xedge.append(vdWel_z[0])
yedge.append(0)

plt.fill(xedge, yedge, alpha=0.1, color='b')
plt.xlabel("z",fontsize=18)
plt.ylabel("n ($Mpc^{-3}$)",fontsize= 18)
plt.yscale( 'log' )
plt.legend()
        
plt.show()