#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 00:04:02 2020

@author: dexter
"""


import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import astro_func_list as func

import numpy as np

box=[[10.8,11.0],[11.0,11.3],[11.3,11.6]]
volume = [29526.108,97704.819,357422.506]
#####
Kalvin14_all ={'M_star':[],
                 'M_star_lerr':[],
                 'M_star_uerr':[],
                 'alpha':[],
                 'alpha_lerr':[],
                 'alpha_uerr':[],
                 'phi_0':[],
                 'phi_0_lerr':[],
                 'phi_0_uerr':[],
                 'line_style':[],
                 'line_colour':[],
                 'label':[],
                 'satu':[]}

Kalvin14_morph ={'M_star':[10**(10.94), 10**(10.25),10**(10.09)],
                 'M_star_lerr':[],
                 'M_star_uerr':[],
                 'alpha':[-0.79,0.87,-0.01],
                 'alpha_lerr':[],
                 'alpha_uerr':[],
                 'phi_0':[ 0.85*1e-3,2.38*1e-3],
                 'phi_0_lerr':[],
                 'phi_0_uerr':[],
                 'line_style':['solid','solid','solid'],
                 'line_colour':['#b13322','#a2579b','#58906d'],
                 'label':['E','S0-Sa','Sab-Scd'],
                 'satu':[0.2,0.2,0.2]}

Kalvin14_dichotomy ={'M_star':[10.6,10.7],
                 'M_star_lerr':[],
                 'M_star_uerr':[],
                 'alpha':[-0.27,-1.37],
                 'alpha_lerr':[],
                 'alpha_uerr':[],
                 'phi_0':[3.96,0.42],
                 'phi_0_lerr':[],
                 'phi_0_uerr':[],
                 'line_style':['solid','solid','solid'],
                 'line_colour':[],
                 'label':['Sph','Disk'],
                 'satu':[0.2,0.2]}
#############
M = np.linspace(10**8, 10**11.6, 10**3)

# = M0*(2e30)

M_star = 10**(10.64)
alpha1, alpha2 = -0.43, -1.53
phi1_0, phi2_0 = 4.18 *1e-3, 0.74*1e-3

Phi1 = func.Schechter_func(M, alpha1, M_star, phi1_0)
Phi2 = func.Schechter_func(M, alpha2, M_star, phi2_0)
Phi = Phi1 + Phi2

M_star_E = 10**(10.94)
alpha_E= -0.79
phi_0E = 0.85 * 1e-3

Phi_E = func.Schechter_func(M, alpha_E, M_star_E, phi_0E)

M_star_S0 = 10**(10.25)
alpha_S0 = 0.87
phi_0S0 = 2.38 * 1e-3

Phi_S0 = func.Schechter_func(M, alpha_S0, M_star_S0, phi_0S0)

M_star_S = 10**(10.09)
alpha_S = -0.01
phi_0S = 3.57 * 1e-3


Phi_S = func.Schechter_func(M, alpha_S, M_star_S, phi_0S)

#############

import matplotlib.pyplot as plt

fig, ax = plt.subplots()

#plt.fill_between(M,Phi_minus,Phi_plus, color='black', alpha=0.3)

plt.plot(M,Phi,'-', color='black', label ="All")
plt.plot(M,Phi_E,'-', color = 'red', label="E")
plt.plot(M,Phi_S0,'-', color = 'purple', label= "S0-Sa")
plt.plot(M,Phi_S,'-', color = 'green', label = "Sab-cd")

plt.xscale( 'log' )
plt.yscale( 'log' )

plt.xlim(10**7.9,10**11.6)
plt.ylim(7*10**-5,3*10**-2 )
plt.legend()
plt.show()