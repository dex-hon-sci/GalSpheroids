#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 01:10:28 2021

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

Bin1V_file = "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin1V_Kcorr.dat"
Bin2V_file = "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin2V_Kcorr.dat"
Bin3V_file= "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin3V_Kcorr.dat"

Bin1V = SRead.read_table(Bin1V_file)
Bin2V = SRead.read_table(Bin2V_file)
Bin3V = SRead.read_table(Bin3V_file)

Bin1V_g, Bin1V_i = Bin1V[:,17], Bin1V[:,16]
Bin2V_g, Bin2V_i = Bin2V[:,17], Bin2V[:,16]
Bin3V_g, Bin3V_i = Bin3V[:,17], Bin3V[:,16]

Bin1V_g_psf, Bin1V_i_psf = Bin1V[:,21], Bin1V[:,22]
Bin2V_g_psf, Bin2V_i_psf = Bin2V[:,21], Bin2V[:,22]
Bin3V_g_psf, Bin3V_i_psf = Bin3V[:,21], Bin3V[:,22]

g = np.concatenate((Bin1V_g, Bin2V_g, Bin3V_g))
i = np.concatenate((Bin1V_i, Bin2V_i, Bin3V_i))
g_psf = np.concatenate((Bin1V_g_psf, Bin2V_g_psf, Bin3V_g_psf))
i_psf =  np.concatenate((Bin1V_i_psf, Bin2V_i_psf, Bin3V_i_psf))

BinV_file= "/home/dexter/result/stat/completeness/gal_central_psf_cpt.txt"
BinV = SRead.read_table(BinV_file)
BinV_n = SRead.read_table(BinV_file,dtype="str")

name = BinV_n[:,0]
BinV_g_psf, BinV_i_psf = BinV[:,1], BinV[:,2]
nuc_disk_class = BinV[:,3]
strong_bar_class = BinV[:,4]
nucpt_class = BinV[:,5]
mix_class = BinV[:,6]

#put the special ones into list
g_n_nucdisk, i_n_nucdisk = [], []
g_n_bar, i_n_bar = [], []
g_n_nuccpt, i_n_nuccpt = [], []
g_n_mix, i_n_mix = [], []

g_o_nucdisk, i_o_nucdisk = [], []
g_o_bar, i_o_bar = [], []
g_o_nuccpt, i_o_nuccpt = [], []
g_o_mix, i_o_mix = [], []

for k in range(len(g)):
    if BinV_g_psf[k] > 0 and nuc_disk_class[k] == 1.0:
        g_n_nucdisk.append(BinV_g_psf[k])
        i_n_nucdisk.append(BinV_i_psf[k])  
        
        g_o_nucdisk.append(g[k])
        i_o_nucdisk.append(i[k])
        
    elif BinV_g_psf[k] > 0 and strong_bar_class[k] == 1.0:
        g_n_bar.append(BinV_g_psf[k])
        i_n_bar.append(BinV_i_psf[k]) 
        
        g_o_bar.append(g[k])
        i_o_bar.append(i[k])
        
    elif BinV_g_psf[k] > 0 and nucpt_class[k] == 1.0:
         g_n_nuccpt.append(BinV_g_psf[k])
         i_n_nuccpt.append(BinV_i_psf[k])    
         
         g_o_nuccpt.append(g[k])
         i_o_nuccpt.append(i[k])
         
    elif BinV_g_psf[k] > 0 and mix_class[k] == 1.0:
        g_n_mix.append(BinV_g_psf[k])
        i_n_mix.append(BinV_i_psf[k])  
        
        g_o_mix.append(g[k]) 
        i_o_mix.append(i[k])
    else:
        pass

g_n_nucdisk, i_n_nucdisk = np.array(g_n_nucdisk), np.array(i_n_nucdisk)
g_n_bar, i_n_bar = np.array(g_n_bar), np.array(i_n_bar)
g_n_nuccpt, i_n_nuccpt = np.array(g_n_nuccpt), np.array(i_n_nuccpt)
g_n_mix, i_n_mix = np.array(g_n_mix), np.array(i_n_mix)

g_o_nucdisk, i_o_nucdisk = np.array(g_o_nucdisk), np.array(i_o_nucdisk)
g_o_bar, i_o_bar = np.array(g_o_bar), np.array(i_o_bar)
g_o_nuccpt, i_o_nuccpt = np.array(g_o_nuccpt), np.array(i_o_nuccpt)
g_o_mix, i_o_mix = np.array(g_o_mix), np.array(i_o_mix)

#rearrange the array, g and i. get rid of all the nan
g_n, i_n =[], []
for k in range(len(g)):
    if g_psf[k] > 0:
        g_n.append(g[k])
        i_n.append(i[k])        
    else:
        pass
g_n, i_n = np.array(g_n), np.array(i_n)

# remove nan, shorten the array

g_psf = g_psf[np.logical_not(np.isnan(g_psf))]
i_psf = i_psf[np.logical_not(np.isnan(i_psf))]


# convert psf from asinh magnitude to pogson magntiude
b_g, b_i = 0.9e-10, 1.8e-10

f0 = 10**(22.5/2.5) #SDSS zp = 22.5, f0 is define as the flux when conventional mag = 0

def asinh_mag_to_flux(m,b):
    return 2*b*np.sinh(((np.log(10)*m)/(-2.5))-np.log(b))

def flux_to_pogson_mag(f,zp):
    return zp-2.5*np.log10(f)


f_g_psf = asinh_mag_to_flux(g_psf,b_g)
f_i_psf = asinh_mag_to_flux(i_psf,b_i)


g_psf_n = flux_to_pogson_mag(f_g_psf,22.5)
i_psf_n = flux_to_pogson_mag(f_i_psf,22.5)

print(g_psf_n, g_psf, g_psf_n-i_psf_n)

print('median', np.median((g_psf_n-i_psf_n)-(g_psf-i_psf)))
##
s_gi = np.std(g_n-i_n) #/ np.sqrt(len(g_n))
s_gi_psf = np.std(g_psf-i_psf) #/ np.sqrt(len(g_n))

median_gi =np.median(g_n-i_n)
median_gi_psf =np.median(g_psf-i_psf)

print(s_gi, s_gi_psf)

print(g_n-i_n, g_psf-i_psf)

print("median colour",np.median(g_n-i_n))
#error = np.sqrt(s_gi**2+s_gi_psf**2)


error =np.sqrt(s_gi**2+s_gi_psf**2)
line_xy = np.array([0,0.5,1,1.5,2,2.5,3])

print('error',error)


# linear fit
def linear_func(x,m,b):
    return m*x+b

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

popt, pcov = curve_fit(linear_func, np.array(g_n-i_n), np.array(g_psf_n-i_psf_n))

perr = np.sqrt(np.diag(pcov))
print('fit_para',*popt)
print('fit_error',perr)
def plot_gi_pasf_compare():
    fig = plt.figure()

    plt.plot(line_xy,line_xy,'k--', lw=3,label=r"")
    #plt.plot([-s_gi,3],[0,3+s_gi_psf],'b--')
    #plt.plot([s_gi,3],[0,3-s_gi_psf],'b--')
    #plt.plot([-0.8665,3],[0,3+0.8665],'b--')
    #plt.plot([0.8665,3],[0,3-0.8665],'b--')
    #plt.plot([-error,3],[0,3+error],'b--', lw=3)
    #plt.plot([error,3],[0,3-error],'b--', lw=3)
    plt.plot(g_n-i_n,g_psf_n-i_psf_n,'go',ms=10,label='all')
    
    plt.plot(g_o_nucdisk-i_o_nucdisk, g_n_nucdisk-i_n_nucdisk,'bo',ms=10,label='nuc Disc')
    #plt.plot(g_o_bar-i_o_bar, g_n_bar- i_n_bar,'yo',ms=10,label='strong bar')
    plt.plot(g_o_nuccpt-i_o_nuccpt, g_n_nuccpt-i_n_nuccpt,'o',color='cyan',
             ms=10,label = 'nuc cpt')
    plt.plot(g_o_mix- i_o_mix , g_n_mix - i_n_mix ,'ko',ms=10, label='mix')

    #plt.plot(line_xy,linear_func(line_xy,*popt),'r-',lw=3,label= r'$\rm linear~fit$')
    #plt.plot(line_xy,linear_func(line_xy,popt[0]+perr[0],popt[1]+perr[1]),'r--')
    #plt.plot(line_xy,linear_func(line_xy,popt[0]-perr[0],popt[1]-perr[1]),'r--')
    
    plt.ylim(0.58,2.1)
    plt.xlim(0.58,2.1)
    plt.xlabel(r"$(g-i)_\mathrm{gal}~(\rm mag)$", fontsize=22)
    plt.ylabel(r"$(g-i)_\mathrm{psf}~(\rm mag)$", fontsize=22)
    plt.legend(loc="lower right")
    plt.show()
plot_gi_pasf_compare()

SPlot.ShowcaseCompare2.plot_compare_generic(g_n-i_n, g_psf-i_psf,
                                            para_name="(g-i)",name=[],
                                            label=["","psf"])

print('min',linear_func(1,*popt))
print('mid',linear_func(1.22,*popt))
print('max',linear_func(1.8,*popt))

