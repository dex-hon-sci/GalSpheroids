#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 12:14:07 2021

@author: dexter

This is a script dedicated to
visualise the spheroid parameters.

This script produce the following plots:
1) mu_0 - n plots (Done) 
2) mu_0 - Re plots (Done)
3) size-mass plot with curved fit

"""
import SphSort as SSort
import SphRead as SRead
import SphPlot as SPlot
import SphAnalysis as SAna
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.style
import matplotlib as mpl
import matplotlib.gridspec as gridspec

from scipy.optimize import curve_fit

plt.style.use('classic')

mpl.rcParams['grid.linewidth'] = 1.0

## Read in the data

# Read my data
# Read the name of the ISOFIT output from a list
outlist = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_output_list_all.dat",
                               dtype='str')

# Read the geometry file for the sample
geom_file = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_geom_all.dat")
geom_file_n = SRead.read_table(
    "/home/dexter/result/stat/completeness/gal_geom_all.dat",dtype='str')

name = geom_file_n[:,0]
mu0 = geom_file[:,6]


# Read the data from the galaxy bundle
sph_mag = SRead.grab_mag("F_Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"])
Re = SRead.grab_parameter("F_Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"], 1) 
Sersic_n = SRead.grab_parameter("F_Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"], 2) 
mu_e = SRead.grab_parameter("F_Gal_bundle_equvi_cpt", ["Bulge","CoreBulge"], 0) 

core_sersic_mu_p = SRead.grab_parameter("F_Gal_bundle_equvi_cpt", ["CoreBulge"], 0)

total_mag = SRead.grab_total_mag("/home/dexter/SphProject/F_Gal_bundle_equvi_cpt")

# extrapolate the central surface brightness
R_gen = np.linspace(0,300,300*2)

# Get the distance from the parent sample
D0_all_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/vel_disp_list_all_mag_NEW.txt")

K_table = SRead.read_table(
    "/home/dexter/result/stat/completeness/diagonal_selection_bag3_Bin_all_Kcorr.dat")

mag_g, mag_i = D0_all_table[:,11], D0_all_table[:,10] #*Galaxy colours
mag_g_kcorr, mag_i_kcorr = K_table[:,19], K_table[:,18] #Galaxy colours with K correction

D, D_lerr, D_uerr = D0_all_table[:,29], D0_all_table[:,30], D0_all_table[:,31]

# Calculate the absoulte magnitude and the stellar mass
ML_select_T11 = SPlot.MLRelationIband(mag_g_kcorr,mag_i_kcorr).Taylor11_MassRatio
M = SPlot.MassCalculation(sph_mag, D, 4.53,mag_g_kcorr,mag_i_kcorr)
Abs_sph_mag = M.cal_abs_mag()


print(len(mu0), len(Sersic_n), len(Abs_sph_mag))
# Read in the data from the others.


bundle_name ="/home/dexter/SphProject/F_Gal_bundle_equvi_cpt"
bundle = SRead.read_list(bundle_name)

core = SRead.grab_parameter_whole(bundle_name, ["CoreBulge"])

##############################################
#core vs not core seperation here
S = SSort.seperator_label_generic(bundle, ["Bulge","CoreBulge"])

# Extract the info from Bulge only
sph_bulge_mag = list(SRead.grab_mag(S[0]["Bulge"], ["Bulge"]))
Re_bulge = list(SRead.grab_parameter(S[0]["Bulge"], ["Bulge"], 1) )
Sersic_n_bulge = list(SRead.grab_parameter(S[0]["Bulge"], ["Bulge"], 2))
mu_e_bulge = list(SRead.grab_parameter(S[0]["Bulge"], ["Bulge"], 0))
total_mag_bulge = list(SRead.grab_total_mag(S[0]["Bulge"]))

# cherry pick info based on the index list
name_bulge = SSort.cherry_pick(S[0]["Bulge_index"], name)
mu0_bulge = SSort.cherry_pick(S[0]["Bulge_index"], mu0)

mag_g_bulge = SSort.cherry_pick(S[0]["Bulge_index"], mag_g)
mag_i_bulge = SSort.cherry_pick(S[0]["Bulge_index"], mag_i)
mag_g_kcorr_bulge = SSort.cherry_pick(S[0]["Bulge_index"], mag_g_kcorr)
mag_i_kcorr_bulge = SSort.cherry_pick(S[0]["Bulge_index"], mag_i_kcorr)

D_bulge = SSort.cherry_pick(S[0]["Bulge_index"], D)
D_lerr_bulge =  SSort.cherry_pick(S[0]["Bulge_index"], D_lerr)
D_uerr_bulge =  SSort.cherry_pick(S[0]["Bulge_index"], D_uerr)

Abs_sph_mag_bulge = SSort.cherry_pick(S[0]["Bulge_index"], Abs_sph_mag)

# Extract the info from Core Bulge only
sph_corebulge_mag = list(SRead.grab_mag(S[1]["CoreBulge"], ["CoreBulge"]))
Re_corebulge = list(SRead.grab_parameter(S[1]["CoreBulge"], ["CoreBulge"], 1))
Sersic_n_corebulge = list(SRead.grab_parameter(S[1]["CoreBulge"], ["CoreBulge"], 2)) 
mu_e_corebulge = list(SRead.grab_parameter(S[1]["CoreBulge"], ["CoreBulge"], 0))
total_mag_corebulge = list(SRead.grab_total_mag(S[1]["CoreBulge"]))

# cherry pick info based on the index list (output list)
name_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], name)
mu0_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], mu0)

mag_g_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], mag_g)
mag_i_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], mag_i)
mag_g_kcorr_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], mag_g_kcorr)
mag_i_kcorr_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], mag_i_kcorr)

D_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], D)
D_lerr_corebulge =  SSort.cherry_pick(S[1]["CoreBulge_index"], D_lerr)
D_uerr_corebulge =  SSort.cherry_pick(S[1]["CoreBulge_index"], D_uerr)

Abs_sph_mag_corebulge = SSort.cherry_pick(S[1]["CoreBulge_index"], Abs_sph_mag)
# put the two array together  
n_combine = np.array([Sersic_n_bulge,Sersic_n_corebulge])
mu0_combine = np.array([mu0_bulge,mu0_corebulge])
Mag_combine = np.array([Abs_sph_mag_bulge,Abs_sph_mag_corebulge])

############End reading 
#%%
def plot_stack_surface_brightness_profile(r):
    bundle_name ="/home/dexter/SphProject/F_Gal_bundle_equvi_cpt"
    
    bundle=SRead.read_list(bundle_name)
    
    fig = plt.figure()

    #plot individual curve
    for i in range(len(bundle)):        
        for j in range(len(bundle[i])):
            identifier = bundle[i][j]
            if identifier == "Bulge":
                para = bundle[i][j+1]
                line = SAna.AnalyticFunctions.mu_sersic_func(r,*list(para))
                plt.plot(r,line,"r--",lw =3)
                pass
            elif identifier == "CoreBulge":
                para = bundle[i][j+1]
                line = SAna.AnalyticFunctions.mu_core_sersic_func(r,*list(para))
                plt.plot(r,line,"k-",lw= 3)
                pass
    plt.gca().invert_yaxis()
    #plt.ylim(np.log10(24),np.log(10))
    plt.ylim(24,12)
    plt.xlim(0,120)
    plt.show()
    
    return fig

#%%
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
    n_line = np.linspace(-5, 20, 100)
    mu_line = np.linspace(5, 30, 100)

    #Fitting the relevant data
    popt_mu,pcov_mu = curve_fit(SAna.AnalyticFunctions.linear_func_1D, 
                                mu0[fit_instruc], Mag[fit_instruc])
    popt_n,pcov_n = curve_fit(SAna.AnalyticFunctions.linear_func_1D, 
                              np.log10(n[fit_instruc]), Mag[fit_instruc])
  
    #Plotting
    fig = plt.figure()
    gs = gridspec.GridSpec(ncols=2, nrows=1,
                               hspace=0, wspace=0.0) 
    # Define colour and markers
    markers = ["o","X", "s","+","*"]
    colour = ["k","#d20d0d","#ebb800","#0f920f","#0f4592"]
    
    # The Mag vs Sersic n plot
    axt0 = plt.subplot(gs[0])
    axt0.plot(10**n_line, 
              SAna.AnalyticFunctions.linear_func_1D(n_line,*popt_n),'b--',lw=4)
    # Check the dimension of the input array 
    if n.shape[0] != len(label):
        axt0.plot(n,Mag,'ko',ms=10)
    elif n.shape[0] == len(label):
        for i in range(len(n)):
            axt0.plot(n[i],Mag[i], linestyle="None",
                      marker=markers[i], color = colour[i],
                      ms=10,label=label[i])
    
    #axt0.legend(loc="lower right")
    axt0.set_ylabel(r"$\mathfrak{M}_{i}$", fontsize=22)
    axt0.set_xlabel(r"$\mathrm{Sersic}~n $", fontsize=22)
    axt0.set_xscale('log')
    axt0.set_xlim(3.1*10**-1,20)
   
    # The Mag vs mu_0 plot
    axt1 = plt.subplot(gs[1],sharey=axt0)
    axt1.plot(mu_line, 
              SAna.AnalyticFunctions.linear_func_1D(mu_line,*popt_mu),'b--',lw=4)
   
    # Check the dimension of the input array  
    if mu0.shape[0] != len(label):
        axt1.plot(mu0,Mag,'ko',ms=10)
    elif mu0.shape[0] == len(label):
        for i in range(len(mu0)):
            axt1.plot(mu0[i],Mag[i], linestyle="None",
                      marker=markers[i],color = colour[i],
                      ms=10,label=label[i])    
    print("Linear fit mu",*popt_mu)
    print("Linear fit n",*popt_n)
    
    axt1.legend(loc="lower right")
    axt1.set_xlabel(r"$ \mu_\mathrm{0,i}$", fontsize=22)
    axt1.set_xlim(19.2,12.2)

    plt.setp(axt1.get_yticklabels(), visible=False)
    plt.ylim(-25.5,-15)
    plt.gca().invert_yaxis()
    plt.show()
    return (fig,*popt_n,*popt_mu)
#%%
# plot the Mag vs n and mu0 plot
#plot_n_mu0_Mag_2plot(n,mu0,Mag,label=[r"$type~1$",r"$type~2$"])

#plot_n_mu0_Mag_2plot(Sersic_n,mu0,Abs_sph_mag)
B = plot_n_mu0_Mag_2plot(n_combine,mu0_combine,
                     Mag_combine,label=[r"$\rm S\'{e}rsic$",
                                        r"$\rm Core-S\'{e}rsic$"])
#,label=[r"$type~1$",r"$type~2$"])
print(B)
#plot_stack_surface_brightness_profile(R_gen)