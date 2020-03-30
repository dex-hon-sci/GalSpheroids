#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 14:31:24 2020

@author: Dexter Hon
"""
#%%
"""
import packages
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import scipy as scipy
from scipy import stats

from astropy.table import Table, Column, MaskedColumn
import matplotlib.patches as mpatches


import SphRead as SRead
#%% tested
class MLRelationIband(object):
    """
    A class used to store the various Mass/Light relation

    ...

    Attributes
    ----------
    g : float
        The g-band magnitude 
    i : float
        The i-band magnitude 
    Into13_MassRatio : float
        The mass-to-light ratio (Kr) given by Into et al. 2013
    Roediger15BC03_MassRatio : float
        The mass-to-light ratio (BC03) given by Roediger et al. 2015
    Roediger15FSPS_MassRatio : float
        The mass-to-light ratio (FSPS) given by Roediger et al. 2015
    Zibetti09_MassRatio : float
        The mass-to-light ratio given by Zibetti et al. 2009
    Taylor11_MassRatio : float
        The mass-to-light ratio given by Taylor et al. 2011

    Methods
    -------
    None
    
    """
    def __init__(self,g,i):
        self.g=g
        self.i=i
        self.Into13_MassRatio = 10**(0.985*(g-i)-0.669)
        self.Roediger15BC03_MassRatio = 10**(0.979*(g-i)-0.831)
        self.Roediger15FSPS_MassRatio = 10**(0.831*(g-i)-0.597)
        self.Zibetti09_MassRatio = 10**(1.032*(g-i)-0.963)
        self.Taylor11_MassRatio =10**(0.70*(g-i)-0.68)
        
#%% tested
class MassCalculation(MLRelationIband):
    """
    A class for calculating stellar mass
    
    The class is updatable for more advance calculation

    ...

    Attributes
    ----------
    m_gal : float
        The apparant magnitude of the galaxy (in mag)
    dist : float
        The Distance of the galaxy (in Mpc)
    M_sun: float
        The reference magnitude. In this case use the relevant absolute 
        magnitude from the sun
        

    Methods
    -------
    cal_mass
        Calculate the mass with given distance and M/L ratio. 
    """
    
    def __init__(self, m_gal, dist, M_sun, *args, **kwargs):
        self._m_gal = m_gal
        self.dist = dist
        self.M_sun = M_sun
        super().__init__(*args, **kwargs)
             
    @property
    def m_gal(self):
        return(self._m_gal)

    def cal_Mass(self,ML_ratio):
        M_gal=self._m_gal-25-5*np.log10(self.dist)  
        #Absoulte Magnitude of the galaxy
        #print('M_gal',M_gal)
        #ML_ratio = getattr(self, ' ML_ratio', 0.8)
        Mass_gal = ML_ratio*(10**((self.M_sun-M_gal)/2.5))
        #print('Mass_gal',Mass_gal)
        return Mass_gal

#%% tested
class SelectionCut(object):
    """
    Class for selection cut 
    
    input by mass, generate and plot selection cut
    
    """
    """
    A class for selection cut. 
    
    It contains the functions of defining compact massive quiescent galaxies,
    and the plotting function for these criteria. 
    ...

    Attributes
    ----------
    mass : float
        The stellar mass of the galaxy 


    Methods
    -------
    cal_mass(sound=None)
        Calculate the mass with given distance and M/L ratio. 
    """
    
    def __init__(self, mass):
        self.mass = mass
        
    def Barro13_cut(self):
        cut=np.zeros(np.size(self.mass))
        i=0
        for i in range(np.size(self.mass)):
            if (np.log10(self.mass[i])>=10.0):
                cut[i] = 10**((np.log10(self.mass[i])-10.3)/1.5)
            elif (np.log10(self.mass[i])<10.0):
                cut[i] = np.nan
        return cut
    
    def vDokkum15_cut(self):
        cut=np.zeros(np.size(self.mass))
        i=0
        for i in range(np.size(self.mass)):
            if (np.log10(self.mass[i])>=10.6):
                cut[i] = 10**(np.log10(self.mass[i])-10.7)
            elif (np.log10(self.mass[i])<10.6):
                cut[i] = np.nan
        return cut

    def vdWel14_cut(self):
        cut=np.zeros(np.size(self.mass))
        i=0
        for i in range(np.size(self.mass)):
            if (np.log10(self.mass[i])>=10.7):
                cut[i] = 2.5*((self.mass[i]/1e11)**0.75)
            elif (np.log10(self.mass[i])<10.7):
                cut[i] = np.nan
        return cut

    def Graham15_broad_cut(self):
        cut=np.zeros(np.size(self.mass))
        i=0
        for i in range(np.size(self.mass)):
            if (np.log10(self.mass[i])>=10.845):
                cut[i] = 2
            elif (np.log10(self.mass[i])<10.845):
                cut[i] = np.nan
        return cut

    def plot_cut(self):

        plt.plot(self.mass,self.Barro13_cut(),"g--" , linewidth=3,
                 label="Barro et al.2013" )
        plt.vlines(1e10, 0, 10**((np.log10(1e10)-10.3)/1.5), 
                   linestyle="dashed", linewidth=3, color='g' )

        plt.plot(self.mass,self.vDokkum15_cut(),"y--" , linewidth=3, 
                 label="van Dokkum et al.2015" )
        plt.vlines(10**10.6, 0, 10**(np.log10(10**10.6)-10.7), 
                   linestyle="dashed", linewidth=3, color='y' )

        plt.plot(self.mass,self.vdWel14_cut(),"b--" , linewidth=3, 
                 label="van der Wel et al.2014" )
        plt.vlines(10**10.7, 0, 2.5*(((10**10.7)/1e11)**0.75), 
                   linestyle="dashed",linewidth=3, color='b' )


        plt.plot(self.mass,self.Graham15_broad_cut(),"k--" , linewidth=3 )
        plt.vlines(7e10, 0, 2, linestyle="dashed",linewidth=3, color='k' )

        plt.xscale( 'log' )
        plt.yscale( 'log' )
        plt.legend()
        
#%%  
class ShowcaseIndi(SelectionCut, MassCalculation):
    """
    Class for visualizing data, assuming a singular bundle input 
    
    It contains plotting methods for a population of galaxies' properties.
    
    ...

    Attributes
    ----------
    input_list : list
        The input galaxy bundle, in cpt form.
    dist: 
        A list of distance in Mpc 


    Methods
    -------
    plot_hist_percentage()
        Plot a histogram of individual galaxies with their components 
        luminosity as the length of the bar.
        Each components are stacked on top of each other with designated 
        colour code.
        
    Mass_Re_plot()
        Plot size-mass relation of the spheroid.
    """
    
    def __init__(self):
        #self.input_list = input_list
        #self.dist = dist
        super().__init__(*args, **kwargs)
        
       
    def Mass_Re_plot(x,y,name,colour,legend,alpha0): #tested
        """
        Plot the size-mass diagram given the input x (radius) and y (Mass).
        
        Parameters
        ----------
        x : numpy array
            The effective radius of a galaxy, in kpc unit.
        y : 
            The stellar mass of a galaxy, in solar mass unit.
        colour : str
            The indicator of the colour and style of the data point.
            e.g.: "bo", blue point; "r^", red triangle.
        legend : str
            The discription of the data point
        alpha0: float
            The opacity level for the highlight zone, from 0 to 1.0.
            
        ------

            
        """
        ###on off switch for names, colour
        x_edge,y_edge= [0,0,2,2], [7e10,90e11,90e11,7e10]

        plt.plot(x,y,"%s"%colour,label="%s"%legend, markersize=10, 
                 alpha=alpha0)

        #for i in range(np.size(name)):
        #    plt.text(x[i],y[i], name[i],fontsize=12)

        plt.xlabel("$Mass$ / $Mass_{sun}$",fontsize=16)
        plt.ylabel("$R_e$ (kpc)",fontsize=16)
        plt.vlines(7e10,0,2,color='k', linestyle="dashed", linewidth=3)
        plt.hlines(2,7e10,90e11, color='k', linestyle="dashed", linewidth=3)
        plt.fill(y_edge,x_edge, alpha=0.1, color='#ade0b9')

        plt.xscale( 'log' )
        plt.yscale( 'log' )
        plt.legend()
        #return fig, ax

    def plot_hist_percentage(input_list, dist): #tested
        
        """
        Plot indiividual galaxies components in histogram form.
        

        Parameters
        ----------
        input_list : list
            The input galaxy bundle, in cpt form.
        dist: 
            A list of distance in Mpc.
            
        Return
        ------
        fig, ax
            The histogram of the luminosity porportion of each galaxies.
            
        """
        category = {
                    "nucDisk" : "#165180",
                    "SecBar" : "#a88e25",
                    "Bulge" :"#ee2b10",
                    "CoreBulge" : "#b62c19",
                    "Lens": "#d67a6d",
                    "IntDisk": "#6d88d6",
                    "PrimBar": "#dad035",
                    "Disk" : "#3535da", 
                    "Ring" : "#12e3e7", 
                    "Point Source": "#3cf017"}

        category_names = list(category.keys())
        catagory_colour = list(category.values())
        
        input_list = SRead.read_list(input_list)
    
        gal_dict, gal_cpt = {},{} 
        #gal_dict with gal name as key, array of mag as data, 
        #gal_cpt same but with cpt
    
        for row in range(len(input_list)):
            label, data = [], []
            #zdist = Q[8][row]
            #D = Q[11][row]
            #print(input_list[row][0],zdist,D)
            D = dist[row]
            for item in range(len(category_names)):  #re-shuffle
                #print(input_list[row][0])
                for index in range(len(input_list[row])):
                    if input_list[row][index] == category_names[item]:  
                        #matching cpt
                    
                        #print(input_list[row][index])
                        m = input_list[row][index+2]
                        M = m - 25 -5*np.log10(D)
                        
                        Lum = 10**((M-(4.53))/(-2.5))
                        #Lum = input_list[row][index+2]
                        label.append(category_names[item])
                        data.append(Lum)
                    else:
                        pass
                gal_dict[input_list[row][0]], gal_cpt[input_list[row][0]] =\
                    data, label
        
        #print("gal_dict",gal_dict,"gal_cpt",gal_cpt)
    
        #shuffle-to-order-then-plot
        #################################################  
        #gal_labels = list(gal_dict.keys())
        #gal_data = np.array(list(gal_dict.values()))
        #gal_data_cum = gal_data.cumsum()
        
        gal_labels = list(gal_dict.keys())
        gal_data = np.array(list(gal_dict.values()))
        gal_data_cum = gal_data.cumsum()

        fig, ax = plt.subplots(figsize=(10, 20))
    
        ax.invert_yaxis()
        #ax.xaxis.set_visible(False)
        #ax.set_xlim(0, np.sum(gal_data, axis=1).max())
    
        for i in range(len(input_list)): # loop through the galaxy list
            gal_name = input_list[i][0]
            #print(gal_name)
            widths, starts = 0,0  
        
            gal_data_cum = np.array(gal_data[i]).cumsum()
            
            for j in range(len(gal_data[i])): 
                # loop through the property & cpt list
            
                widths = gal_data[i][j]                
                starts = gal_data_cum[j] - widths
            
                #print(gal_cpt[gal_name][j], category[gal_cpt[gal_name][j]],widths,starts)
                ax.barh(gal_labels[i], widths, left=starts, height=0.9, 
                        label='',#gal_cpt[gal_name][j], 
                        color=category[gal_cpt[gal_name][j]])
    
        #ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),
        #loc='lower left', fontsize='small')
        plt.xlabel( '$L_*/L_{\odot}$' ,fontsize=16)
        #plt.xscale( 'log' )
        
        for dict_row in category:
            plt.plot([],[], color = category[dict_row], linestyle='-', 
                     linewidth=13, label = dict_row)

            plt.legend(loc=1, fontsize = 13)
        return fig, ax

    def vdis_mass_plot(mass, vdis, colour, legend):
        """
        Plot central velocity dispersion versus stellar mass.
        

        Parameters
        ----------
        mass : numpy array, list
            The stellar mass of the spheroid.
            
        vdis: numpy array, list
            The central velocity dispersion
            
        colour: str, list
            The colour and the point style for the data points.
            
        legend: str, list
            The description of the 
        
        Return
        ------
        fig, ax
            The histogram of the luminosity porportion of each galaxies.
        """

            
        fig, ax = plt.subplots()        

        if type(colour) == list and type(legend) == list:

            i = 0
            for i in range(len(mass)):
                ax.plot(mass[i], vdis[i], 
                        "%s" %colour[i], label="%s" %legend[i], markersize=6)
            
        elif type(colour) == str and type(legend) == str: 
            ax.plot(mass, vdis, "%s"%colour,label="%s"%legend, markersize=6) 
        
        plt.xlabel("$M_*/M_{\odot}$",fontsize=16)
        plt.ylabel("$\sigma$ (km/s)",fontsize=16)
            
        plt.xscale( 'log' )
        plt.yscale( 'log' )
        plt.legend()
        return fig, ax

    def n_Re_plot():
        return None
    def mu0_Re_plot():
        return None
    def z_Re_plot():
        return None
        
#%%


#Mass_Re_plot_hist

#hist_4plot
#hist_4plot_norm
#hist_4plot_cum_KS
#selection_distant_Mass

#%% tested Structurally
class ShowcaseCompare2(ShowcaseIndi):

    """
    Class for visualizing data, assuming two bundle inputs
    
    It contains plotting methods for a population of galaxies' properties,
    for the sake of comparing the two data sets.
    
    ...

    Attributes
    ----------
    input_list1,2 : list
        The input galaxy bundle, in cpt form.



    Methods
    -------
    dist_compatre
        To plot the difference in distance estimation from the two input
    plot_compare_rms
        To plot the difference in rms of two decomposition
    scat_arrow
        To plot the shift in the size-mass diagram 

    """
    
    def __init__(self,input_list1,input_list2):
        #self.input_list1 = input_list1
        #self.input_list2 = input_list2
        super().__init__(*args, **kwargs)
        
    
    def plot_distdist(DD,scale,dc,sc,name,limit):
        """
        A method to plot the difference in distance estimation 
        ----------                        
        DD : 1D numpy array
            Distance by redshift independent measurement.
                
        dc : 1D numpy array         
            Distance by redshift dependent measurement.

            
        name : str
            The name of each galaxy. 
        
        scale:
    
        limit : float
            Marker of selection limit.
        
        Return
        ------
        Plot    
    
        """
        dc, DD = np.array(dc), np.array(DD)
        index = np.linspace(0.5,len(DD),len(DD))
        x_edge,y_edge= [limit,limit,120,120], \
            [min(index),max(index),max(index),min(index)]
    
        fig = plt.figure()
        gs = gridspec.GridSpec(1, 2, wspace=0.0,width_ratios=[3,1]) 
    
        # fig, axs = plt.subplot(1,2, sharey= 'row', gridspec_kw={'hspace': 0,'wspace': 0})
    
        axs1 = plt.subplot(gs[1])
        axs0 = plt.subplot(gs[0], sharey = axs1)
    
        for i in DD:
            axs0.hlines(index,0, 120, linestyle="dashed", 
                        linewidth = 1, color= 'k')
            axs1.hlines(index,0, 0.6, linestyle="dashed", 
                        linewidth = 1, color= 'k')

        axs0.plot(DD,index, "bo",label="z-independent")
        axs0.plot(dc,index, "go",label="z-dependent")
    
        axs1.plot(scale,index, "bo",label="z-independent")
        axs1.plot(sc,index, "go",label="z-dependent")

    
        axs0.vlines(limit, 0, len(DD), linestyle="dashed",
                    linewidth=3, color='black')
    
        axs0.vlines(np.average(dc), 0, len(dc), linestyle="solid",
                    linewidth=1, color='green', 
                    label='mean Dist(z-independent)')
        axs0.vlines(np.average(DD), 0, len(DD), linestyle="solid",
                    linewidth=1, color='blue', 
                    label='mean Dist(z-dependent)')

        axs0.vlines(np.average(dc)+np.std(dc), 0, len(dc), 
                    linestyle="dashed",linewidth=1, color='green')
        axs0.vlines(np.average(DD)+np.std(DD), 0, len(DD), 
                    linestyle="dashed",linewidth=1, color='blue')
    
        axs0.vlines(np.average(dc)-np.std(dc), 0, len(dc), 
                    linestyle="dashed",linewidth=1, color='green')
        axs0.vlines(np.average(DD)-np.std(DD), 0, len(DD), 
                    linestyle="dashed",linewidth=1, color='blue')
    
        axs0.fill(x_edge,y_edge, alpha=0.3, color='red',
                  label='selection limit')
    
    
        axs0.fill([np.average(DD)-np.std(DD),np.average(DD)-np.std(DD),
               np.average(DD)+np.std(DD),np.average(DD)+np.std(DD)],[min(index),
                         max(index),max(index),min(index)], alpha=0.15, 
                                        color='blue', label = '1 $\sigma$')
        axs0.fill([np.average(dc)-np.std(dc),np.average(dc)-np.std(dc),
               np.average(dc)+np.std(dc),np.average(dc)+np.std(dc)],[min(index),
                         max(index),max(index),min(index)],y_edge, alpha=0.15, 
                                        color='green',label = '1 $\sigma$')


        axs0.invert_yaxis()
        #plt.subplots_adjust(wspace = 0)
        plt.yticks(index, name)
        #ax2.yticks.set_visible(False)
        plt.setp(axs1.get_yticklabels(), visible=False)
        
        
        axs0.set_xlabel("$Distance / Mpc$",fontsize=16)
        axs1.set_xlabel("$Angular   scale$ \n  $kpc/arcsec$",fontsize=12)
        axs0.legend(loc='center left', bbox_to_anchor=(1.4, 0.5))


        return fig, axs0, axs1
    
    
    def plot_compare_rms(input_list1,input_list2): #tested
        """
        Plot comparison of rms between two list of decompistion result.
        

        Parameters
        ----------
        input_list1,2 : list
            The input galaxy bundle, in cpt form.
            
        Return
        ------
        fig, ax
            The histogram of the luminosity porportion of each galaxies.
            
        """
        fig, ax = plt.subplots()
        # same size
        input_list1 = SRead.read_list(input_list1)
        input_list2 = SRead.read_list(input_list2)
        index, rms1, rms2, gal_name= [], [], [], []
    
        for i in range(len(input_list1)):        
            if input_list1[i][0] == input_list2[i][0]: #check name consistency
                index.append(i)
                rms1.append(input_list1[i][1])
                rms2.append(input_list2[i][1])
                gal_name.append(input_list1[i][0])
            
        
        plt.hlines(np.average(rms1), 0, 104, linestyle="solid",linewidth=3, 
                   color='#7e491b')
        plt.hlines(np.average(rms2), 0, 104, linestyle="solid",linewidth=3, 
                   color='#244b87')
               
        plt.hlines(np.average(rms1)+np.std(rms1), 0, 104, linestyle="dashed",
                   linewidth=3, color='#7e491b')
        plt.hlines(np.average(rms2)+np.std(rms2), 0, 104, linestyle="dashed",
                   linewidth=3, color='#244b87')   
               
        plt.hlines(np.average(rms1)-np.std(rms1), 0, 104, linestyle="dashed",
                   linewidth=3, color='#7e491b')
        plt.hlines(np.average(rms2)-np.std(rms2), 0, 104, linestyle="dashed",
                   linewidth=3, color='#244b87')  
    
        s_rms1, s_rms2 = np.std(rms1), np.std(rms2)
        avg_rms1, avg_rms2 = np.average(rms1), np.average(rms2)
        
        #print(avg_rms1, avg_rms2,s_rms1, s_rms2)
    
        x_edge1,y_edge1= [0,0,104,104], [avg_rms1-s_rms1,avg_rms1+s_rms1,\
                         avg_rms1+s_rms1,avg_rms1-s_rms1]
        x_edge2,y_edge2= [0,0,104,104], [avg_rms2-s_rms2,avg_rms2+s_rms2,\
                         avg_rms2+s_rms2,avg_rms2-s_rms2]
           
        ax.fill(x_edge1, y_edge1, alpha=0.2, color='#e87964')
        ax.fill(x_edge2, y_edge2, alpha=0.2, color='#41a9c4')
                        
        ax.plot(index, rms1, 'o', color='#d5721d', label= "multi-cpt")
        ax.plot(index, rms2, 'o', color='#1d64d5', label="B+D")
        #j=0
        #for j in range(len(input_list1)):
        #ax.text(index[j],rms1[j], gal_name[j], fontsize=12)
        
        mark = r"$\langle \Delta rms \rangle$"
        
        ax.text(50, np.average(rms1)+0.005, 
                f"{mark} {round(avg_rms1,3)} $\pm$ {round(s_rms1,2)}",
                fontsize=14, color="black")
        ax.text(50, np.average(rms2)+0.005, 
                r"$\langle \Delta rms \rangle$ = %s $\pm$ %s"  
                %(round(avg_rms2,3),round(s_rms2,2)),fontsize=14, 
                color="black")    

        #plt.xticks(index,gal_name,rotation='vertical', fontsize)
        plt.ylabel("$\Delta$ rms", fontsize=16)
        plt.legend()   
        plt.xticks([1,50,100], [])
        return fig, ax
    
        
    def scat_arrow(x1,y1,name1,colour1,legend1,x2,y2,name2,colour2,legend2): #tested
        """
        Produce arrow plot pointing from x1, y1 to x2, y2.
        
        Parameters
        ----------
        x1, y1 : float, ndarray
            The coordinates of the origin data point
            
        name1 : str
            The name of point (x1,y1)
            
        colour1 :
            The style of the points. e.g., "bo"
            
        legend1 : str            
            The description of point (x1,y1)
        
        x2, y2 : float, ndarray
            The coordinates of the target data point
        
        name2 : str
            The name of point (x2,y1)
            
        colour1 :
            The style of the points. e.g., "bo"
            
        legend1 : str            
            The description of point (x1,y1)
            
        Return
        ------
        fig, ax
            The histogram of the luminosity porportion of each galaxies.
        """
        fig, ax = plt.subplots()

        ShowcaseIndi.Mass_Re_plot(x1,y1,name1,colour1,legend1,1.0)
        ShowcaseIndi.Mass_Re_plot(x2,y2,name2,colour2,legend2,1.0)

        N = len(x1)

        i=0
        for i in range(N):
            #ax.arrow(P2[2][i],P2[0][i], Q2[2][i]-P2[2][i],Q2[0][i]-P2[0][i],
            #head_width=1, head_length=1, fc='k', ec='k')
            arrow = mpatches.FancyArrowPatch((x1[i],y1[i]), ((x2[i]),(y2[i])),
                                             mutation_scale=15)
            ax.add_patch(arrow)
            
#############################
    def comapre_generic():
        fig,ax = plt.subplots()

        delta = para1 - para2
        
        index = np.linspace(0,len(para1),len(para1))

        
        ax.plot(delta, index, 'o')
        ax.vlines(0,min(index),max(index), linestyle="dashed",linewidth=5, color='black')
        ax.vlines(np.average(total_mag- mag_i_p),min(index),max(index), linestyle="solid",linewidth=3, color='blue')

        ax.vlines(np.average(total_mag- mag_i_p)+np.std(total_mag- mag_i_p),min(index),max(index), linestyle="dashed",linewidth=2, color='blue')
        ax.vlines(np.average(total_mag- mag_i_p)-np.std(total_mag- mag_i_p),min(index),max(index), linestyle="dashed",linewidth=2, color='blue')


        ax.invert_yaxis()
        return None
            
    
    def compare_para(list_input1,list_input2, func_index, element_index, 
                     para_name):  
        """
         Compare the parameters of the fitting components 
         Base on indexing
        ----------
        list_input1, list_input2 : list
            The python lists containing the information of ALL galaxy
                        
        func_index :
            The function of interest
            
        element_index :             
            The parameters of interest in the function
            
        para_name : 
            the parameter name
        
        Return
        ------
        delta_para
        
        """
        fig, ax = plt.subplots()
    
        list1, list2 = SRead.read_list(list_input1), SRead.read_list(list_input2)
        delta_para, index, gal_name = [], [], []
        i = 0 
        for i in range(len(list1)):
            delta_para.append(list1[i][func_index][element_index] - list2[i][func_index][element_index])
            index.append(i)
            gal_name.append(list1[i][0])
            #delta_para = np.array(delta_para)

        ax.plot(index, delta_para, 'bo')
    
        j=0
        for j in range(len(list1)):
            ax.text(index[j],delta_para[j], gal_name[j],fontsize=12)
        
        #plt.xlabel(para_name,fontsize=16)
        plt.ylabel("$\Delta$ %s" %para_name,fontsize=16)
        plt.hlines(1.0, 0, 150, linestyle="dashed",linewidth=3, color='b' )
        plt.hlines(-1.0, 0, 150, linestyle="dashed",linewidth=3, color='b' )
        return delta_para       
        
    def compare_para2(list_input1,list_input2, feature, 
                      element_index, para_name):  
        """
         Compare the parameters of the fitting components 
         Base on feature name
        ----------
        list_input1, list_input2 : list
            The python lists containing the information of ALL galaxy
                        
        func_index :
            The function of interest
            
        element_index :             
            The parameters of interest in the function
            
        para_name : 
            the parameter name
        
        Return
        ------
        delta_para
        
        """
        #######################################################################
        ## Compare the parameters of the fitting components
        ## Export an ASCII file containing the  
        ## 
        #######################################################################
        ## list_input1, 2: The python lists containing the information of ALL galaxy
        ## func_index: The function you are inteested in 
        ## feature: The component you are interested in, acceptable input: 
        ## 1) Bulge, 2)Core-depleted Bulge 3) Extended Disk, 4) Nuclear Disk, 5) Intermediate Disk, 
        ## 6) outer dearth (typeII truncation), 7) outer Surplus (typeIII truncation) 
        ## 8) Primary Bar, 9) Secondary Bar 
        ## Extra label to override the selection scheme
        #################################################################################################### 
        fig, ax = plt.subplots()
    
        list1, list2 = SRead.read_list(list_input1), SRead.read_list(list_input2)
        delta_para, index, gal_name = [], [], []
        
        #i = 0 
        for i in range(len(list1)):  #i=row
            j=0
            for j in range(len(list1[i])): #j=index and element
                
                if list1[i][j] == feature and list2[i][j] == feature:
                    delta_para.append(list1[i][j+1][element_index] - list2[i][j+1][element_index])
                    index.append(i)
                    gal_name.append(list1[i][0])
                else:
                    pass
        #delta_para = np.array(delta_para)

        #ax.plot(index, delta_para, 'bo')
        ax.plot(index, delta_para, 'bo')

        k=0
        for k in range(len(index)):
            ax.text(index[k],delta_para[k], gal_name[k],fontsize=12)
       
        print(feature, len(index))
        
        #plt.xlabel(para_name,fontsize=16)
        plt.ylabel("$\Delta$ %s" %para_name,fontsize=16)
        plt.hlines(1.0, 0, 150, linestyle="dashed",linewidth=3, color='b')
        plt.hlines(-1.0, 0, 150, linestyle="dashed",linewidth=3, color='b')
        return delta_para


    def compare_para3(list_input1,list_input2, feature, element_index, para_name1,para_name2):
        """
         Compare the parameters of the fitting components 
        ----------
        list_input1, list_input2 : list
            The python lists containing the information of ALL galaxy
                        
        func_index :
            The function of interest
            
        element_index :             
            The parameters of interest in the function
            
        para_name : 
            the parameter name
        
        Return
        ------
        delta_para
        
        """
        ####################################################################################################
        ## Compare the parameters of the fitting components
        ## Export an ASCII file containing the  
        ## 
        ####################################################################################################
        ## list_input1, 2: The python lists containing the information of ALL galaxy
        ## func_index: The function you are inteested in 
        ## feature: The component you are interested in, acceptable input: 
        ## 1) Bulge, 2)Core-depleted Bulge 3) Extended Disk, 4) Nuclear Disk, 5) Intermediate Disk, 
        ## 6) outer dearth (typeII truncation), 7) outer Surplus (typeIII truncation) 
        ## 8) Primary Bar, 9) Secondary Bar 
        ## Extra label to override the selection scheme
        #################################################################################################### 
        fig, ax = plt.subplots()
        
        list1, list2 = SRead.read_list(list_input1), SRead.read_list(list_input2)
        delta_para, index, gal_name = [], [], []
        
        #i = 0 
        for i in range(len(list1)):  #i=row
            j=0
            for j in range(len(list1[i])): #j=index and element
                
                if list1[i][j] == feature and list2[i][j] == feature:
                    delta_para.append(list1[i][j+1][element_index]/list2[i][j+1][element_index])
                    index.append(i)
                    gal_name.append(list1[i][0])
                else:
                    pass
        #delta_para = np.array(delta_para)
        
        #ax.plot(index, delta_para, 'bo')
        ax.plot(index, delta_para, 'bo')

        #k=0
        #for k in range(len(index)):
        #   ax.text(index[k],delta_para[k], gal_name[k],fontsize=12)
       
        print(feature, len(index))
        
        plt.xlabel(feature,fontsize=16)
        plt.ylabel("$log(%s/%s)$" %(para_name1,para_name2),fontsize=16)
        plt.yscale( 'log' )
        plt.xticks([1,50,100], [])
        
        return delta_para

    def compare_mag(list_input1,list_input2, element_index, para_name):
        """
         Compare the parameters of the fitting components 
        ----------
        list_input1, list_input2 : list
            The python lists containing the information of ALL galaxy
                        
        func_index :
            The function of interest
            
        element_index :             
            The parameters of interest in the function
            
        para_name : 
            the parameter name
        
        Return
        ------
        delta_para
        
        """
        list1, list2 = SRead.read_list(list_input1), SRead.read_list(list_input2)
        delta_para, index, gal_name = [], [], []
        i = 0 
        for i in range(len(list1)):
            delta_para.append(list1[i][element_index] - list2[i][element_index])
            index.append(i)
            gal_name.append(list1[i][0])
            
            ax.plot(index, delta_para, 'bo')
            
        j=0
        for j in range(len(list1)):
            ax.text(index[j],delta_para[j], gal_name[j],fontsize=12)
            
            #plt.xlabel(para_name,fontsize=16)

    
        plt.ylabel("$\Delta$ %s" %para_name,fontsize=16)
        return delta_para

#%%

#%%


# new criteria 
# look at the number of Sersic/Ferrer/Exp profile -> check which one are bulge/extended disk-> select
# what to compare major vs equvi
# Bulge, extended Disk, nuclear Disk
# Bulge = Sersic, with the higher sersic index 
## bulge_can = [array1,array2] (if  label = "sersic" or label = "CoreSersic"), bulge = array1 (if array1_n = max(array_n))
    # maybe Sph not bulge
# nuclear Disk = exp that has lowest h, usually h<7
    ## nucDisk_can = [array1,array2] (if label = "Exp", h<7 ), nucDisk = 
# intermediate Disk = exp lower h than the extended disk
    ## intDisk_can = array1 (if label = "Exp", "incl_Exp", "Broken_Exp", h<7 ) if at R_tot/2, Sersic(x) > Exp (x)
# Extended Disk = exp with the largest h 
    ##  extDisk = array1 (if label = "Exp", "incl_Exp", "Broken_Exp", h<7 ) if at R_tot/2, Sersic(x) < Exp (x) 
## primary Bar = bridge between disk and bulge; secondary Bar = internal bar form out of instability    
    ##  Bar_can = [array1, array2] (if label = "Ferrer" h<7 ), prim_bar = array1, sec_bar = array2 if r_out1 > r_out2
# Anse = the accumulation of stars builded by 
    ## Anse_can = array1 (if label = "Gauss", )
# rings/ spiral arm = 
    ## Ring_can = [array1] (if label = "Gauss", on the disk)
# 
    
# For Bars and other stuff,create a checklist of the components markers [NGC1234, Bar = 2]
# Loop through the marker table and cross matches the feature to the parameters
# override list [NGC1234, Sersic == disk], check override list first then 
    #output cpt label replace func label
    
#cpt_label = [Spheroid, NucDisk, IntDisk , ExtDisk, PrimBar, SecBar] 
    
    


#%%


#%%
#plot old area


def Mass_Re_plot_line(x,y1,y2,name,colour):
    x_edge,y_edge= [0,0,2,2], [7e10,6e11,6e11,7e10]

    yp= (y1+y2)/2.

    ax.errorbar(x,yp,xerr=[abs(y2-yp),abs(y1-yp)], fmt='|', ecolor="%s" %colour)

    for i in range(np.size(name)):
        ax.text(x[i],yp[i], name[i],fontsize=12)

    plt.ylabel("$Mass$ / $Mass_{sun}$")
    plt.xlabel("$R_e$ (kpc)")
    plt.vlines(2,7e10,6e11,color='k', linestyle="dashed", linewidth=3)
    plt.hlines(7e10,0,2, color='k', linestyle="dashed", linewidth=3)
    ax.fill(x_edge,y_edge, alpha=0.3, color='green')
    plt.legend()
#
def n_Re_plot(x,y, name ,colour,legend):
    ax.plot(x,y,"%s"%(colour) ,label="%s" %(legend) )

    #for i in range(np.size(name)):

    #    ax.text(x[i],y[i], name[i],fontsize=12)
    plt.ylabel("n",fontsize=16)
    plt.xlabel("$Mass$ / $Mass_{sun}$",fontsize=16)
    plt.xscale( 'log' )
    plt.legend()

def z_Re_plot(x,y, name ,colour,legend):
    ax.plot(x,y,"%s"%(colour) ,label="%s" %(legend) )

    #for i in range(np.size(name)):

    #    ax.text(x[i],y[i], name[i],fontsize=12)
    plt.xlabel("z",fontsize=16)
    plt.ylabel("$R_e$ (kpc)",fontsize=16)

    plt.legend()

def Mass_Re_plot_hist(x,y,name):
    x_edge,y_edge= [0,0,2,2], [7e10,2.5e11,2.5e11,7e10]

    gs = gridspec.GridSpec(4,4)

    ax_main = fig.add_subplot(gs[1:4,0:3])
    ax_xDist = fig.add_subplot(gs[0,0:3])
    ax_yDist = fig.add_subplot(gs[1:4,3])

    ax_main.scatter(x,y,marker='o')
    ax_main.set_xlabel("$Mass$ / $Mass_{sun}$",fontsize=16)
    ax_main.set_ylabel("$R_e$ (kpc)",fontsize=16)

    ax_xDist.hist(x,bins=50,align='mid')
    ax_xDist.set_ylabel('count',fontsize=16)
    ax_xCumDist = ax_xDist.twinx()
    ax_xCumDist.hist(x,bins=50,cumulative=True,histtype='step',normed=True,color='r',align='mid')
    ax_xCumDist.tick_params('y', colors='r')
    ax_xCumDist.set_ylabel('cumulative',color='r',fontsize=16)

    ax_yDist.hist(y,bins=50,orientation='horizontal',align='mid')
    ax_yDist.set_xlabel('count',fontsize=16)
    ax_yCumDist = ax_yDist.twiny()
    ax_yCumDist.hist(y,bins=50,cumulative=True,histtype='step',normed=True,color='r',align='mid',orientation='horizontal')
    ax_yCumDist.tick_params('x', colors='r')
    ax_yCumDist.set_xlabel('cumulative',color='r',fontsize=16)

    for i in range(np.size(name)):
        ax_main.text(x[i],y[i], name[i],fontsize=12)

    ax_main.axvline(x=7e10, color='k', linestyle='--')
    ax_main.axhline(y=2, color='k', linestyle='--')
    ax_main.fill(y_edge,x_edge, alpha=0.3, color='green')

#%%
################################
def hist_4plot(local_pop,distant_pop,bins, selection, out_median_local , out_median_distant, median_name_local ,median_name_distant,right, top):
    ####################################################################################################
    ## Comparison of two population in histogram form, with a certain value x, e.g. effective radius Re
    ## It seperate data into 4 sections
    ####################################################################################################
    ## local_pop: the local population data in terms of value x; input of an 4x1 numpy array #
    ## distant_pop: the distant population in terms of value x; input of an 4x1 numpy array #
    ## bins: number of bins employed for the histogram; input of an integer #
    ## selection: The selection criteria that seperate the sample; input of a list of string with 4 elements
    ## out_median_local: The array of another quantity of which the corresponding value of median(x) is of interest, e.g. redshift; input of an 1-D numpy array #
    ## out_median_distant: The array of another quantity of which the corresponding value of median(x) is of interest, e.g. redshift; input of an 1-D numpy array #
    ## median_name_local: The name of the quantity of out_median_local; input of a string
    ## median_name_distant: The name of the quantity of out_median_distant; input of a string
    ## right: x-coordinate of the text
    ## top: y-coordinate of the text
    ####################################################################################################
    fig, axs = plt.subplots(nrows=4, ncols=1,sharex=True, sharey=True, gridspec_kw={'hspace': 0})
    fig.suptitle('Population distribution')

    #index0 = np.where( distant_pop[0]== np.median(distant_pop[0]))
    #index1 = np.where( distant_pop[1]== np.median(distant_pop[1]))
   # a = np.where(distant_pop[2] == np.median(distant_pop[2]))
   # b = np.searchsorted(distant_pop[2], np.median(distant_pop[2]))

    #index2 = {True: a, False: b}[np.median(distant_pop[2]) == int(np.median(distant_pop[2]))]
    #print(np.median(distant_pop[2]), int(np.median(distant_pop[2])),np.searchsorted(distant_pop[2], np.median(distant_pop[2])))
    #index2 = np.where( distant_pop[2]== np.median(distant_pop[2]))
   # print(index2, out_median[2][index2])

    #index3 = np.where( distant_pop[3]== np.median(distant_pop[3])) #1.61e11
    index0, index0_local = np.searchsorted(distant_pop[0], np.median(distant_pop[0])), np.searchsorted(local_pop[0], np.median(local_pop[0]))
    index1, index1_local = np.searchsorted(distant_pop[1], np.median(distant_pop[1])), np.searchsorted(local_pop[1], np.median(local_pop[1]))
    index2, index2_local = np.searchsorted(distant_pop[2], np.median(distant_pop[2])), np.searchsorted(local_pop[2], np.median(local_pop[2]))
    index3, index3_local = np.searchsorted(distant_pop[3], np.median(distant_pop[3])), np.searchsorted(local_pop[3], np.median(local_pop[3]))

    local_pop_hist0 = axs[0].hist(local_pop[0],bins,facecolor='blue',label='local spheroids')
    distant_pop_hist0 = axs[0].hist(distant_pop[0],local_pop_hist0[1],facecolor='orange',alpha=0.3, label='high-z galaxies')

    axs[0].text(right, top, '%s' %selection[0], horizontalalignment='right', verticalalignment='top', fontsize=15)
    #axs[0].text(right, top-15, '0 < Dist < 1332.88 Mpc', horizontalalignment='right', verticalalignment='top')
    axs[0].text(right, top-(top/2.5), '%s = %s' %(median_name_distant, format(out_median_distant[0][index0], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)
    axs[0].text(right, top-2*(top/2.5), '%s = %s' %(median_name_local, format(out_median_local[0][index0_local], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)

    axs[1].text(right, top, '%s' %selection[1], horizontalalignment='right', verticalalignment='top' , fontsize=15)
    #axs[1].text(right, top-15, '1332.88 Mpc < Dist < 2665.75 Mpc', horizontalalignment='right', verticalalignment='top')
    axs[1].text(right, top-(top/2.5), '%s = %s' %(median_name_distant,format(out_median_distant[1][index1], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)
    axs[1].text(right, top-2*(top/2.5), '%s = %s' %(median_name_local, format(out_median_local[1][index1_local], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)

    axs[2].text(right, top, '%s' %selection[2], horizontalalignment='right', verticalalignment='top', fontsize=15)
    #axs[2].text(right, top-15, '2665.75 Mpc < Dist < 3997.75 Mpc', horizontalalignment='right', verticalalignment='top')
    axs[2].text(right, top-(top/2.5), '%s = %s' %(median_name_distant,format(out_median_distant[2][index2], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)
    axs[2].text(right, top-2*(top/2.5), '%s = %s' %(median_name_local, format(out_median_local[2][index2_local], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)

    axs[3].text(right, top, '%s' %selection[3], horizontalalignment='right', verticalalignment='top', fontsize=15)
    #axs[3].text(right, top-15, '3997.75 Mpc < Dist < 5331.50 Mpc', horizontalalignment='right', verticalalignment='top')
    axs[3].text(right, top-(top/2.5), '%s = %s' %(median_name_distant, format(out_median_distant[3][index3], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)
    axs[3].text(right, top-2*(top/2.5), '%s = %s' %(median_name_local, format(out_median_local[3][index3_local], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)

    local_pop_hist1 = axs[1].hist(local_pop[1], bins,facecolor='blue')
    distant_pop_hist1 = axs[1].hist(distant_pop[1], local_pop_hist1[1],facecolor='orange',alpha=0.3)

    local_pop_hist2 = axs[2].hist(local_pop[2], bins,facecolor='blue')
    distant_pop_hist2 = axs[2].hist(distant_pop[2], local_pop_hist2[1],facecolor='orange',alpha=0.3)

    local_pop_hist3 = axs[3].hist(local_pop[3], bins,facecolor='blue')
    distant_pop_hist3 = axs[3].hist(distant_pop[3], local_pop_hist3[1],facecolor='orange',alpha=0.3)

    plt.xlabel("Effective Radius $R_e$ (kpc)")
    plt.xlim(0,14)
    plt.yscale( 'log' )
    #plt.xscale( 'log' )

    axs[0].legend(loc='upper left')
    # Hide x labels and tick labels for all but bottom plot.
    for ax in axs:
        ax.label_outer()
    plt.show()
    
# %%
def hist_4plot_norm(local_pop,distant_pop,bins, selection,out_median_local , out_median_distant, median_name_local ,median_name_distant,right, top):
    ####################################################################################################
    ## Comparison of two population in histogram form, in a normalized form, with a certain value x, e.g. effective radius Re
    ## y-axis= y(x)/ sum(y)
    ## It seperate data into 4 sections
    ####################################################################################################
    ## local_pop: the local population data in terms of value x; input of an 4x1 numpy array #
    ## distant_pop: the distant population in terms of value x; input of an 4x1 numpy array #
    ## bins: number of bins employed for the histogram; input of an integer #
    ## selection: The selection criteria that seperate the sample; input of a list of string with 4 elements
    ## out_median_local: The array of another quantity of which the corresponding value of median(x) is of interest, e.g. redshift; input of an 1-D numpy array #
    ## out_median_distant: The array of another quantity of which the corresponding value of median(x) is of interest, e.g. redshift; input of an 1-D numpy array #
    ## median_name_local: The name of the quantity of out_median_local; input of a string
    ## median_name_distant: The name of the quantity of out_median_distant; input of a string
    ## right: x-coordinate of the text
    ## top: y-coordinate of the text
    ####################################################################################################
    fig, axs = plt.subplots(nrows=4, ncols=1,sharex=True, sharey=True, gridspec_kw={'hspace': 0})
    fig.suptitle('Probability P($R_e$)')

    index0, index0_local = np.searchsorted(distant_pop[0], np.median(distant_pop[0])), np.searchsorted(local_pop[0], np.median(local_pop[0]))
    index1, index1_local = np.searchsorted(distant_pop[1], np.median(distant_pop[1])), np.searchsorted(local_pop[1], np.median(local_pop[1]))
    index2, index2_local = np.searchsorted(distant_pop[2], np.median(distant_pop[2])), np.searchsorted(local_pop[2], np.median(local_pop[2]))
    index3, index3_local = np.searchsorted(distant_pop[3], np.median(distant_pop[3])), np.searchsorted(local_pop[3], np.median(local_pop[3]))
    
    #print(index0,index1, index2, index3)
        
    axs[0].text(right, top, '%s' %selection[0], horizontalalignment='right', verticalalignment='top', fontsize=15)
    #axs[0].text(right, top-0.08, '0 < Dist < 1332.88 Mpc', horizontalalignment='right', verticalalignment='top')
    axs[0].text(right, top-(top/2.5), '%s = %s' %(median_name_distant, format(out_median_distant[0][index0], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)
    axs[0].text(right, top-2*(top/2.5), '%s = %s' %(median_name_local, format(out_median_local[0][index0_local], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)

    axs[1].text(right, top, '%s' %selection[1], horizontalalignment='right', verticalalignment='top', fontsize=15)
    #axs[1].text(right, top-0.08, '1332.88 Mpc < Dist < 2665.75 Mpc', horizontalalignment='right', verticalalignment='top')
    axs[1].text(right, top-(top/2.5), '%s = %s' %(median_name_distant, format(out_median_distant[1][index1], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)
    axs[1].text(right, top-2*(top/2.5), '%s = %s' %(median_name_local, format(out_median_local[1][index1_local], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)

    axs[2].text(right, top, '%s' %selection[2], horizontalalignment='right', verticalalignment='top', fontsize=15)
    #axs[2].text(right, top-0.08, '2665.75 Mpc < Dist < 3997.75 Mpc', horizontalalignment='right', verticalalignment='top')
    axs[2].text(right, top-(top/2.5), '%s = %s' %(median_name_distant, format(out_median_distant[2][index2], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)
    axs[2].text(right, top-2*(top/2.5), '%s = %s' %(median_name_local, format(out_median_local[2][index2_local], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)

    axs[3].text(right, top, '%s' %selection[3], horizontalalignment='right', verticalalignment='top', fontsize=15)
    #axs[3].text(right, top-0.08, '3997.75 Mpc < Dist < 5331.50 Mpc', horizontalalignment='right', verticalalignment='top')
    axs[3].text(right, top-(top/2.5), '%s = %s' %(median_name_distant, format(out_median_distant[3][index3], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)
    axs[3].text(right, top-2*(top/2.5), '%s = %s' %(median_name_local, format(out_median_local[3][index3_local], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)

    local_pop_hist0 = np.histogram(local_pop[0],bins)
    distant_pop_hist0 = np.histogram(distant_pop[0],local_pop_hist0[1])

    local_pop_hist1 = np.histogram(local_pop[1], bins)
    distant_pop_hist1 = np.histogram(distant_pop[1], local_pop_hist1[1])

    local_pop_hist2 = np.histogram(local_pop[2], bins)
    distant_pop_hist2 = np.histogram(distant_pop[2], local_pop_hist2[1])

    local_pop_hist3 = np.histogram(local_pop[3], bins)
    distant_pop_hist3 = np.histogram(distant_pop[3], local_pop_hist3[1])

    local_N0_total, local_N1_total, local_N2_total, local_N3_total= np.sum(local_pop_hist0[0]), np.sum(local_pop_hist1[0]), np.sum(local_pop_hist2[0]), np.sum(local_pop_hist3[0])
    distant_N0_total, distant_N1_total, distant_N2_total, distant_N3_total= np.sum(distant_pop_hist0[0]), np.sum(distant_pop_hist1[0]), np.sum(distant_pop_hist2[0]), np.sum(distant_pop_hist3[0])

    local_mid_point0 = np.array([((local_pop_hist0[1][i] + local_pop_hist0[1][i+1])/2) for i in range(np.size(local_pop_hist0[1])-1)])
    local_mid_point1 = np.array([((local_pop_hist1[1][i] + local_pop_hist1[1][i+1])/2) for i in range(np.size(local_pop_hist1[1])-1)])
    local_mid_point2 = np.array([((local_pop_hist2[1][i] + local_pop_hist2[1][i+1])/2) for i in range(np.size(local_pop_hist2[1])-1)])
    local_mid_point3 = np.array([((local_pop_hist3[1][i] + local_pop_hist3[1][i+1])/2) for i in range(np.size(local_pop_hist3[1])-1)])

    distant_mid_point0 = np.array([((distant_pop_hist0[1][i] + distant_pop_hist0[1][i+1])/2) for i in range(np.size(distant_pop_hist0[1])-1)])
    distant_mid_point1 = np.array([((distant_pop_hist1[1][i] + distant_pop_hist1[1][i+1])/2) for i in range(np.size(distant_pop_hist1[1])-1)])
    distant_mid_point2 = np.array([((distant_pop_hist2[1][i] + distant_pop_hist2[1][i+1])/2) for i in range(np.size(distant_pop_hist2[1])-1)])
    distant_mid_point3 = np.array([((distant_pop_hist3[1][i] + distant_pop_hist3[1][i+1])/2) for i in range(np.size(distant_pop_hist3[1])-1)])

    local_prob_point0 = np.array([(float(local_pop_hist0[0][x])/local_N0_total) for x in range(np.size(local_pop_hist0[0]))])
    local_prob_point1 = np.array([(float(local_pop_hist1[0][x])/local_N1_total) for x in range(np.size(local_pop_hist1[0]))])
    local_prob_point2 = np.array([(float(local_pop_hist2[0][x])/local_N2_total) for x in range(np.size(local_pop_hist2[0]))])
    local_prob_point3 = np.array([(float(local_pop_hist3[0][x])/local_N3_total) for x in range(np.size(local_pop_hist3[0]))])

    distant_prob_point0 = np.array([(float(distant_pop_hist0[0][x])/distant_N0_total) for x in range(np.size(distant_pop_hist0[0]))])
    distant_prob_point1 = np.array([(float(distant_pop_hist1[0][x])/distant_N1_total) for x in range(np.size(distant_pop_hist1[0]))])
    distant_prob_point2 = np.array([(float(distant_pop_hist2[0][x])/distant_N2_total) for x in range(np.size(distant_pop_hist2[0]))])
    distant_prob_point3 = np.array([(float(distant_pop_hist3[0][x])/distant_N3_total) for x in range(np.size(distant_pop_hist3[0]))])

    axs[0].plot(local_mid_point0,local_prob_point0,color='blue', marker='o', label='local spheroids')
    axs[1].plot(local_mid_point1,local_prob_point1,color='blue', marker='o')
    axs[2].plot(local_mid_point2,local_prob_point2,color='blue', marker='o')
    axs[3].plot(local_mid_point3,local_prob_point3,color='blue', marker='o')

    axs[0].plot(distant_mid_point0,distant_prob_point0,color='orange', marker='o',linestyle='dashed', label='high-z galaxies')
    axs[1].plot(distant_mid_point1,distant_prob_point1,color='orange', marker='o',linestyle='dashed')
    axs[2].plot(distant_mid_point2,distant_prob_point2,color='orange', marker='o',linestyle='dashed')
    axs[3].plot(distant_mid_point3,distant_prob_point3,color='orange', marker='o',linestyle='dashed')

    plt.xlabel("Effective Radius $R_e$ (kpc)")

    plt.xlim(0,14)
    #plt.xscale( 'log' )

    axs[0].legend(loc='upper left')
    # Hide x labels and tick labels for all but bottom plot.
    for ax in axs:
        ax.label_outer()
    plt.show()
# %%
################################
def hist_4plot_cum_KS(local_pop,distant_pop,bins, selection, out_median, median_name,right, top):
    ####################################################################################################
    ## Comparison of two population in histogram form, with a certain value x, e.g. effective radius Re
    ## It seperate data into 4 sections
    ####################################################################################################
    ## local_pop: the local population data in terms of value x; input of an 4x1 numpy array #
    ## distant_pop: the distant population in terms of value x; input of an 4x1 numpy array #
    ## bins: number of bins employed for the histogram; input of an integer #
    ## selection: The selection criteria that seperate the sample; input of a list of string with 4 elements
    ## out_median_local: The array of another quantity of which the corresponding value of median(x) is of interest, e.g. redshift; input of an 1-D numpy array #
    ## out_median_distant: The array of another quantity of which the corresponding value of median(x) is of interest, e.g. redshift; input of an 1-D numpy array #
    ## median_name_local: The name of the quantity of out_median_local; input of a string
    ## median_name_distant: The name of the quantity of out_median_distant; input of a string
    ## right: x-coordinate of the text
    ## top: y-coordinate of the text
    ####################################################################################################
    fig, axs = plt.subplots(nrows=4, ncols=1,sharex=True, sharey=True, gridspec_kw={'hspace': 0})
    fig.suptitle('Cumulative population distribution')

    #index0 = np.where( distant_pop[0]== np.median(distant_pop[0]))
    #index1 = np.where( distant_pop[1]== np.median(distant_pop[1]))
   # a = np.where(distant_pop[2] == np.median(distant_pop[2]))
   # b = np.searchsorted(distant_pop[2], np.median(distant_pop[2]))

    #index2 = {True: a, False: b}[np.median(distant_pop[2]) == int(np.median(distant_pop[2]))]
    #print(np.median(distant_pop[2]), int(np.median(distant_pop[2])),np.searchsorted(distant_pop[2], np.median(distant_pop[2])))
    #index2 = np.where( distant_pop[2]== np.median(distant_pop[2]))
   # print(index2, out_median[2][index2])

    #index3 = np.where( distant_pop[3]== np.median(distant_pop[3])) #1.61e11
    index0, index0_local = np.searchsorted(distant_pop[0], np.median(distant_pop[0])), np.searchsorted(local_pop[0], np.median(local_pop[0]))
    index1, index1_local = np.searchsorted(distant_pop[1], np.median(distant_pop[1])), np.searchsorted(local_pop[1], np.median(local_pop[1]))
    index2, index2_local = np.searchsorted(distant_pop[2], np.median(distant_pop[2])), np.searchsorted(local_pop[2], np.median(local_pop[2]))
    index3, index3_local = np.searchsorted(distant_pop[3], np.median(distant_pop[3])), np.searchsorted(local_pop[3], np.median(local_pop[3]))

    local_pop_hist0 = axs[0].hist(local_pop[0],bins,facecolor='blue',label='local spheroids',cumulative=True, histtype="step",linewidth=3)
    distant_pop_hist0 = axs[0].hist(distant_pop[0],local_pop_hist0[1],facecolor='orange',alpha=0.3, label='high-z galaxies',cumulative=True, histtype="step",linewidth=3)

    axs[0].text(right, top, '%s' %selection[0], horizontalalignment='right', verticalalignment='top', fontsize=15)
    #axs[0].text(right, top-15, '0 < Dist < 1332.88 Mpc', horizontalalignment='right', verticalalignment='top')
    axs[0].text(right, top-(top/1.5), '%s = %s' %(median_name, format(out_median[0], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)

    axs[1].text(right, top, '%s' %selection[1], horizontalalignment='right', verticalalignment='top' , fontsize=15)
    #axs[1].text(right, top-15, '1332.88 Mpc < Dist < 2665.75 Mpc', horizontalalignment='right', verticalalignment='top')
    axs[1].text(right, top-(top/1.5), '%s = %s' %(median_name,format(out_median[1], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)

    axs[2].text(right, top, '%s' %selection[2], horizontalalignment='right', verticalalignment='top', fontsize=15)
    #axs[2].text(right, top-15, '2665.75 Mpc < Dist < 3997.75 Mpc', horizontalalignment='right', verticalalignment='top')
    axs[2].text(right, top-(top/1.5), '%s = %s' %(median_name,format(out_median[2], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)

    axs[3].text(right, top, '%s' %selection[3], horizontalalignment='right', verticalalignment='top', fontsize=15)
    #axs[3].text(right, top-15, '3997.75 Mpc < Dist < 5331.50 Mpc', horizontalalignment='right', verticalalignment='top')
    axs[3].text(right, top-(top/1.5), '%s = %s' %(median_name, format(out_median[3], '.2e')), horizontalalignment='right', verticalalignment='top', fontsize=15)

    local_pop_hist1 = axs[1].hist(local_pop[1], bins,facecolor='blue',cumulative=True,histtype="step",linewidth=3)
    distant_pop_hist1 = axs[1].hist(distant_pop[1], local_pop_hist1[1],facecolor='orange',alpha=0.3,cumulative=True,histtype="step",linewidth=5)

    local_pop_hist2 = axs[2].hist(local_pop[2], bins,facecolor='blue',cumulative=True, histtype="step",linewidth=3)
    distant_pop_hist2 = axs[2].hist(distant_pop[2], local_pop_hist2[1],facecolor='orange',alpha=0.3,cumulative=True, histtype="step",linewidth=3)

    local_pop_hist3 = axs[3].hist(local_pop[3], bins,facecolor='blue',cumulative=True, histtype="step",linewidth=3)
    distant_pop_hist3 = axs[3].hist(distant_pop[3], local_pop_hist3[1],facecolor='orange',alpha=0.3,cumulative=True, histtype="step",linewidth=3)

    plt.xlabel("Effective Radius $R_e$ (kpc)")
    plt.xlim(0,14)
    plt.yscale( 'log' )
    #plt.xscale( 'log' )

    axs[0].legend(loc='upper left')
    # Hide x labels and tick labels for all but bottom plot.
    for ax in axs:
        ax.label_outer()
    plt.show()
    
# %% 
    
def selection_distant_Mass(Mass_list,Re_list,z_list):
    z_list_1, Re_list_1 = np.array([]), np.array([])
    z_list_2, Re_list_2 = np.array([]), np.array([])
    z_list_3, Re_list_3 = np.array([]), np.array([])
    z_list_4, Re_list_4 = np.array([]), np.array([])
    i=0
    for i in range(np.size(z_list)):

        if Mass_list[i] <  6e10 and  Mass_list[i] > 4e10:
            z_list_1 = np.append(z_list_1,z_list[i])
            Re_list_1 = np.append(Re_list_1,Re_list[i])

        elif Mass_list[i] < 8e10 and Mass_list[i] > 6e10 :
            z_list_2 = np.append(z_list_2,z_list[i])
            Re_list_2 = np.append(Re_list_2,Re_list[i])

        elif Mass_list[i] < 1e11 and Mass_list[i] > 8e10 :
            z_list_3 = np.append(z_list_3,z_list[i])
            Re_list_3 = np.append(Re_list_3,Re_list[i])

        elif Mass_list[i] < 2e11 and Mass_list[i] > 1e11:
            z_list_4 = np.append(z_list_4,z_list[i])
            Re_list_4 = np.append(Re_list_4,Re_list[i])

        i=i+1

    distant_pop_Re, distant_pop_z=[],[]
    distant_pop_Re.append(Re_list_1)
    distant_pop_Re.append(Re_list_2)
    distant_pop_Re.append(Re_list_3)
    distant_pop_Re.append(Re_list_4)

    distant_pop_z.append(z_list_1)
    distant_pop_z.append(z_list_2)
    distant_pop_z.append(z_list_3)
    distant_pop_z.append(z_list_4)

    distant_pop_Re,    distant_pop_z = np.array(distant_pop_Re), np.array(distant_pop_z)
    return {"Re":distant_pop_Re,"z":distant_pop_z}


################################
