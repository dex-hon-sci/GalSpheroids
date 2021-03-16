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
from scipy.signal import find_peaks
import matplotlib.image as mpimg

import SphRead as SRead
import SphSort as SSort


__all__=["ImageProcessing","MLRelationIband","MassCalculation", "SelectionCut", 
         "ShowcaseIndi","ShowcaseCompare2", "PlotHist","Plot2D"]

__author__="Dexter S.H. Hon"

fontname_choice = "Times New Roman"
## make histogram class inheriting Shocaseindi and 
#%%
class ImageProcessing():
    
    def __init__(self,a):
        self._a = a


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
        
    #@property
    #def dist(self):
    #    return(self.dist)
        
    #@property
    #def M_sun(self):
    #    return(self.M_sun)
    def cal_abs_mag(self):
        """
        Calculate the absoulte magntiude of a galaxy

        Returns
        -------
        None.

        """
        M_gal=self._m_gal-25-5*np.log10(self.dist)  
        
        return M_gal

    def cal_Mass(self,ML_ratio):
        """
        A method to calculate the stellar mass of a galaxy
        
        Parameters
        ----------
        ML_ratio : float
            The mass-light ratio of the galaxy given by the class
            MLRelationIband.

            
        Return
        ------
        """
        
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
    A class for selection cut. 
    
    It contains the functions of defining compact massive quiescent galaxies,
    and the plotting function for these criteria. 
    ...

    Attributes
    ----------
    mass : float
        The stellar mass of the galaxy 
    Dist: float

    Methods
    -------
    Barro13_cut
    """
    
    def __init__(self, mass, Dist):
        self.mass = mass
        self.Dist = Dist
        
    def parent_sample_cut(self):
        cut=np.zeros(np.size(self.Dist))
        i=0
        for i in range(np.size(self.Dist)):
            if ((self.Dist[i])<=115.0):
                cut[i] = 10**((0.01414*self.Dist[i])+10)
            elif ((self.Dist[i])>115.0):
                cut[i] = np.nan
        return cut

        
    def Barro13_cut(self):
        cut=np.zeros(np.size(self.mass))
        i=0
        a = 10.0 #mass limit
        for i in range(np.size(self.mass)):
            if (np.log10(self.mass[i])>=a):
                cut[i] = 10**((np.log10(self.mass[i])-10.3)/1.5)
            elif (np.log10(self.mass[i])<a):
                cut[i] = np.nan
        return cut
    
    def vDokkum15_cut(self):
        cut=np.zeros(np.size(self.mass))
        i=0
        a = 10.6 # mass limit
        for i in range(np.size(self.mass)):
            if (np.log10(self.mass[i])>=a):
                cut[i] = 10**(np.log10(self.mass[i])-10.7)
            elif (np.log10(self.mass[i])<a):
                cut[i] = np.nan
        return cut

    def vdWel14_cut(self):
        cut=np.zeros(np.size(self.mass))
        i=0
        a = 10.7 #mass limit
        for i in range(np.size(self.mass)):
            if (np.log10(self.mass[i])>=a):
                cut[i] = 2.5*((self.mass[i]/1e11)**0.75)
            elif (np.log10(self.mass[i])<a):
                cut[i] = np.nan
        return cut
    
    def Cassata11_cut(self):
        cut=np.zeros(np.size(self.mass))
        i=0
        a = 10.0 #mass limit
        for i in range(np.size(self.mass)):
            if (np.log10(self.mass[i])>=a):
                cut[i] = 10**((np.log10(self.mass[i])*0.4)-5.5)
            elif (np.log10(self.mass[i])<a):
                cut[i] = np.nan
        return cut
    
    def Damjanov14_cut(self):
        cut=np.zeros(np.size(self.mass))
        i=0
        a = 10.0 # mass limit
        for i in range(np.size(self.mass)):
            if (np.log10(self.mass[i])>=a):
                cut[i] = 10**((np.log10(self.mass[i])*0.568)-5.74)
            elif (np.log10(self.mass[i])<a):
                cut[i] = np.nan
        return cut
    
    def Graham15_broad_cut(self):
        cut=np.zeros(np.size(self.mass))
        i=0
        a = 10.845
        for i in range(np.size(self.mass)):
            if (np.log10(self.mass[i])>=a):
                cut[i] = 2
            elif (np.log10(self.mass[i])<a):
                cut[i] = np.nan
        return cut
    

    def plot_cut(self):
        """
        Plot all th cut avaliable in this class.

        Returns
        -------
        None.

        """
        
        plt.plot(self.mass,self.Cassata11_cut(),
                 ls = "dashed", color="cyan", linewidth=3,
                 label="Cassata et al. 2011" )
        plt.plot(self.mass,self.Damjanov14_cut(),"r--" , linewidth=3,
                 label="Damjanov et al. 2014" )

        
        plt.plot(self.mass,self.Barro13_cut(),"g--" , linewidth=3,
                 label="Barro et al. 2013" )
        plt.vlines(1e10, 0, 10**((np.log10(1e10)-10.3)/1.5), 
                   linestyle="dashed", linewidth=3, color='g' )

        plt.plot(self.mass,self.vDokkum15_cut(),"y--" , linewidth=3, 
                 label="van Dokkum et al. 2015" )
        plt.vlines(10**10.6, 0, 10**(np.log10(10**10.6)-10.7), 
                   linestyle="dashed", linewidth=3, color='y' )

        plt.plot(self.mass,self.vdWel14_cut(),"b--" , linewidth=3, 
                 label="van der Wel et al. 2014" )
        plt.vlines(10**10.7, 0, 2.5*(((10**10.7)/1e11)**0.75), 
                   linestyle="dashed",linewidth=3, color='b' )


        plt.plot(self.mass,self.Graham15_broad_cut(),"k--" , linewidth=3 )
        plt.vlines(7e10, 0, 2, linestyle="dashed",linewidth=3, color='k' )

        plt.xscale( 'log' )
        plt.yscale( 'log' )
        plt.legend()
        
    
    def plot_cut_specific(self,input_cut,label,alpha0= 0.2, AX=plt):
        """
        Plot a specific cut in this class.

        Parameters
        ----------
        input_cut : str
            The keyword for the selection cut. The options are:
            "Barro": Barro
        AX:
            
            
        Returns
        -------
        None.

        """
        a = 1e10 # initial guess of the left edge
        
        
        
        if input_cut == "Barro":
            AX.plot(self.mass,self.Barro13_cut(),"g--" , linewidth=3,
                 label= label )
            AX.vlines(1e10, 0, 10**((np.log10(1e10)-10.3)/1.5), 
                   linestyle="dashed", linewidth=3, color='g' )
            
            # define the left edge (mass limit)
            a = 1e10
            xedge, yedge = [5e12,a], [1e-3,1e-3]
            #define the upper edge
            xedge.append(a)
            yedge.append(10**((np.log10(a)-10.3)/1.5))
                         
            xedge.append(xedge[0])
            yedge.append(10**((np.log10(xedge[0])-10.3)/1.5))
            
            AX.fill(xedge, yedge, alpha=0.1, color='g')

            
        elif input_cut == "vDokkum":
            AX.plot(self.mass,self.vDokkum15_cut(),"y--" , linewidth=3, 
                 label=label )
            AX.vlines(10**10.6, 0, 10**(np.log10(10**10.6)-10.7), 
                   linestyle="dashed", linewidth=3, color='y' )
            
            # define the left edge (mass limit)
            a = 10**10.6
            xedge, yedge = [5e12,a], [1e-3,1e-3]
            #define the upper edge
            xedge.append(a)
            yedge.append(10**(np.log10(a)-10.7))
                         
            xedge.append(xedge[0])
            yedge.append(10**(np.log10(xedge[0])-10.7))
                         
            AX.fill(xedge, yedge, alpha=0.1, color='y')

            
        elif input_cut == "vdWel":
            AX.plot(self.mass,self.vdWel14_cut(),"b--" , linewidth=3, 
                 label=label )
            AX.vlines(10**10.7, 0, 2.5*(((10**10.7)/1e11)**0.75), 
                   linestyle="dashed",linewidth=3, color='b' )
            
            # define the left edge (mass limit)
            a = 10**10.7
            xedge, yedge = [5e12,a], [1e-3,1e-3]
            #define the upper edge
          
            xedge.append(a)
            yedge.append(2.5*(((a)/1e11)**0.75))
                         
            xedge.append(xedge[0])
            yedge.append(2.5*(((xedge[0])/1e11)**0.75))
                         
            AX.fill(xedge, yedge, alpha=0.1, color='b')

            
        elif input_cut == "Damjanov":
            AX.plot(self.mass,self.Damjanov14_cut(),"r--" , linewidth=3,
                 label=label )
            AX.vlines(1e10, 0, 10**((np.log10(1e10)*0.568)-5.74), 
                   linestyle="dashed",linewidth=3, color='r' )         
            # define the left edge (mass limit)
            a = 1e10
            xedge, yedge = [5e12,a], [1e-3,1e-3]
            xedge.append(a)
            yedge.append(10**((np.log10(a)*0.568)-5.74))
                         
            xedge.append(xedge[0])
            yedge.append(10**((np.log10(xedge[0])*0.568)-5.74))
                         
            AX.fill(xedge, yedge, alpha=0.1, color='r')

        elif input_cut == "Graham":
            AX.plot(self.mass,self.Graham15_broad_cut(),"k--" , linewidth=3,
                     label=label  )
            AX.vlines(7e10, 0, 2, linestyle="dashed",linewidth=3, color='k' )
            
            # define the left edge (mass limit)
            a = 7e10
            xedge, yedge = [5e12,a], [1e-3,1e-3]       
            xedge.append(a)
            yedge.append(2)
                         
            xedge.append(xedge[0])
            yedge.append(2)
            
            AX.fill(xedge, yedge, alpha=0.1, color='k')

            
        else:
            raise KeyboardInterrupt 
        
        AX.set_xscale( 'log' )
        AX.set_yscale( 'log' )
        AX.legend()
        
        
    def selection_subsample(self, input_list, 
                            direction="down"):
        """
        A generic method to select a set of sample base on a cut.
        
        Parameters
        ----------
        input_list: list
        
        
        Returns
        -------
        Bag: Dict
            {index:[1,4,5],
             bag:[[],[]...]}

        """
        index_list, bag_list=[],[]
        func = self.parent_sample_cut
        
        if direction == "up":
            for i in range(len(input_list)):
                if self.mass[i] > func[i]:
                    index_list.append(i)
                    bag_list.append(input_list[i])
                    
        elif direction == "down":
            for i in range(len(input_list)):
                if self.mass[i] < func[i]:
                    index_list.append(i)
                    bag_list.append(input_list[i])
        
        Bag = {'index':index_list,
               'bag':bag_list}
        return Bag  
    
    
#%%
class ShowcaseIndi(SelectionCut, MassCalculation):
    """
    Class for visualizing data, assuming a singular bundle input.
    
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
    show_name():
        Show the name of the data point in the plot.
        
    err_I_ratio_plot:    
    
    Mass_Re_plot():
        Plot size-mass relation of the spheroid.
        
    plot_hist_percentage():
        Plot a histogram of individual galaxies with their components 
        luminosity as the length of the bar.
        Each components are stacked on top of each other with designated 
        colour code.
        
    vdis_mass_plot():
    """
    
    def __init__(self):
        #self.input_list = input_list
        super().__init__(*args, **kwargs)
        
    def show_name(x,y,name,size=12,A=plt):
        """
        Show the name of a data in a plot

        Parameters
        ----------
        x, y : list, array
            The x,y coordinate of the point.

        name : TYPE
            The name of the data point.

        Returns
        -------
        None.

        """
        for j in range(len(x)):
            A.text(x[j],y[j], name[j], fontsize=size)
            
    def err_I_ratio_plot(list1, cutlist=None):
        """
        Plot the ratio of error to pixel value as a function of radius.
        It reads a list of output file from ISOFIT.
        The function extract the SMA, NPIX_DATA, VAR and I2 from the dat file.
        
        error = sigma/sqrt(N)
        

        Parameters
        ----------
        list1 : str
            ASCII filname for a list of ISOFIT output file.

        Returns
        -------
        plot.

        """
        #data extraction part
        o_list = SRead.read_table(list1,dtype="str")
                
        m_SMA, m_error_I = [], []    
        
        for row in range(len(o_list)):
            isofit_output = SRead.read_table(o_list[row])
            
            sigma = np.array(isofit_output[:,3])
            N = np.array(isofit_output[:,34])
            I = np.array(isofit_output[:,1])
            
            SMA = np.array(isofit_output[:,0])*0.4
            
            SMA = SMA/max(SMA)
            
            #print((sigma/np.sqrt(N))/I)  
            Ierr = np.array(isofit_output[:,2])
            error_I = Ierr/ I          
            #error_I = np.sqrt(sigma/N)
                    
            m_SMA.append(SMA)
            m_error_I.append(error_I)
            
        final_storage = {"SMA": m_SMA,
                         "error_I": m_error_I}    
        # loop through, cut base if SMA> cut value
        
        if cutlist == None:
            cutlist=np.zeros(np.size(final_storage['SMA']))
        else:
            pass
            
        cut_list = SRead.read_table(cutlist)
        cut_value = cut_list[:,12]
        
        m_SMA2, m_error_I_2 = [],[]
        for row in range(len(m_SMA)):
            
            m_index = []
            for index in range(len(m_SMA)):
                #print(m_SMA[index],cut_value[row])
                if m_SMA[row][index] > cut_value[row]:
                    m_index.append(index)
                else:
                    pass
                SMA_2 = np.delete(m_SMA[row],m_index)
                error_I_2 = np.delete(m_error_I[row],m_index)
                
            #print(SMA_2[0], m_SMA[row][0])
                
            m_SMA2.append(SMA_2)
            m_error_I_2.append(error_I_2)
        
        for i in range(len(m_SMA2)):
            #print((m_SMA[i]), (m_error_I[i]))

            #print(np.size(m_SMA[i]), np.size(m_error_I[i]))
            #plt.plot(m_SMA2[i], m_error_I_2[i], color='blue', linestyle="solid",
            #         linewidth=0.4)
            plt.plot(m_SMA2[i], m_error_I_2[i], '-bo', ms= 0.5, alpha=0.3)
            plt.plot(m_SMA2[i], m_error_I_2[i]*-1.0, '-bo', ms= 0.5, alpha=0.3)

            file_name =i

        plt.xlabel("$R/ R_{max}$",fontsize=20,fontname = fontname_choice)
        plt.ylabel("$I_{err}(R)/I$",fontsize=20, fontname = fontname_choice)
            #plt.yscale( 'log' )
            #plt.xscale( 'log' )
        plt.xlim(0,1)
        plt.ylim(-2.0,2.0)

        plt.hlines(0,0,200, color='k', linestyle="dashed", linewidth=3)
        #plt.legend()
        plt.savefig("./indi_err_img2/overall_err.png", dpi=300)
        plt.show()
#            plt.savefig("./indi_err_img2/%s.png"%file_name, dpi=200)
#            plt.close()
            #plt.plot(m_SMA[i], -1.0*m_error_I[i],'b-')
            #plt.fill(m_SMA[i], x_edge, alpha=0.1, color='#ade0b9')
 
        #plt.xlabel("$R/\,arcsec$",fontsize=16)
        #plt.ylabel("$I_{err}(R)/I$",fontsize=16)
        ##plt.yscale( 'log' )
        ##plt.xscale( 'log' )


        #plt.hlines(0,0,200, color='k', linestyle="dashed", linewidth=3)
        ##plt.legend()

        
        return None
    
    def Mass_Re_plot(x,y,xerr=None,yerr=None,name=None,
                     ms= 4 , colour='#e6a31e',
                     legend=None,alpha0=None,
                     marker ='o',lw= 12): #tested
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
            
        Return
        ------
        The size-mass plot.
            
        """
        #x_edge,y_edge= [0,0,2,2], [7e10,90e11,90e11,7e10]
        #plt.plot(x,y,fmt='o')
        
        plt.plot(x,y, marker, color=colour,label=legend,markersize=ms, 
                 alpha=alpha0+0.35)
        plt.errorbar(x,y,xerr=xerr,yerr=yerr,ls='none',linewidth=lw,
                  ecolor=colour,capsize=0,
                  alpha=alpha0
                  )
        #ShowcaseIndi.show_name(x, y, name)



        #plt.xlabel("$M_{*}$ / $M_{\odot}$",fontsize=16)
        #plt.ylabel("$R_e$ (kpc)",fontsize=16)
        
        #plt.vlines(7e10,0,2,color='k', linestyle="dashed", linewidth=3)
        #plt.hlines(2,7e10,90e11, color='k', linestyle="dashed", linewidth=3)
       #plt.fill(y_edge,x_edge, alpha=0.1, color='#ade0b9')
        plt.xlim(1.87e9,5e12)
        plt.ylim(0,300)

        plt.xscale( 'log' )
        plt.yscale( 'log' )
        
        plt.legend()
        #return fig, ax

    def plot_hist_percentage(input_list, dist, norm=False,scale='log'): #tested
        
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
    

        plt.xlabel( '$L_*/L_{\odot}$' ,fontsize=16)
        plt.xscale( scale )
        
        for dict_row in category:
            plt.plot([],[], color = category[dict_row], linestyle='-', 
                     linewidth=13, label = dict_row)
            ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),
                  loc='lower left', fontsize='small')
            #plt.legend(loc=1, fontsize = 13)
        return fig, ax

    def vdis_mass_plot(mass, vdis, colour, legend): #tested
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

    #
    def mass_function_plot(mass,box, volume=None, colour=None, label=None, 
                           trim=True, plot_yes=True):
        """
        Plot the mass function 

        Parameters
        ----------
        mass : 1-D numpy array
            DESCRIPTION.
        box : list
            Box parameters. It contains the boundaries of the bin.
            e.g. [[0,10],[10,100],[100,200]]
                3 bins, from 0-10, 10-100 and 100-200
        bg : list, optional
            Background curve. 
            
            The default is None.
        volume: float
            The volume to a bin, in Mpc^{3}
        Returns
        -------
        None.

        """
        # Grouping
        master_bin= []
        n=0
        for n in range(len(box)):
            box_bin=[]
            for i in range(np.size(mass)):
                if mass[i] > box[n][0] and mass[i] < box[n][1]:
                    box_bin.append(mass[i])                    
                else:
                    pass
            box_bin = np.array(box_bin)
            print(n,'number of sampe in the box', np.size(box_bin))

            #if box_bin.size == 0:
            #    master_bin.append(np.array([np.nan]))
            #else:
            master_bin.append(box_bin) #output an 
        #print('master_bin', master_bin)
        
        # calculating
        if volume == None:
            vol = 1
        else:
            vol = volume
            

        # the number density, default volume to be 1 Mpc^{-3} dex^{-1}
        
        dex_factor = np.array([(box[j][1]-box[j][0]) for j in range(len(box))])
        dex_factor=np.average(dex_factor)
        
        
        # The number of sample within each bin
        N = [np.size(master_bin[j]) for j in range(len(master_bin))]
        N_err = [(np.sqrt(np.size(master_bin[j]))/vol)/dex_factor 
                   for j in range(len(master_bin))]
        
        # The number density of each bin
        nu_dens = [(np.size(master_bin[j])/vol)/dex_factor 
                   for j in range(len(master_bin))]
        nu_dens = np.array(nu_dens)
        
        # the mid point for each x-axis bin
        mid_pt =  [10**(sum(box[j])/2) for j in range(len(box))]
        mid_pt = np.array(mid_pt)
        
        #Trim away the zeros, and the corresponding mid_pt
        if trim == True:
            zero_i = [i for i in range(len(list(nu_dens))) if nu_dens[i] == 0]

            nu_dens = SSort.trim_value(nu_dens, 0)
            N = SSort.trim_value(N, 0)
            N_err = SSort.trim_value(np.array(N_err), 0)

            mid_pt = SSort.remove_value(mid_pt, zero_i)
        elif trim == False:
            pass
        else:
            raise Exception("Trim option has to be either ture of false")
        
        # Assume Poisson error sqrt(N)/vol/dex_factor
        nu_dens_err = N_err
            
        
        N =  [np.sqrt(np.size(master_bin[j])/vol/dex_factor) for j in 
                   range(len(master_bin))]
        N = np.array(N)

        #plotting 
        #print('dex_factor', dex_factor)
        #print('nudens',nu_dens)
        #print('nudens_err',nu_dens_err)
        #print('mid_pt',mid_pt)
        if plot_yes == True:
            plt.plot(mid_pt, nu_dens, 'o', color = colour, ms = 14, label=label)
            plt.errorbar(mid_pt,nu_dens,yerr=nu_dens_err,ls='none',
                     linewidth=3, ecolor=colour,zorder=20,mew=1,capsize=3)        
        
            plt.xlabel("$M_*/M_{\odot}$",fontsize=16)
            plt.ylabel("$\Phi(Mpc^{-3}dex^{-1})$",fontsize=16)
            
            plt.xscale( 'log' )
            plt.yscale( 'log' )
        elif plot_yes == False:
            pass
        else:
            raise Exception("You need to decide either show plot or not")
        return nu_dens, mid_pt

    def cpt_to_total_by_type_plot(bundle, morph_name, morph, cpt=["Bulge","CoreBulge"], 
                                  show_plot=True, 
                                  mode="average",AX=plt):
        """
        Plotting the magntiude cpt-to-total ratio

        Parameters
        ----------
        bundle : str
            The name of the input galaxy bundle.
        morph : TYPE
            DESCRIPTION.
        cpt : str, optional
            DESCRIPTION. The default is "Bulge".
        show_plot : bol, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        The dictionary sortted by morphology.

        """
        # input: Bundle, what component, morphology type
        name = SRead.grab_name(bundle)
        cpt_mag = SRead.grab_mag(bundle, cpt)
        total_mag = SRead.grab_total_mag(bundle)
        
        # calculare flux ratio
        mag_ratio = 10**((cpt_mag-total_mag) / 2.5)
        
        
        #mag_ratio = 10**(cpt_mag)/10**total_mag

        mag_ratio = np.log10(mag_ratio)
        # check name alignment
        if len(morph_name) == len(name):
            for i in range(len(name)):
                #print(morph_name[i],name[i])

                if morph_name[i] == name[i]:
                    pass
                else:
                    raise Exception("name doesn't match!")
                
        
        # Bin by morphology
        type_dict = ["E","0","S"]
        morph_list = ["EAS","EABS","EBS","SA0","SAB0","SB0","SA","SAB","SB"]
        
        gal_bin = []
        
        # loop for morphology type

        for morph_name in morph_list:
            morph_bin =[] # bin that store based on morph type
            for i in range(len(morph)):
                if morph[i] == morph_name:
                    morph_bin.append(mag_ratio[i])
            
            gal_bin.append(np.array(morph_bin))
            #gal_bin = np.array(gal_bin)
        
        mag_dict={"EAS": gal_bin[0], "EABS": gal_bin[1], "EBS": gal_bin[2],
                  "SA0": gal_bin[3],"SAB0": gal_bin[4],"SB0": gal_bin[5],
                  "SA" : gal_bin[6], "SAB": gal_bin[7], "SB": gal_bin[8]}
                                   
        average_bin, std_bin = [], []        
        
        if mode=="average":
            for i in range(len(morph_list)):
                avg_ratio = np.average(gal_bin[i]), 
                std_ratio = np.std(gal_bin[i])/np.sqrt(len(gal_bin[i]))
                average_bin.append(avg_ratio)
                std_bin.append(std_ratio)
                print(average_bin)
        else:
            pass
        average_bin = np.array(average_bin)
        std_bin = np.array(std_ratio)
        # plot it        
        
        x_index =[]
        
        for i in range(len(mag_dict.keys())):
            b = [i]*len(gal_bin[i])
            x_index.append(np.array(b))
            
        # make x axis
        #AX.plot(x_index)
        
        #plot the individual galaxies
        
        if show_plot == True:
            for i in range(len(x_index)):
                AX.plot(x_index[i],gal_bin[i],'o',ms=15,color="#905a97",alpha =0.4
                    )
            AX.plot([],[],'o',ms=15,color="#905a97",alpha =0.4, label="galaxies")
        #print(list(np.linspace(0,len(mag_dict.keys()),num=len(mag_dict.keys()))))
        
            AX.plot([0,1,2,3,4,5,6,7,8],average_bin, 'o',color="k", ms = 12,
                label = "average")
            AX.errorbar([0,1,2,3,4,5,6,7,8], average_bin, yerr = std_bin, 
                     ms = 12, linewidth=3, ls='none',
                     color='k',zorder=20,mew=1,capsize=3)
        
            AX.set_ylabel(r"$\rm log_{10}(B/T)$",fontsize=22)
        #AX.set_xlabel(r"$ \rm Morphology$",fontsize=18)
        #AX.set_yscale( 'log' )
            AX.set_xlim(-1,9)
            AX.set_ylim(-0.2,2)
        #AX.set_xticks(mag_)
            plt.xticks([0,1,2,3,4,5,6,7,8], list(mag_dict.keys()), fontsize=18, 
                   rotation=70)
            plt.legend(loc=1)
            #[0,1,2,3,4,5,6,7,8], 
            #          list(mag_dict.keys()))
        else:
            pass

        return mag_dict
        
    def n_Re_plot():
        return None
    def mu0_Re_plot():
        return None
    def z_Re_plot():
        return None
        
#%% tested Structurally
class ShowcaseCompare2(ShowcaseIndi):

    """
    Class for visualizing data, assuming two bundle inputs.
    
    It contains plotting methods for a population of galaxies' properties,
    for the sake of comparing the two data sets.
    
    ...

    Attributes
    ----------
    input_list1,2 : list
        The input galaxy bundle, in cpt form.



    Methods
    -------
    plot_distdist
        To plot the difference in distance estimation from the two input
    plot_compare_rms
        To plot the difference in rms of two decomposition
    scat_arrow
        To plot the shift in the size-mass diagram 
    comapre_generic
    
    plot_seperation_generic
    """
    
    def __init__(self,input_list1,input_list2):
        #self.input_list1 = input_list1
        #self.input_list2 = input_list2
        super().__init__(*args, **kwargs)
        
    
    def plot_distdist(DD,scale,dc,sc,name, l_limit, u_limit, 
                      DD_err=None,dc_err=None):
        """
        A method to plot the difference in distance estimation 
        ----------                        
        DD : 1D numpy array
            Comoving distance by redshift independent measurement.
        DD_err:
            
        scale: 1D numpy array
            The kpc/arcsec scale by redshift independent measurement.
            
        dc : 1D numpy array         
            Distance by redshift dependent measurement.

        sc: 1D numpy array
            The kpc/arcsec scale by redshift dependent measurement.

            
        name : str
            The name of each galaxy. 
            
        l_limit: float
            Marker of selection lower limit.
        
        u_limit : float
            Marker of selection upper limit.
        
        Return
        ------
        Plot    
    
        """
        dc, DD = np.array(dc), np.array(DD)
        DD_err = np.array(DD_err)
        index = np.linspace(0,len(DD),num=len(DD),dtype=int)
        
        x_edge_l,y_edge_l = [0,0,l_limit,l_limit], \
            [min(index),max(index),max(index),min(index)]        
        x_edge_u,y_edge_u = [u_limit,u_limit,120,120], \
            [min(index),max(index),max(index),min(index)]
    
        fig = plt.figure()
        gs = gridspec.GridSpec(1, 2, wspace=0.0,width_ratios=[3,1]) 

        # fig, axs = plt.subplot(1,2, sharey= 'row', gridspec_kw={'hspace': 0,'wspace': 0})
    
        axs1 = plt.subplot(gs[1])
        axs0 = plt.subplot(gs[0], sharey = axs1)
    
        for i in DD:
            axs0.hlines(index,0, 120, linestyle="dashed", 
                        linewidth = 0.5, color= 'k')
            axs1.hlines(index,0, 0.6, linestyle="dashed", 
                        linewidth = 0.5, color= 'k')

        axs0.plot(DD,index, "bo",label="z-independent")
        axs0.errorbar(DD,index,xerr=DD_err,yerr=None,ls='none',linewidth=1,
                     ecolor='b',zorder=20,mew=1,capsize=3)
        
        axs0.plot(dc,index, "go",label="z-dependent")
        #axs0.errorbar(dc,index,xerr=dc_err,yerr=None,ls='none',linewidth=1,
        #             ecolor='g',zorder=20,mew=1,capsize=3)

        axs1.plot(scale,index, "bo",label="z-independent")
        axs1.plot(sc,index, "go",label="z-dependent")

    
        axs0.vlines(u_limit, 0, len(DD), linestyle="dashed",
                    linewidth=3, color='black')
        axs0.vlines(l_limit, 0, len(DD), linestyle="dashed",
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
    
        axs0.fill(x_edge_u,y_edge_u, alpha=0.3, color='red',
                  label='selection limit')
        axs0.fill(x_edge_l,y_edge_l, alpha=0.3, color='red',
                  label='selection limit')    
    
        axs0.fill([np.average(DD)-np.std(DD),np.average(DD)-np.std(DD),
               np.average(DD)+np.std(DD),np.average(DD)+np.std(DD)],
                  [min(index),max(index),max(index),min(index)], alpha=0.15, 
                                        color='blue', label = '1 $\sigma$')
        axs0.fill([np.average(dc)-np.std(dc),np.average(dc)-np.std(dc),
               np.average(dc)+np.std(dc),np.average(dc)+np.std(dc)],
                  [min(index), max(index),max(index),min(index)], alpha=0.15, 
                                        color='green',label = '1 $\sigma$')


        axs0.invert_yaxis()
        #plt.subplots_adjust(wspace = 0)
        plt.yticks(index, name)
        #ax2.yticks.set_visible(False)
        plt.setp(axs1.get_yticklabels(), visible=False)
        
        
        axs0.set_xlabel("$Distance/ \, Mpc$",fontsize=12)
        axs1.set_xlabel("Angular \n  scale \n $kpc/arcsec$",fontsize=12)
        axs0.legend(loc='center left', bbox_to_anchor=(1.4, 0.5))


        return fig, axs0, axs1
    
    
    def plot_distdist_3points(DD,dc,d, name, l_limit, u_limit, 
                              DD_err=None,dc_err=None, d_err=None, 
                               d_spec = None, d_spec_name = None , d_spec_err =None,
                               decision=[]):
        """
        A method to plot the difference in distance estimation 
        ----------                        
        DD : 1D numpy array
            Comoving distance by redshift dependent measurement. Usually 
            reserved for Mould 2000 distance
            
            
        dc : 1D numpy array         
            Distance by redshift dependent measurement. Usually reserved for 
            Willick 1997 distance

        d: 1D numpy array
            Distance by redshift dependent measurement. Usually reserved for 
            Cosmicflow-3 distance

        name : str
            The name of each galaxy. 
        
        limit : float
            Marker of selection limit.
        
        
        Optional
        --------
        DD_err: 1D numpy array
            default: None
        
        dc_err: 1D numpy array
            default: None

        d_err: 1D numpy array
            default: None
            
        d_spec: 1D numpy array
            Special distance. Usually reserved for the z-independent distance.
            default: None
            
        d_spec_name: 1D numpy array
            The name of the galaxy with assigned special distance. 
            Usually reserved for the z-independent distance.
            default: None
            
        d_spec_err: 1D numpy array
            The error for the special distance. 
            Usually reserved for the z-independent distance.
            default: None
            
        decision: list
            A list of str containing the decision for the chosen distance 
            source.
            default: []

        Return
        ------
        Plot    
    
        """    
        d, dc, DD = np.array(d), np.array(dc), np.array(DD)
        d_err, dc_err, DD_err = np.array(d_err), np.array(dc_err), np.array(DD_err)
        
        decision_color = {"Cosmicflow-3":"#ab005a" ,
                          "Mould2000": "blue",
                          "Willick1997": "green",
                          "z-independent":"#d79734",
                          "out_of_bound":"black",
                          "out_of_bound(1)":"black",
                          "out_of_bound(2)":"black",
                          "o":"black"}
        ############################create d_spec####
        #lookup name -> make index table - > plot on the right index
        d_spec_index = []
        d_spec_d = []
        d_spec_d_err = []
        
        for i in range(len(d_spec_name)):
            for j in range(len(name)):
                if d_spec_name[i] == name[j]: 
                    print(name[j])
                    d_spec_index.append(float(j))
                    d_spec_d.append(d_spec[i])
                    d_spec_d_err.append(d_spec_err[i])
                    
        print(d_spec_index, d_spec_d, d_spec_d_err)       
        #############################################
        index = np.linspace(0,len(DD),len(DD),dtype=int)
        
        x_edge_l,y_edge_l = [0,0,l_limit,l_limit], \
            [min(index)-10,max(index)+10,max(index)+10,min(index)-10]        
        
        x_edge_u,y_edge_u = [u_limit,u_limit,120,120], \
            [min(index)-10,max(index)+10,max(index)+10,min(index)-10]
    
        fig = plt.figure()
        gs = gridspec.GridSpec(1, 1) 

        # fig, axs = plt.subplot(1,2, sharey= 'row', gridspec_kw={'hspace': 0,'wspace': 0})
    
        axs0 = plt.subplot(gs[0])
    
        #for i in DD:
        #    axs0.hlines(index,0, 120, linestyle="dashed", 
        #                linewidth = 0.5, color= 'k')

        ms0, lw = 12, 4        
        
        axs0.plot(dc,index, "o", ms = ms0, color="green",
                  label="Willick et al. 1997")
        axs0.errorbar(dc,index,xerr=dc_err,yerr=None,ls='none',linewidth=lw,
                     ecolor='g',zorder=20,mew=1,capsize=ms0)

        axs0.plot(DD,index, "o", ms = ms0, color="blue",
                  label="Mould et al. 2000")
        axs0.errorbar(DD,index,xerr=DD_err,yerr=None,ls='none',linewidth=lw,
                     ecolor='b',zorder=20,mew=1,capsize=ms0)



        axs0.plot(d,index, "o", ms = ms0,  color="#ab005a",
                  label="Cosmicflow-3")
        axs0.errorbar(d,index,xerr=d_err,yerr=None,ls='none',linewidth=lw,
                     ecolor='#ab005a',zorder=20,mew=1,capsize=ms0)
        
        
        axs0.plot(d_spec_d,d_spec_index, "o", ms = ms0, color="#d79734" ,
                  label="z-independent")
        axs0.errorbar(d_spec_d,d_spec_index,xerr=d_spec_d_err,yerr=None,
                      ls='none',linewidth=lw, ecolor='#d79734', zorder=20,
                      mew=1,capsize=ms0)
        
        axs0.vlines(l_limit, 0 -10, len(DD)+10, linestyle="dashed",
                    linewidth=5, color='black')       
        
        axs0.vlines(u_limit, 0 -10, len(DD)+10, linestyle="dashed",
                    linewidth=5, color='black')
        
#        axs0.vlines(np.average(d), 0, len(d), linestyle="solid",
#                    linewidth=1, color='#ab005a', 
#                    label=' ')
#        axs0.vlines(np.average(d)+np.std(d), 0, len(d), 
#                    linestyle="dashed",linewidth=1, color='#ab005a')
#        axs0.vlines(np.average(d)-np.std(d), 0, len(d), 
#                    linestyle="dashed",linewidth=1, color='#ab005a')    
#
#        axs0.vlines(np.average(DD), 0, len(DD), linestyle="solid",
#                    linewidth=1, color='blue', 
#                    label='mean Dist(z-dependent)')
#        axs0.vlines(np.average(DD)+np.std(DD), 0, len(DD), 
#                   linestyle="dashed",linewidth=1, color='blue')
#        axs0.vlines(np.average(DD)-np.std(DD), 0, len(DD), 
#                    linestyle="dashed",linewidth=1, color='blue')
#
#           
#        axs0.vlines(np.average(dc), 0, len(dc), linestyle="solid",
#                    linewidth=1, color='green', 
#                    label='mean Dist(z-independent)')
#        axs0.vlines(np.average(dc)+np.std(dc), 0, len(dc), 
#                    linestyle="dashed",linewidth=1, color='green')
#        axs0.vlines(np.average(dc)-np.std(dc), 0, len(dc), 
#                    linestyle="dashed",linewidth=1, color='green')
#
#    
        axs0.fill(x_edge_l,y_edge_l, alpha=0.3, color='red',
                  label='selection limit')

        axs0.fill(x_edge_u,y_edge_u, alpha=0.3, color='red',
                  label='selection limit')
#    
#        axs0.fill([np.average(d)-np.std(d),np.average(d)-np.std(d),
#               np.average(d)+np.std(d),np.average(d)+np.std(d)],
#                  [min(index),max(index),max(index),min(index)], alpha=0.05, 
#                                        color='#ab005a')
#        
#        axs0.fill([np.average(DD)-np.std(DD),np.average(DD)-np.std(DD),
#               np.average(DD)+np.std(DD),np.average(DD)+np.std(DD)],
#                  [min(index),max(index),max(index),min(index)], alpha=0.05, 
#                                        color='blue')
#        
#        axs0.fill([np.average(dc)-np.std(dc),np.average(dc)-np.std(dc),
#               np.average(dc)+np.std(dc),np.average(dc)+np.std(dc)],
#                  [min(index), max(index),max(index),min(index)], alpha=0.05, 
#                                        color='green')
        


        ## Write the name of the distance source
        #for i in range(len(decision)):
        #    for j in  decision_color.keys():
        #        if decision[i] == j:
        #            color_choice = decision_color[j]
        #            axs0.text(limit+10, index[i], decision[i], fontsize=12, 
        #                      color = color_choice)                       
            
          #      else: 
           #         color_choice = "black"
          #          axs0.text(limit+10, index[i], decision[i], fontsize=12, 
          #                    color = color_choice)


        axs0.invert_yaxis()
        plt.subplots_adjust(wspace = 0)
        plt.yticks(index, name)
        #ax2.yticks.set_visible(False)
        
        ##################
        #make the decision making on the right
        ##################
        axs0.set_xlabel(r"$Distance/ \, M p c$",fontsize=22)
        axs0.legend(loc='upper left')


        return fig, axs0

    def plot_dist_difference3_summary(DD,dc,d, name, 
                              DD_err=None,dc_err=None, d_err=None, 
                               d_spec = None, d_spec_name = None , d_spec_err =None,
                               decision=[]):
        """
        Same as above, but plotting all three bins in the same time.

        Parameters
        ----------
        DD : TYPE
            DESCRIPTION.
        dc : TYPE
            DESCRIPTION.
        d : TYPE
            DESCRIPTION.
        name : TYPE
            DESCRIPTION.
        l_limit : TYPE
            DESCRIPTION.
        u_limit : TYPE
            DESCRIPTION.
        DD_err : TYPE, optional
            DESCRIPTION. The default is None.
        dc_err : TYPE, optional
            DESCRIPTION. The default is None.
        d_err : TYPE, optional
            DESCRIPTION. The default is None.
        d_spec : TYPE, optional
            DESCRIPTION. The default is None.
        d_spec_name : TYPE, optional
            DESCRIPTION. The default is None.
        d_spec_err : TYPE, optional
            DESCRIPTION. The default is None.
        decision : TYPE, optional
            DESCRIPTION. The default is [].

        Returns
        -------
        None.

        """

        d, dc, DD = np.array(d), np.array(dc), np.array(DD)
        d_err, dc_err, DD_err = np.array(d_err), np.array(dc_err), np.array(DD_err)
        
        decision_color = {"Cosmicflow-3":"#ab005a" ,
                          "Mould2000": "blue",
                          "Willick1997": "green",
                          "z-independent":"#d79734",
                          "out_of_bound":"black",
                          "out_of_bound(1)":"black",
                          "out_of_bound(2)":"black",
                          "o":"black"}
        ############################create d_spec####
        #lookup name -> make index table - > plot on the right index
        d_spec_index = []
        d_spec_d = []
        d_spec_d_err = []
        
        for i in range(len(d_spec_name)):
            for j in range(len(name)):
                if d_spec_name[i] == name[j]: 
                    print(name[j])
                    d_spec_index.append(float(j))
                    d_spec_d.append(d_spec[i])
                    d_spec_d_err.append(d_spec_err[i])
                    
        print(d_spec_index, d_spec_d, d_spec_d_err)       
        #############################################
        index = np.linspace(0,len(DD),len(DD),dtype=int)
        
    
        fig = plt.figure()
        gs = gridspec.GridSpec(1, 1) 

    
        axs0 = plt.subplot(gs[0])
    


        ms0, lw = 12, 4        
        
        axs0.plot(dc,index, "o", ms = ms0, color="green",
                  label="Willick et al. 1997")
        axs0.errorbar(dc,index,xerr=dc_err,yerr=None,ls='none',linewidth=lw,
                     ecolor='g',zorder=20,mew=1,capsize=ms0)

        axs0.plot(DD,index, "o", ms = ms0, color="blue",
                  label="Mould et al. 2000")
        axs0.errorbar(DD,index,xerr=DD_err,yerr=None,ls='none',linewidth=lw,
                     ecolor='b',zorder=20,mew=1,capsize=ms0)



        axs0.plot(d,index, "o", ms = ms0,  color="#ab005a",
                  label="Cosmicflow-3")
        axs0.errorbar(d,index,xerr=d_err,yerr=None,ls='none',linewidth=lw,
                     ecolor='#ab005a',zorder=20,mew=1,capsize=ms0)
        
        
        axs0.plot(d_spec_d,d_spec_index, "o", ms = ms0, color="#d79734" ,
                  label="z-independent")
        axs0.errorbar(d_spec_d,d_spec_index,xerr=d_spec_d_err,yerr=None,
                      ls='none',linewidth=lw, ecolor='#d79734', zorder=20,
                      mew=1,capsize=ms0)

        # Define highlighted boxes
        
        #x_edge1,y_edge1= [75,75,110,110], [5e11,max_mass,max_mass,5e11]
        #x_edge2,y_edge2= [45,45,75,75], [2e11,max_mass,max_mass,2e11]
        #x_edge3,y_edge3= [0,0,45,45], [1e11,max_mass,max_mass,1e11]
        
        x_edge_1,y_edge_1 = [75,75,110,110], \
            [min(index)-10,max(index)+10,max(index)+10,min(index)-10]        
        
        x_edge_2,y_edge_2 = [45,45,75,75], \
            [min(index)-10,max(index)+10,max(index)+10,min(index)-10]

        x_edge_3,y_edge_3 = [0,0,45,45], \
            [min(index)-10,max(index)+10,max(index)+10,min(index)-10]        

        x_edge_limit,y_edge_limit = [110,110,120,120], \
            [min(index)-10,max(index)+10,max(index)+10,min(index)-10]         

        
    
        
        axs0.vlines(110, 0 -10, len(DD)+10, linestyle="dashed",
                    linewidth=5, color='black')
        
        axs0.vlines(75, 0 -10, len(DD)+10, linestyle="dashed",
                    linewidth=5, color='black')
        axs0.vlines(45, 0 -10, len(DD)+10, linestyle="dashed",
                    linewidth=5, color='black')           
        
        axs0.fill(x_edge_1,y_edge_1, alpha=0.3, color='green',
                  label='Bin 1')

        axs0.fill(x_edge_2,y_edge_2, alpha=0.2, color='green',
                  label='Bin 2')
        
        axs0.fill(x_edge_3,y_edge_3, alpha=0.1, color='green',
                  label='Bin 3')

        axs0.fill(x_edge_limit,y_edge_limit, alpha=0.3, color='red',
                  label='selection limit')
    
        axs0.invert_yaxis()
        plt.subplots_adjust(wspace = 0)
        plt.yticks(index, name)
        
        ##################
        #make the decision making on the right
        ##################
        axs0.set_xlabel(r"$\rm Distance/ \, M p c$",fontsize=22)
        axs0.legend(loc='upper left')


        return fig, axs0   
    
    
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
        plt.ylabel(r"$\rm \Delta$ rms", fontsize=16)
        plt.legend()   
        plt.xticks([1,50,100], [])
        return fig, ax
    
        
    def plot_scat_arrow(x1,y1,name1,colour1,legend1, 
                        x2,y2,name2,colour2,legend2): #tested
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
                                             mutation_scale=4)
            ax.add_patch(arrow)
            
    def plot_seperation_generic(para1, para2, limit, para_name="", 
                                colour1="blue", colour2="orange", name=[], 
                                label=[1,2]):
        """
        A generic method to plot the seperation between 
        the same parameter from different source.
        ----------
        para1: list, numpy array
            parameter 1
            
        para2: list, numpy array
            parameter 2
            
        limit:
        
        Optional
        ---------
        para_name: str
        
        colour: str

        name: list
            A list of galaxy names.
            
        colour: str
            The colour of the lines.
            
        name: list

        label: list
        
        Return
        ------
        delta: numpy array
            parameter1 - parameter2
            
        """
        fig, ax = plt.subplots()
        
        index = np.linspace(0,len(para1),len(para1))
        

        
        min_index, max_index = min(index), max(index)
        
        avg_para1, std_para1 = np.average(para1), np.std(para1)
        avg_para2, std_para2 = np.average(para2), np.std(para2)


        range_can = [avg_para1+std_para1, avg_para1-std_para1, 
                     avg_para2+std_para2, avg_para2-std_para2]
        
        para = np.concatenate([para1,para2])
        
        #for i in range(len(index)):
        #    ax.hlines(index, 2*min(para), 2*max(para), linestyle="dashed", linewidth = 0.2, color= 'k')
        
        
        
        ax.plot(para1,index,'o',color = colour1, label=label[0])
        ax.plot(para2,index,'o',color = colour2, label=label[1])



        ax.vlines(limit,min_index,max_index,linestyle="dashed",
                  linewidth=3, color="black", label= "limit")
        
        ax.vlines(avg_para1,min_index,max_index,linestyle="solid",
                  linewidth=3, color=colour1, label= "mean %s"%para_name)
        ax.vlines(avg_para1+std_para1,min_index,max_index, 
                  linestyle="dashed",linewidth=2, color=colour1)
        ax.vlines(avg_para1-std_para1,min_index,max_index, 
                  linestyle="dashed",linewidth=2, color=colour1)
        
        
        ax.vlines(avg_para2,min_index,max_index,linestyle="solid",
                  linewidth=3, color=colour2, label= "mean %s"%para_name)
        ax.vlines(avg_para2+std_para2,min_index,max_index, 
                  linestyle="dashed",linewidth=2, color=colour2)
        ax.vlines(avg_para2-std_para2,min_index,max_index, 
                  linestyle="dashed",linewidth=2, color=colour2)
        
        
        
        ax.invert_yaxis()
        #ax.set_ytick(index,name)
        
        plt.yticks(index, name)
        #ax2.yticks.set_visible(False)
        plt.setp(ax.get_yticklabels(), visible=False)
        ax.set_xlabel(para_name, fontsize=16)
        plt.legend()
        return fig, ax
        
        
    def plot_compare_generic(para1, para2, sub=True , div=False,
                            para_name="", colour="blue", 
                            name=None, label=[1,2]): 
        #tested
        """
        A generic method to plot the difference between the same parameter from
        different source.
        ----------
        para1: list, numpy array
            parameter 1
        para2: list, numpy array
            parameter 2

        
        Optional
        ---------
        sub: bool
            Subtraction indicator, default: True
        div: bool
            Division indicator, default: False
        para_name: str
            The name of the parameter
        colour: str
            The colour of the lines.
        name: list
            A list of galaxy names.
        
        Return
        ------
        delta: numpy array
            parameter1 - parameter2
            
        """
        
        fig,ax = plt.subplots()

        delta_sub = para1-para2
        delta_div = para1/para2
        
        if sub == div:
            raise Exception("Mode has to be either subtraction or division")
        elif sub == True:
            delta = delta_sub
            mid_line = 0
            xlabel = "$ %s_{%s} \,-\, %s_{%s}$"  %(para_name, label[0] , 
                                                   para_name, label[1])
            
        elif div == True:
            delta = delta_div
            mid_line = 1.0
            xlabel = "$ %s_{%s} \,/\, %s_{%s}$"  %(para_name, label[0] , 
                                                   para_name, label[1])
            
        else: 
            raise Exception("Impossible, you will \
                            never reach this error message.")

        index = np.linspace(0,len(para1),len(para1))
        min_index, max_index = min(index),max(index)
        
        
        avg_delta = np.average(delta)
        std_delta = np.std(delta)
        
        ax.plot(delta, index, 'o')
        
        ax.vlines(mid_line,min_index,max_index, linestyle="dashed",
                  linewidth=5, color='black')
        ax.vlines(avg_delta,min_index,max_index,linestyle="solid",
                  linewidth=3, color=colour, label= "mean %s"%para_name)

        ax.vlines(avg_delta+std_delta,min_index,max_index, 
                  linestyle="dashed",linewidth=2, color=colour)
        ax.vlines(avg_delta-std_delta,min_index,max_index, 
                  linestyle="dashed",linewidth=2, color=colour)

        x_edge,y_edge=  [avg_delta-std_delta, avg_delta+std_delta,\
                         avg_delta+std_delta,avg_delta-std_delta], \
                         [min_index,min_index,max_index,max_index]
      
        
        ax.fill(x_edge, y_edge, alpha=0.1, color=colour, label="$1\sigma$")
        
        ax.invert_yaxis()
        
        #plt.setp(ax.get_yticklabels(), visible=False)
        
        plt.yticks(index, name)

        ax.set_xlabel(xlabel, fontsize=16)
        ax.legend(loc='center left', bbox_to_anchor=(0.5, 0.5))
        return fig, ax, delta
    
    def plot_compare_index_para(list_input1, list_input2, 
                                  func_index, number, sub =True , div = False, 
                                  para_name="para", colour="green", name=[],
                                  label=[1,2]):
        #tested
        """
        Compare the parameters of the fitting components base on indexing
        ----------
        list_input1, list_input2: str
            The galaxy bundle containing the information of ALL galaxy


        func_index: float
            The index of the function.
            
        number: float
            The index of the parameter.
        
        Optional
        ---------
        sub: bool
            Subtraction indicator, default: True
        div: bool
            Division indicator, default: False
        para_name: str
            The name of the parameter
        colour: str
            The colour of the lines.
        name: list
            A list of galaxy names.
        
        Return
        ------
        plot:
            
            
        """
        list1 = SRead.read_list(list_input1)
        list2 = SRead.read_list(list_input2)
        
        
        para1, para2, gal_name = [], [], []
        
        
        for i in range(len(list1)):
            para1.append(list1[i][func_index][number]) 
            para2.append(list2[i][func_index][number])

            gal_name.append(list1[i][0])
            #delta_para = np.array(delta_para)
        para1 = np.array(para1)
        para2 = np.array(para2)

        
        plot = ShowcaseCompare2.plot_compare_generic(para1, para2, 
                                                     sub=sub , div=div, 
                                                     para_name=para_name, 
                                                     colour=colour,
                                                     name=gal_name,
                                                     label = label)
        return plot

    def plot_compare_feature_para(list_input1, list_input2, 
                                  keyword, number, sub =True , div = False, 
                                  para_name="para", colour="green", name=[], 
                                  label=[1,2]):
        #tested
        """
        Compare the parameters of the fitting components, base on feature.
        ----------
        list_input1, list_input2 : list
            The galaxy bundle containing the information of ALL galaxy
                        
        keyword : list
            The feature name.
            e.g. ["Bulge","CoreBulge"]
        
        number: float
            The index of the parameter .
            
        Optional
        ---------
        sub: bool
            Subtraction indicator, default: True
        div: bool
            Division indicator, default: False
        para_name: str
            The name of the parameter
        colour: str
            The colour of the lines.
        name: list
            A list of galaxy names.
        
        Return
        ------
        plot
        
        """
        para1 = SRead.grab_parameter(list_input1, keyword, number)
        para2 = SRead.grab_parameter(list_input2, keyword, number)
        
        plot = ShowcaseCompare2.plot_compare_generic(para1, para2, 
                                                     sub=sub , div=div, 
                                                     para_name=para_name, 
                                                     colour=colour,name=name,
                                                     label = label)
        return plot

    def plot_compare_feature_mag(list_input1, list_input2, 
                                  keyword, sub =True, div = False, 
                                  para_name="mag", colour="green", name=[],
                                  label = [1,2]):
        #tested
        """
        Compare the magnitude of the components, base on feature.
        ----------
        list_input1, list_input2 : list
            The galaxy bundle containing the information of ALL galaxy
                        
        keyword : list
            The feature name.
            e.g. ["Bulge","CoreBulge"]
            
            
        Optional
        ---------
        sub: bool
            Subtraction indicator, default: True
        div: bool
            Division indicator, default: False
        para_name: str
            The name of the parameter
        colour: str
            The colour of the lines.
        name: list
            A list of galaxy names.
        
        Return
        ------
        plot
        
        """
        mag1 = SRead.grab_mag(list_input1, keyword)
        mag2 = SRead.grab_mag(list_input2, keyword)
        
        plot = ShowcaseCompare2.plot_compare_generic(mag1, mag2, 
                                                     sub=sub , div=div, 
                                                     para_name=para_name, 
                                                     colour=colour,name=name,
                                                     label =label)
        return plot
    
    def plot_compare_total_mag(list_input1,list_input2, 
                               sub =True, div = False, 
                               para_name="Total mag", colour="green", 
                               name=[], label=[1,2]):
        #tested
        """
        Compare the total magnitude of the galaxy..
        ----------
        list_input1, list_input2 : list
            The galaxy bundle containing the information of ALL galaxy
                                    
            
        Optional
        ---------
        sub: bool
            Subtraction indicator, default: True
        div: bool
            Division indicator, default: False
        para_name: str
            The name of the parameter
        colour: str
            The colour of the lines.
        name: list
            A list of galaxy names.
        
        Return
        ------
        plot
        
        """
        
        total_mag1 = SRead.grab_total_mag(list_input1)
        total_mag2 = SRead.grab_total_mag(list_input2)
        
        plot = ShowcaseCompare2.plot_compare_generic(total_mag1, total_mag2, 
                                                     sub=sub , div=div, 
                                                     para_name=para_name, 
                                                     colour=colour,name=name,
                                                     label=label)
        return plot

#%%
#Mass_Re_plot_hist

#hist_4plot
#hist_4plot_norm
#hist_4plot_cum_KS
#selection_distant_Mass
    
    
#%% Construction area
## warning, no judgement, fuck off
from astropy.io import fits
import cmasher as cmr

class PlotHist(ShowcaseCompare2, ShowcaseIndi):
    
    #def __init__(file):
    #    self.file = file
        
    def A(mass, d, colour, legend):
        vdis_mass_plot(mass, d, colour, legend)

#%%
class Plot2D(object):
    """
    
    """
    def __init__(self, file):
        self.file = file    
    
    
    def read_fits_img(file):
        """
        Read 2D image from Fits files
        
        """
        hdul = fits.open(file)
        #hdul.info()
        data = np.array(hdul[0].data)
        return data
    
    def edge_evaluation():
        return None
    
    def trunk_window(array,centre,r_max):
        """
        Resize the window of plots to focus on one specfic galaxy.

        Parameters
        ----------
        array : TYPE
            DESCRIPTION.
        centre : TYPE
            DESCRIPTION.
        r_max : TYPE
            DESCRIPTION.

        Returns
        -------
        window : TYPE
            DESCRIPTION.

        """

        x0,y0 = centre[0],centre[1]        
        #print('edges', x0-r_max,x0+r_max,y0-r_max,y0+r_max)       
        window = array[y0-r_max:y0+r_max,x0-r_max:x0+r_max]
        return window
    
    
    def flatten2D(data,condition = 0, value = 1):
        """
        Normalize the image base on some condition
        
    
        Parameters
        ----------
        data : 2D list, 2D array
            The input image.
        
        condition : float, optional
            The condition to be evaluated. 
            The default is 0.
            
        value : float, optional
            The value to be replace if the condition is not met.
            The default is 1.

        Returns
        -------
        data : 2D list, 2D array
            The new image.

        """
        for i in range(len(data)):
            for j in range(len(data[i])):
                if data[i][j] ==condition:
                    pass
                else:
                    data[i][j] = value
        return data
                    
    # tested
    def x_average(matrix):
        """
        

        Parameters
        ----------
        matrix : TYPE
            DESCRIPTION.

        Returns
        -------
        array : TYPE
            DESCRIPTION.

        """
        
        array = np.array([np.average(matrix[:,x]) 
                          for x in range(len(matrix[:,0]))])

        return array
    # tested
    def y_average(matrix):
        """
        

        Parameters
        ----------
        matrix : TYPE
            DESCRIPTION.

        Returns
        -------
        array : TYPE
            DESCRIPTION.

        """
        
        array = np.array([np.average(matrix[x,:]) 
                          for x in range(len(matrix[0,:]))])

        return array
    
    #
    def find_mode_sky(file_name,showplot=False, saveplot= True):
        """
        

        Parameters
        ----------
        file_name : TYPE
            DESCRIPTION.
        showplot : TYPE, optional
            DESCRIPTION. The default is False.
        saveplot : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        None.

        """
        #N=np.loadtxt(file_name)
        N = Plot2D.read_fits_img(file_name)
        q=10000
        N = N.flatten() # turn 2D image into 1D array
                
        fig , ax = plt.subplots()
        n, bins, patches = plt.hist(N, q, facecolor='g', alpha=1.0)
        w=0.5 #width
        peaks, _ = find_peaks(n,height=max(n))
        
        #check for peaks array size, pick the middle
        if len(peaks) == 1:
            pass
        #elif (len(peaks) % 2 ) == 0:
        #    peaks = peaks[(len(peaks)/2)+1]
        else:
            peaks = int(np.median(peaks))
        
        xmin, xmax, ymin , ymax = abs(np.median(n))-w, abs(np.median(n))+w, \
            0, max(n)+max(n)/8        
        
        plt.plot((bins[peaks]+bins[peaks+1])/2, (n[peaks]), "x")
        plt.xlabel(r'$\rm Pixel Value$')
        plt.ylabel(r'$\rm Count$')
        #plt.text((bins[peaks]+bins[peaks+1])/2+(bins[peaks]+bins[peaks+1])*3, 
        #         (n[peaks]), r'%s' %(file_name))
        plt.text((bins[peaks]+bins[peaks+1])/2+(bins[peaks]+bins[peaks+1])*3, 
                 (n[peaks]-1), r'peak value= {:.3f}'.format(float(bins[peaks])))
        plt.title('Histogram of {}'.format(file_name))
        plt.axis([xmin, xmax, ymin , ymax])
        plt.grid(True)
        #plt.show()

        plt.savefig('{}_sky.png'.format(file_name))
        #np.savefig('foo.pdf')
        K,KK=(n[peaks]+n[peaks+1])/2, np.argmax(n)
        #print(K, KK)
        print(file_name, n[peaks], bins[peaks])
        return None

    
    def plot_galaxy(data,val_min,val_max, centre):    
        """
        Plot a single galaxy image
        
        (x y inversion)

        ----------
        file : str
                                    
        val_min, val_max: float
            
            
        centre: tuple
           
        
        Return
        ------
        plot
        
        """
        ## mode to centre on the average value
        
        
        avg_data, s_data = np.average(data),np.std(data)
        a =15
        data_log = np.log(a*data+1) / np.log(a)
        #plt.imshow(data,cmap=cmr.heat)
        index = np.linspace(0,len(data[:,0]),len(data[:,0])/10)
        new_scale = (index - np.full(index.shape, centre[1]))*0.4
        new_scale = np.array([format(new_scale[x], '.2f') for x in range(
            len(new_scale))])

                
        plt.imshow(data_log,cmap=cmr.heat, vmin= val_min, vmax = val_max)


        #plt.yticks(index, new_scale)
        #plt.hlines(centre[1],0,centre[0],'green')
        #plt.vlines(centre[0],0,centre[1],'green')
        plt.xlabel(r"$\rm arcsec$")
        plt.ylabel(r"$\rm arcsec$")
        plt.show()
        
    def plot_galaxy_3plot(file_name,md_file_name, res_file_name,
                         centre,r_max=400, a=15):
        """
        Plot individual galaxy 
        ----------
        file_name : str
            The fits file of the original galaxy image 
            
        md_file_name: str
            The fits file of the galaxy model image 

        res_file_name: str
        
        centre_x: tuple
            
        Optional
        ---------
        r_max:
            
        a:
               
        
        Return
        ------
        plot
        
        """
        x0,y0 = centre[0],centre[1]

        data0 = Plot2D.read_fits_img(file_name)   
        data1 = Plot2D.read_fits_img(md_file_name)        
        data2 = Plot2D.read_fits_img(res_file_name)    
        
        index = np.linspace(0,len(data0[:,0]),len(data0[:,0]))


        data_trunk0 = Plot2D.trunk_window(data0,centre,r_max)
        data_log0 = np.log(a*data_trunk0+1) / np.log(a)

        data_trunk1 = Plot2D.trunk_window(data1,centre,r_max)
        data_log1 = np.log(a*data_trunk1+1) / np.log(a)
         
        data_trunk2 = Plot2D.trunk_window(data2,centre,r_max)
        data_log2 = np.log(a*data_trunk2+1) / np.log(a)

        #############diagnoisis mode
        #data_trunk1_mask = Plot2D.flatten2D(data_trunk1)
        #
        #data_trunk0 = data_trunk0*data_trunk1_mask
        #data_log0= np.log(a*data_trunk0+1) / np.log(a)
        #
        #data_trunk2 = data_trunk2*data_trunk1_mask
        #data_log2 = np.log(a*data_trunk2+1) / np.log(a)
        ################

        l= len(data_trunk0[:,0])
                
        a,b,c,d,e = 0, int((l-1)/2)-int((l-1)/4),int((l-1)/2),\
                        int((l-1)/2)+int((l-1)/4),l-1 
        C=np.array([a,b,c,d,e])
        
        B=(C-np.full(5,c))*0.4
        B = np.array([format(B[x], '.2f') for x in range(len(B))])    
        
        index_x0 = np.linspace(0,len(data_log0[:,0]),len(data_log0[:,0]))
        index_y0 = np.linspace(0,len(data_log0[0,:]),len(data_log0[0,:]))
        index_x1 = np.linspace(0,len(data_log1[:,0]),len(data_log1[:,0]))
        index_y1 = np.linspace(0,len(data_log1[0,:]),len(data_log1[0,:]))
        index_x2 = np.linspace(0,len(data_log2[:,0]),len(data_log2[:,0]))
        index_y2 = np.linspace(0,len(data_log2[0,:]),len(data_log2[0,:]))  
        
        avg_data0, std_data0 = np.average(data_trunk0), np.std(data_trunk0)
        avg_data1, std_data1 = np.average(data_trunk1), np.std(data_trunk1)
        avg_data2, std_data2 = np.average(data_trunk2), np.std(data_trunk2)
        
        val_min0, val_max0 =  avg_data0 , avg_data0+1.5*std_data0
        
            
        fig = plt.figure()

        gs = gridspec.GridSpec(ncols=3, nrows=2,height_ratios=[1,3], 
                               hspace=0, wspace=0.05)
    
        ax_main0 = fig.add_subplot(gs[3])      
        ax_xDist0 = fig.add_subplot(gs[0])

        ax_main1 = fig.add_subplot(gs[4])      
        ax_xDist1 = fig.add_subplot(gs[1])
        
        ax_main2 = fig.add_subplot(gs[5])      
        ax_xDist2 = fig.add_subplot(gs[2])


        ax_main0.imshow(data_log0, cmap=cmr.heat, vmin=val_min0, vmax=val_max0)
        #ax_main0.set_xlabel("arcsec",fontsize=16)
        ax_main0.set_ylabel(r"$\rm arcsec$",fontsize=24)
        ax_main0.set_xticks([])
        ax_main0.set_yticks(C)
        ax_main0.set_yticklabels(B,fontsize=12)
        
        ax_xDist0.plot(index_x0,Plot2D.x_average(data_trunk0))
        #ax_xDist0.hlines(avg_data0,0,len(index_x0),linestyle="dashed") 
        #ax_xDist0.hlines(0,0,len(index_x0),linestyle="dashed")        

        ax_xDist0.set_ylabel(r'$\rm count$',fontsize=24)
        ax_xDist0.set_xticks([])
        ax_xDist0.set_ylim(bottom=0, top=val_max0)    
        #ax_xDist0.set_yticklabels(,fontsize=18)

        #val_min0-0.3*std_data0

        ax_main1.imshow(data_log1, cmap=cmr.heat, vmin=val_min0, vmax=val_max0)
        ax_main1.set_yticks([])
        ax_main1.set_xticks([])

        ax_xDist1.plot(index_x1,Plot2D.x_average(data_trunk1))
        #ax_xDist1.hlines(avg_data1,0,len(index_x1),linestyle="dashed") 
        #ax_xDist1.hlines(0,0,len(index_x1),linestyle="dashed")        

        ax_xDist1.set_xticks([])
        ax_xDist1.set_yticks([])
        ax_xDist1.set_ylim(bottom=0, top=val_max0)       
#val_min0-0.3*std_data0
        ax_main2.imshow(data_log2, cmap=cmr.heat, vmin=val_min0, vmax=val_max0)
        ax_main2.set_yticks([])
        ax_main2.set_xticks([])

        ax_xDist2.plot(index_x2,Plot2D.x_average(data_trunk2))
        #ax_xDist2.hlines(avg_data2,0,len(index_x2),linestyle="dashed")   
        #ax_xDist2.hlines(0,0,len(index_x2),linestyle="dashed")        

        ax_xDist2.set_xticks([])
        ax_xDist2.set_yticks([])

        ax_xDist2.set_ylim(bottom=0, top=val_max0)       
#
#val_min0-0.3*std_data0
        
        plt.savefig("%s.pdf"%file_name, dpi=200)
        return fig
    
    def plot2_import_profiler(image1,image2):
        
        return None
    
    def res_analysis():
        return None




img ="/home/dexter/result/image_plot/fit_example/NGC2872.fits"
md = "/home/dexter/result/image_plot/fit_example/md1_NGC2872.fits"
res= "/home/dexter/result/image_plot/fit_example/res1_NGC2872.fits"

#Plot2D.plot_galaxy_3plot(img,md,
#                         res, (252, 432),r_max=207, a =15)
#plt.show() #2872


#Plot2D.plot_galaxy_3plot(img,md,
  #                       res, (1088, 789),r_max=255, a =15)
#plt.show()4045
#
#Plot2D.plot_galaxy_3plot(img,md,
#                         res, (271, 1080),r_max=250, a =15)
#plt.show() #3675

# need to fix the edge problem
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
