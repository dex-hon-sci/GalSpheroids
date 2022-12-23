 |arXiv|

*GalSpheroirds*, an utility package that analyse galaxy structures.
===================================================================

Introduction
============
*GalSpheroids* is a pacakage desgined to analyse galaxy structure from multi-componenets structural decomposition.

Galaxies are complex systems that are made of trillion of stars. 
While it is difficult to model such a system by tracking the dynamics of individaul stars, it is possible to use a series of mathematical function to depict the shape and feature of its surface brightness profile. 
This practise is known as multi-components structural decomposition. 

These are a few common components astronomer describe a galaxy with,

Galaxy Bundle
-------------
*GalSpheroid* operate on what I called a *Galaxy Bundle*, a list-based data format that contains all the model parameters of galaxy structural components.


data format
Galaxy Bundle are formatted as follows,

[[NGC0001, residual_rms,"Component_1", [parameter_1,parameter_2,..], component_1_magnitude, "Component_2",...,"Total_mag", [Total_magnitude]], 
 [NGC0002, residual_rms,"Component_1", [parameter_1,..], component_1_magnitude,...,"Total_mag", [Total_magnitude], ...]
 
Each row contains the information of one galaxy. 


Modules
-------
This package contains four main modules, each with a specific functionality:

 -SphRead.py: It contains the functions used in reading the galaxy bundle,
 
 -SphSort.py: It contains the functions used in sorting and manipulate the data in the galaxy bundle,
 
 -SphPlot.py: It contains the functions used in plotting the data in the galaxy bundle,
 
 -SphAnalysis.py: It contains the functions used in analyzing, fitting, and modelling of the data in the galaxy bundle.
    

Utility example

Getting started
===============
Installation
------------
*GalSpheroids* can be easily installed by either cloning the `repository`_ and installing it manually::

    $ git clone https://github.com/dex-hon-sci/GalSpheroids
    $ cd GalSpheroids
    $ pip install .
    

Example usage
=============
To create a Galaxy Bundle
-------------------------
To add, read, replace, or remove data from a galaxy bundle.
-----------------------------------------------------------
.. code:: python

    # Import 
    import SphRead as SRead
    import SphSort as SSort
    
    # Read a galaxy bundle
    
    # L
  
To make a plot with the data
----------------------------
.. code:: python

    # Import 
    import SphPlot as Splot
    
    # input files
    image_name = "./NGC3805.fits" # Raw telescope image file (in fits format)
    md_name = "./md1_NGC3805.fits" # Isophotal model file (in fits format)
    res_name= "./res1_NGC3805.fits" # The residual image file (in fits format) by sutracting image with model
    
    centre = (1828.328, 996.8774) # The pixel indices for the centre of the galaxy 
    
    # Plot galaxy image, model, and residual
    SPlot.Plot2D.plot_galaxy_3plot(image_name,md_name, res_name, 
                               centre,r_max=400, alp=15)  

>>> 
.. image:: 

To make a Galaxy Bundle
-----------------------
.. code:: python

    # Imports
    import SphRead as SRead
    import SphSort as SSort
    
    # Define an override list
    override_list_equ = ["NGC2862",2,"Disk","NGC2872",5,"Point Source", "NGC3805",5,"IntDisk","NGC3805",8,"Point Source","NGC3812",2,"Background",\
 "NGC3872",5,"Halo","NGC3940",5,"Point Source", "NGC4065",5,"Point Source", "NGC4555",5,"Point Source",\
 "UGC8736",2,"Disk", "NGC5382", 5,"IntDisk","NGC4772", 14 ,"Point Source", "NGC4845", 14, "Point Source",
 "NGC5375",11 , "Point Source","NGC2832",5,"cD Halo"]
 
    # Read an ASCII file, "F_equvi_list_Bin1V.txt", which contain the name and location of "Profiler" output log file on each row for each galaxy.
    # Record and transform all the raw information into galaxy bundle format, named ""F_Gal_bundle_equvi_Bin1V".
    # "F_Gal_bundle_equvi_Bin1V" contains the label of each mathematical functions used to model galaxy components.     
    SRead.run_list("F_equvi_list_Bin1V.txt","F_Gal_bundle_equvi_Bin1V",True)
    
    # Separate components by the analytical function types for further assessement.
    C2 = SSort.cpt_seperator_demo('F_Gal_bundle_equvi_Bin1V')
    
    # Run diagnosis on the analytical functions and assign proper component name for each galaxies
    # Output a new galaxy bundle, 'F_Gal_bundle_equvi_Bin1V_cpt', with each component named.
    # Read an override list, "override_list_equ", to manually assign new name for special components.
    SSort.cpt_classifier_demo('F_Gal_bundle_equvi_Bin1V',C2,'F_Gal_bundle_equvi_Bin1V_cpt',override_list_equ)
    
    
>>> [['IC00983',0.0576173775624,'Bulge', array([19.20348468,  3.57098939,  1.70405518]), 13.4886302019,
  'Disk', array([21.6, 50.34337508]), 11.0948540004,'Ring',  array([23.40194399, 38.58954655,  8.01304625]), 15.1130784521,
  'Ring',array([24.28387892, 25.44196011,  6.30355541]),  16.7078413371,
  'Ring', array([23.50085385, 55.58167158, 12.02601045]), 14.3750236666,
  'Ring',  array([24.34763182, 67.56223984,  5.24776207]),15.9102391492,
  'PrimBar',array([20.77187614, 13.02105021,  2.00679133,  0.21960181]),14.6499487471,'Total_mag',[10.861464728]],
   ...
 ['NGC2796', 0.06062014856064842, 'CoreBulge', array([13.44195809,  4.01365186,  2.78400144,  0.66268723,  3.84146707, 0.12541042]),
  12.958529155512657,'Disk',  array([20.81851904, 15.93231547]), 12.811674436254169,'Total_mag', [12.128826141806979]]]
  
Community guidelines
====================
PRISM is an open-source and free-to-use software package (and it always will be), provided under the BSD-3 license.

Users are highly encouraged to make contributions to the package or request new features by opening a GitHub issue. 
If you would like to contribute to the package, but do not know what, then there are quite a few ToDos in the code that may give you some inspiration. 
As with contributions, if you find a problem or issue with PRISM, please do not hesitate to open a GitHub issue about it or post it on Gitter.

To acknowledge this work and reference the original galaxy structure data, please cite the following:

::

    @ARTICLE{2022MNRAS.514.3410H,
        author = {{Hon}, Dexter S. -H. and {Graham}, Alister W. and {Davis}, Benjamin L. and {Marconi}, Alessandro},
        title = "{Disc cloaking: Establishing a lower limit to the number density of local compact massive spheroids/bulges and the potential fate of some high-z red nuggets}",
        journal = {\mnras},
        keywords = {galaxies: abundances, galaxies: bulges, galaxies: discs, galaxies: elliptical and lenticular, cD, galaxies: evolution, galaxies: structure, Astrophysics - Astrophysics of Galaxies},
        year = 2022,
        month = aug,
        volume = {514},
        number = {3},
        pages = {3410-3451},
        doi = {10.1093/mnras/stac1171},
        archivePrefix = {arXiv},
        eprint = {2204.13408},
        primaryClass = {astro-ph.GA},
        adsurl = {https://ui.adsabs.harvard.edu/abs/2022MNRAS.514.3410H},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
      }
.. _paper: https://arxiv.org/abs/2204.13408


Acknowledgements
================
The project is made possible by using the following software packages: AstroPy (Astropy Collaboration et al. 2013, 2018), Cmasher
(van der Velden 2020), IRAF (Tody 1986, 1993), ISOFIT (Ciambur2015), Linmix (Kelly 2007),  Matplotlib (Hunter 2007), pandas (Mckinney et al. 2010), pickle (Van Rossum, G. 2020), NumPy (Harris et al. 2020), Profiler (Ciambur 2016), SAOImageDS9 (Joye & Mandel 2003), and SciPy (Virtanen et al. 2020)

.. |arXiv| image:: https://img.shields.io/badge/arXiv-1901.08725-brightgreen
    :target: https://arxiv.org/abs/2204.13408
    :alt: arXiv - Paper
 
