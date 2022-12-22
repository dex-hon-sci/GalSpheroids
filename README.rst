
*GalSpheroirds*, an utility package that analyse galaxy structures.
==============================================================================


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

[[NGC0001, residual_rms,"Component_1", [parameter_1,parameter_2,..], component_1_magnitude, "Component_2",...,Total_magnitude], 
 [NGC0002, residual_rms,"Component_1", [parameter_1,..], component_1_magnitude,...,Total_magnitude], ...]
 
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
    
    
  
To add, read, replace, or remove data from a galaxy bundle.
-----------------------------------------------------------
.. code:: python

    # Import 
    import SphRead as SRead
    import SphSort as SSort
    
    # Read a galaxy bundle
  
To make a plot with the data
----------------------------
.. code:: python

    # Import 
    import SphPlot as Splot

Community guidelines
====================
To cite the use of the original gaalxy structure data, please cite the following

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
         
.._paper: https://arxiv.org/abs/2204.13408


Acknowledgements
================
