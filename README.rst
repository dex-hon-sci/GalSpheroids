
*GalSpheroirds*, an utility package that analyse galaxy photometric structures.
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
    *SphRead.py: It contains the functions used in reading the galaxy bundle,
    *SphSort.py: It contains the functions used in sorting and manipulate the data in the galaxy bundle,
    *SphPlot.py: It contains the functions used in plotting the data in the galaxy bundle,
    *SphAnalysis.py: It contains the functions used in analyzing, fitting, and modelling of the data in the galaxy bundle.
    

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
-------------
.. code:: python

    # Imports
    import SphRead as SRead
    import SphSort as SSort
    
    import SphPlot as SPlot
    
    import SphAnalysis as SAna
  
To make a plot with the data
-------------
.. code:: python

    # Import 
    import SphAnalysis as SAna

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
