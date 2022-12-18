
*GalSpheroirds*, an utility package that analyse galaxy photometric structures.
==============================================================================


Introduction
============
*GalSpheroids* is a pacakage desgined to analyse galaxy structure from multi-componenets sttructural decomposition.
It contains the utility functions to manipulate


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
