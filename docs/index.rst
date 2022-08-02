.. This file should at least contain the root `toctree` directive.

.. image:: _static/Special_logo.png
   :align: center
   :width: 825px

Welcome to the `special` documentation
======================================

Recent technological progress in high-contrast imaging has allowed the spectral 
characterization of directly imaged giant planet and brown dwarf companions at 
ever shorter angular separation from their host stars, hence opening a new 
avenue to study their formation, evolution, and composition.
 
In this context, ``special`` is a Python package for the SPEctral Characterization of directly ImAged Low-mass companions. While some tools are specific to the characterisation of low-mass (M, L, T) dwarfs down to giant planets at optical/IR wavelengths, the main routines of ``special`` (MCMC and nested samplers) can also be used in a more general way for the characterisation of any type of object with a measured spectrum, provided a relevant input model grid, regardless of the observational method used to obtain the spectrum (direct imaging or not) and regardless of the format of the spectra (multi-band photometry, low-resolution or medium-resolution spectrum, or a combination thereof).

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   features
   trimmed_readme

.. toctree::
   :maxdepth: 1
   :caption: Tutorial
   :hidden:

   tutorials/walkthrough.ipynb

.. toctree::
   :maxdepth: 3
   :caption: About

   about

.. toctree::
   :maxdepth: 3
   :caption: Package content

   special



API
---
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

