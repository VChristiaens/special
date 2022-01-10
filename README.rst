specfit
=======

.. image:: https://badge.fury.io/py/specfit.svg
    :target: https://pypi.python.org/pypi/specfit

.. image:: https://img.shields.io/badge/Python-3.5%2C%203.6%2C%203.7%2C%203.8%2C%203.9-brightgreen.svg
    :target: https://pypi.python.org/pypi/specfit

.. image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat
    :target: https://github.com/vortex-exoplanet/VIP/blob/master/LICENSE


`specfit` is a package for the spectral characterization of (sub-)stellar objects.

This package provides the tools for the analysis of measured spectra (e.g. low-res spectra of directly imaged companions, but not exclusively). It enables:
- fitting of input spectra to different grids of models (the code is grid-independent, download of the grid is left to the user);  
- MCMC sampling of the model grid parameter space;
- searching for the best-fit template spectrum within a given library.

In particular, the MCMC sampler routine makes use of the `emcee` package (Foreman-Mackey et al. 2013):
- is flexible, as it is usable on any grid of models downloaded by the user (provided a snippet function specifying the format of the input);
- can fit (additional) blackbody components;
- can sample the effect of extinction; 
- can sample different extinction laws than ISM (parametrised using RV);
- accepts either uniform or Gaussian priors for each model parameter; 
- accepts a prior on the mass of the object (if surface gravity is one of the model parameters);
- considers convolution with the spectral PSF and/or resampling of the model for consistency with the input spectrum;

Furthermore, the log-likelihood expression has the option to include:
- spectral correlation between measurements of adjacent channels of a given instrument;
- weights that are proportional to the relative spectral span of each measurement, in case these are obtained from different instruments.

More details are available in `Christiaens et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021MNRAS.502.6117C/abstract>`_.
Please cite this publication if you use `specfit` for your research, along with `Foreman-Mackey et al. (2013) <https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract>`_ if you use the MCMC sampler.