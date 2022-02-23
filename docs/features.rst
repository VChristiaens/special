What makes it `special`?
------------------------

This package provides the following tools for the analysis of measured spectra:

* calculation of the spectral correlation between channels of an IFS datacube;
* calculation of empirical spectral indices for MLT-dwarfs;
* fitting of input spectra to different grids of models, including additional parameters such as (extra) black body component(s) and extinction;
* using either the MCMC (`emcee <https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract>`_) or nested (`nestle <http://github.com/kbarbary/nestle>`_) sampler to infer posterior distributions on spectral model parameters;
* searching for the best-fit template spectrum within a given template library, with up to two free parameters (flux scaling and relative extinction).


The MCMC and nested sampler routines have been adapted to:

* be flexible, as they are usable on any grid of models provided by the user (along with a snippet function specifying the format of the input);
* sample the effect of (additional) blackbody components;
* sample the effect of extinction; 
* sample different extinction laws than ISM (parametrised using RV);
* accept either uniform or Gaussian priors for each model parameter;
* accept a prior on the mass of the object (if surface gravity is one of the model parameters, and for the MCMC sampler only);
* consider convolution with the spectral PSF, photometric filters transmission and/or resampling of the model for consistency with the input spectrum.
* use a log-likelihood expression that can include i) spectral correlation between measurements of adjacent channels of a given instrument, and ii) additional weights that are proportional to the relative spectral bandwidth of each measurement, in case these are obtained from different instruments (e.g. photometry+spectroscopy).

More details are available in `Christiaens et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021MNRAS.502.6117C/abstract>`_ (note it was originally implemented as ``specfit``, a former module of the ``VIP`` package).