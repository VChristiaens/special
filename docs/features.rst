What makes it `special`?
------------------------

The ``special`` package provides a number of tools for the analysis of spectra from any (sub)stellar object, regardless of the observational method used to obtain the spectra (direct imaging or not) and the format of the spectra (multi-band photometry, low-resolution or medium-resolution spectrum, or a combination thereof). That being said, the main routines in the package (e.g. Bayesian retrieval of model parameters though MCMC or nested samplers, or best-fit template search) can also be applied to the spectrum of any object, provided a relevant grid of models or library of templates for the fit.

The main available features of the package are listed below:

* calculation of the spectral correlation between channels of an IFS datacube (relevant to directly imaged companions with an IFS, where the uncertainty reflects spectrally correlated residual speckle noise);
* calculation of empirical spectral indices for MLT-dwarfs;
* fitting of input spectra to either photo-/atmospheric model grids or a blackbody model, including additional parameters such as (extra) black body component(s), extinction and total-to-selective extinction ratio;
* using either MCMC (`emcee <https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract>`_) or nested (`nestle <http://github.com/kbarbary/nestle>`_ or `UltraNest <https://johannesbuchner.github.io/UltraNest/>`_) samplers to infer posterior distributions on spectral model parameters in a Bayesian framework;
* searching for the best-fit template spectrum within a given template library, with up to two free parameters (relative flux scaling and extinction).


The MCMC and nested sampler routines have been adapted to:

* be flexible, as they are usable on any grid of models provided by the user (along with a snippet function specifying how to read the format of the input files);
* sample the effect of (additional) blackbody components;
* sample the effect of extinction (AV); 
* sample different extinction laws than ISM (parametrised using the total-to-selective extinction ratio RV);
* sample a list of potential emission lines;
* accept either uniform or Gaussian priors for each model parameter;
* accept a prior on the mass of the object (if surface gravity is one of the model parameters, and for the MCMC sampler only);
* consider convolution with the line spread function, photometric filters transmission and/or resampling of the model for consistency with the input spectrum - in particular convolution and resampling are done in two consecutive steps, and multiple resolving powers can be provided as input;
* use a log-likelihood expression that can include i) spectral correlation between measurements of adjacent channels of a given instrument, and ii) additional weights that are proportional to the relative spectral bandwidth of each measurement, in case these are obtained from different instruments (e.g. photometry+spectroscopy).

``special`` was originally developed for the specific need of characterizing (sub)stellar companion CrA-9b/B at a time when no other public package achieving this task was available (`Christiaens et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021MNRAS.502.6117C/abstract>`_ ). First implemented as ``specfit``, a former module of the ``VIP`` package, it then underwent significant expansion, with the current package now detailed in the following JOSS publication `Christiaens et al. (2022) <https://joss.theoj.org/papers/10.21105/joss.04456>`_ (preferred citation if you use this package).