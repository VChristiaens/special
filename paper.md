---
title: 'special: A Python package for the spectral characterization of directly imaged low-mass companions'
tags:
  - Python
  - astronomy
  - exoplanets
  - high-contrast
  - spectroscopy
  - direct imaging
authors:
  - name: Valentin Christiaens
    orcid: 0000-0002-0101-8814
    affiliation: 1


affiliations:
  - name: Space sciences, Technologies & Astrophysics Research Institute, Université de Liège, Belgium
    index: 1


date: 19 May 2022
bibliography: paper.bib
---

# Summary

Recent technological progress in high-contrast imaging has allowed the spectral 
characterization of directly imaged giant planet and brown dwarf companions at 
ever shorter angular separation from their host stars, hence opening a new 
avenue to study their formation, evolution, and composition. In this context, 
``special`` is a Python package that was developed to provide the tools to 
analyse the low- to medium-resolution optical/IR spectra of these directly 
imaged low-mass companions.

# Statement of need

``special`` provides the following tools for the analysis of measured spectra:

* calculation of the spectral correlation between channels of an integral field
spectrograph (IFS) datacube [@Greco:2016; @Delorme:2017];

* calculation of empirical spectral indices for MLT-dwarfs 
[@Gorlova:2003; @Slesnick:2004; @Allers:2007], enabling their 
classification;

* fitting of input spectra to different (user-provided) grids of models, with 
the possibility to include additional parameters such as extra blackbody 
component(s) and extinction;

* estimating most likely model parameters in a Bayesian framework, using either 
MCMC [@Goodman:2010] or nested 
[@Skilling:2004; @Mukherjee:2006; @Feroz:2009; @Buchner:2021a] samplers to infer 
their posterior distributions;

* searching for the best-fit template spectrum within a given template library, 
with up to two free parameters (flux scaling and relative extinction).

The MCMC sampler relies on ``emcee`` 
[@Foreman-Mackey:2013; @Foreman-Mackey:2019], while two options are 
available for the nested sampling: ``nestle`` [@nestle] and ``ultranest`` 
[@Buchner:2021b]. The samplers have been adapted for flexibility - 
they are usable on any grid of models provided by the user, simply requiring a 
snippet function specifying the format of the input. Moreover they can sample 
the effect of blackbody component(s) (either as a separate model or as extra 
components to an atmospheric model), extinction, and different extinction laws 
than ISM. The samplers can accept either uniform or Gaussian priors for each 
model parameter. In the case of the MCMC sampler, a prior on the mass of the 
object can also be provided (if surface gravity is one of the model parameters). 
The code also considers convolution and resampling of model spectra to match 
the observed spectrum. Either spectral resolution or photometric filter 
transmission (or combinations thereof for compound input spectra) can be 
provided as input to the algorithm, for appropriate convolution/resampling of
different parts of the model spectrum. The adopted log-likelihood expression 
can include i) spectral correlation between measurements of adjacent channels 
of a given instrument, and ii) additional weights that are proportional to the 
relative spectral bandwidth of each measurement, in case these are obtained 
from different instruments (e.g. photometry+spectroscopy).

Finally, a jupyter notebook tutorial is made available on the, which covers a 
significant fraction of the available features in ``special`` applied to the 
case of the composite spectrum of CrA-9 B/b [@Christiaens:2021].

# Acknowledgements

VC acknowledges financial support from the Belgian F.R.S.-FNRS.

# References
