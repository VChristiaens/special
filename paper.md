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

``special`` provides a number of tools for the analysis of spectra from any (sub)stellar 
object, regardless of the observational method used to obtain the spectra (direct imaging 
or not) and the format of the spectra (multi-band photometry, low-resolution or 
medium-resolution spectrum, or a combination thereof). Although implemented with
the characterization of directly imaged substellar companions in mind, the main routines 
in ``special`` (e.g. Bayesian retrieval of model parameters though MCMC or nested 
samplers, or best-fit template search) can also be applied to the spectrum of any type of 
object, provided a relevant grid of models or library of templates for the fit.

``special`` shares similar basic utilities as offered in ``splat`` [@Burgasser:2017], such as 
dereddening, spectral indices calculation, model grid fitting through MCMC and template 
fitting. However, a number of features are currently unique to ``special``, such as (i) Bayesian 
inference through nested samplers; (ii) inclusion of non-grid parameters for model fits (e.g. 
extinction, extra blackbody components, specific emission lines); (iii) inclusion of relative 
extinction and flux scaling, and handling of spectral coverage mismatches when searching 
for the best-fit template in a library; (iv) empirical estimation of spectral correlation between 
channels of an integral field spectrograph, which is relevant to the directly imaged companions 
for which uncertainties in the spectrum capture correlated residual speckle noise [@Greco:2016]; 
and (v) compatibility of all ``special`` fitting routines with combined spectra (i.e. obtained with 
multiple instruments with potentially different resolving powers or photometric filters).

The main available features of the package are listed below:

* calculation of the spectral correlation between channels of an integral field
spectrograph (IFS) datacube [@Greco:2016; @Delorme:2017];

* calculation of empirical spectral indices for MLT-dwarfs 
[@Gorlova:2003; @Slesnick:2004; @Allers:2007], enabling their 
classification;

* fitting of input spectra to either photo-/atmospheric model grids or a blackbody model, 
including additional parameters such as (extra) black body component(s), extinction,
 total-to-selective extinction ratio or specific emission lines.

* estimating most likely model parameters in a Bayesian framework, using either 
MCMC [@Goodman:2010] or nested 
[@Skilling:2004; @Mukherjee:2006; @Feroz:2009; @Buchner:2021a] samplers to infer 
their posterior distributions;

* searching for the best-fit template spectrum within a given template library, 
with up to two free parameters (flux scaling and relative extinction).

The MCMC sampler relies on ``emcee`` 
[@Foreman-Mackey:2013; @Foreman-Mackey:2019], while two options are available 
for nested sampling: ``nestle`` [@nestle] and ``ultranest`` [@Buchner:2021b]. 
The samplers have been adapted for flexibility - they are usable on any grid of 
input models provided by the user, simply requiring a snippet function 
specifying the format of the input. Moreover they can sample the effect of 
blackbody component(s) (either as a separate model or as extra components to an 
atmospheric model), extinction, and different extinction laws than ISM. The 
samplers can accept either uniform or Gaussian priors for each model parameter. 
In the case of the MCMC sampler, a prior on the mass of the object can also be 
provided if surface gravity is one of the model parameters. The code also 
considers convolution and resampling of model spectra to match the observed 
spectrum. Either spectral resolution or photometric filter transmission (or 
combinations thereof for compound input spectra) can be provided as input to 
the algorithm, for appropriate convolution/resampling of different parts of the 
model spectrum. The adopted log-likelihood expression can include i) spectral 
covariance between measurements of adjacent channels of a given instrument, 
and ii) additional weights that are proportional to the relative spectral 
bandwidth of each measurement, in case these are obtained from different 
instruments (e.g. photometry+spectroscopy):

\begin{equation}
\label{Eq:logL}
\log \mathcal{L}(D|M) = - \frac{1}{2} \big[\mathbf{W}\odot(\mathbf{F_{\rm obs}}-\mathbf{F_{\rm mod}})\big]^T \mathbf{C^{-1}} \big[\mathbf{W}\odot(\mathbf{F_{\rm obs}}-\mathbf{F_{\rm mod}})\big]
\end{equation}

where $D$ are the data at hand (measured fluxes and spectral covariance), 
$M$ is the considered model, $\mathbf{F_{\rm obs}}$ and $\mathbf{F_{\rm mod}}$ 
are the fluxes of the observed and model spectra respectively (both are vectors of 
length $n_z$, the number of spectro-/photometric points), $\mathbf{C}$ is 
the spectral covariance matrix ($n_z$ x $n_z$), $\odot$ stands for the Hadamard product, and 
$\mathbf{W}$ is the vector of weights $w_i \propto \Delta\lambda_i/\lambda_i$ (length $n_z$), 
with $\Delta\lambda_i$ the width of spectral channels (for integral field spectrograph 
points) or the FWHM of photometric filters.

A Jupyter notebook tutorial illustrates most available features in ``special`` 
through their application for the analysis of the composite spectrum of CrA-9 B/b 
[@Christiaens:2021]. It is available on 
[GitHub](https://github.com/VChristiaens/special_extras), 
[Binder](https://mybinder.org/v2/gh/VChristiaens/special_extras/main) and the 
[documentation](https://special.readthedocs.io/en/latest/) of ``special``.

# Acknowledgements

VC acknowledges financial support from the Belgian F.R.S.-FNRS.

# References
