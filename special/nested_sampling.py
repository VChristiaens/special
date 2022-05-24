#! /usr/bin/env python

"""
Module with functions for posterior sampling of the model spectra parameters 
using nested sampling (either ``nestle`` or ``ultranest`` samplers).

.. [nes13]
   | Barbary 2013
   | **nestle**
   | *GitHub*
   | `https://github.com/kbarbary/nestle
     <https://github.com/kbarbary/nestle>`_

.. [SKI04]
   | Skilling 2004
   | **Bayesian Inference and Maximum Entropy Methods in Science and Engineering: 
     24th International Workshop on Bayesian Inference and Maximum Entropy 
     Methods in Science and Engineering**
   | *American Institute of Physics Conference Series, Volume 735, pp. 395-405*
   | `https://ui.adsabs.harvard.edu/abs/2004AIPC..735..395S
     <https://ui.adsabs.harvard.edu/abs/2004AIPC..735..395S>`_

.. [MUK06]
   | Mukherjee et al. 2006
   | **A Nested Sampling Algorithm for Cosmological Model Selection**
   | *ApJL, Volume 638, Issue 2, pp. 51-54*
   | `https://arxiv.org/abs/astro-ph/0508461
     <https://arxiv.org/abs/astro-ph/0508461>`_

.. [FER09]
   | Feroz et al. 2009
   | **MULTINEST: an efficient and robust Bayesian inference tool for cosmology 
     and particle physics**
   | *MNRAS, Volume 398, Issue 4, pp. 1601-1614*
   | `https://arxiv.org/abs/astro-ph/0809.3437
     <https://arxiv.org/abs/astro-ph/0809.3437>`_
     
.. [BUC21]
   | Buchner 2021
   | **UltraNest - a robust, general purpose Bayesian inference engine**
   | *JOSS, Volume 6, Issue 60, p. 3001*
   | `https://arxiv.org/abs/astro-ph/2101.09604
     <https://arxiv.org/abs/astro-ph/2101.09604>`_
"""



__author__ = 'V. Christiaens',
__all__ = ['nested_spec_sampling',
           'show_nestle_results',
           'show_ultranest_results']

from astropy import constants as con
import nestle
import ultranest
from ultranest.plot import cornerplot, PredictionBand
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from .config import time_ini, timing
from .fits import open_fits, write_fits
from .mcmc_sampling import lnlike, confidence, show_walk_plot, show_corner_plot
from .model_resampling import make_resampled_models, make_model_from_params
from .utils_nested import un_burning
from os.path import isfile
from scipy.special import ndtri


def nested_spec_sampling(init, lbda_obs, spec_obs, err_obs, dist, 
                         grid_param_list, labels, bounds, resamp_before=True, 
                         model_grid=None, model_reader=None, em_lines={}, 
                         em_grid={}, dlbda_obs=None, instru_corr=None, 
                         instru_res=None, instru_idx=None, use_weights=True,
                         filter_reader=None, AV_bef_bb=False, units_obs='si', 
                         units_mod='si', interp_order=1, priors=None, 
                         physical=True, interp_nonexist=True, 
                         output_dir='special/', grid_name='resamp_grid.fits', 
                         sampler='ultranest', method='single', npoints=100, 
                         dlogz=0.2, verbose=True, **kwargs):
                            
    
    """ Runs a nested sampling algorithm in order to retrieve the best-fit 
    parameters for given spectral model and observed spectrum.. The
    result of this procedure is either a ``nestle`` [nes13] or ``ultranest`` 
    [BUC21] object (depending on the chosen sampler) containing the samples 
    from the posterior distributions of each of the parameters.
    For the ``nestle`` sampler, several methods are available corresponding to
    MCMC [SKI04], single ellipsoid [MUK06] or multiple ellipsoid [FER09].

    Parameters
    ----------
    init: numpy ndarray or tuple
        Initial guess on the best fit parameters of the spectral fit. Length of 
        the tuple should match the total number of free parameters. 
        - first all parameters related to loaded models (e.g. 'Teff', 'logg')
        - then the planet photometric radius 'R', in Jupiter radius
        - (optionally) the intensity of emission lines (labels must match \
        those in the em_lines dict), in units of the model spectrum (x mu)
        - (optionally) the optical extinction 'Av', in mag
        - (optionally) the ratio of total to selective optical extinction 'Rv'
        - (optionally) 'Tbb1', 'Rbb1', 'Tbb2', 'Rbb2', etc. for each extra bb \
        contribution
    lbda_obs : numpy 1d ndarray or list
        Wavelength of observed spectrum. If several instruments, should be 
        ordered per instrument, not necessarily as monotonically increasing 
        wavelength. Hereafter, n_ch = len(lbda_obs).
    spec_obs : numpy 1d ndarray or list
        Observed spectrum for each value of lbda_obs.
    err_obs : numpy 1d/2d ndarray or list
        Uncertainties on the observed spectrum. If 2d array, should be [2,n_ch]
        where the first (resp. second) column corresponds to lower (upper) 
        uncertainty, and n_ch is the length of lbda_obs and spec_obs.
    dist :  float
        Distance in parsec, used for flux scaling of the models.
    grid_param_list : list of 1d numpy arrays/lists OR None
        - If list, should contain list/numpy 1d arrays with available grid of \
        model parameters. 
        - Set to None for a pure n-blackbody fit, n=1,2,...
        - Note1: model grids should not contain grids on radius and Av, but \
        these should still be passed in initial_state (Av optional).
        - Note2: for a combined grid model + black body, just provide \
        the grid parameter list here, and provide values for 'Tbbn' and 'Rbbn' \
        in initial_state, labels and bounds.
    labels: Tuple or list of strings
        Tuple of labels in the same order as initial_state, that is:
        - first all parameters related to loaded models (e.g. 'Teff', 'logg')
        - then the planet photometric radius 'R', in Jupiter radius
        - (optionally) the flux of emission lines (labels should match those \
        in the em_lines dictionary), in units of the model spectrum (times mu)
        - (optionally) the optical extinction 'Av', in mag
        - (optionally) the ratio of total to selective optical extinction 'Rv'
        - (optionally) 'Tbb1', 'Rbb1', 'Tbb2', 'Rbb2', etc. for each extra bb \
        contribution.      
    bounds: dictionary
        Each entry should be associated with a tuple corresponding to lower and 
        upper bounds respectively. Bounds should be provided for ALL model
        parameters, including 'R' (planet photometric radius). 'Av' (optical 
        extinction) is optional. If provided here, Av will also be fitted.
        Example for BT-SETTL: bounds = {'Teff':(1000,2000), 'logg':(3.0,4.5),
        'R':(0.1,5), 'Av':(0.,2.5)}
        'M' can be used for a prior on the mass of the planet. In that case the
        corresponding prior log probability is computed from the values for 
        parameters 'logg' and 'R' (if both exist).
    resamp_before: bool, optional
        Whether to prepare the whole grid of resampled models before entering 
        the MCMC, i.e. to avoid doing it at every MCMC step. Recommended.
        Only reason not to: model grid is too large and individual models 
        require being opened and resampled at each step.
    model_grid : numpy N-d array, optional
        If provided, should contain the grid of model spectra for each
        free parameter of the given grid. I.e. for a grid of n_T values of Teff 
        and n_g values of Logg, the numpy array should be n_T x n_g x n_ch x 2, 
        where n_ch is the number of wavelengths for the observed spectrum,
        and the last 2 dims are for wavelength and fluxes respectively.
        If provided, takes precedence over filename/file_reader.
    model_reader : python routine, optional
        External routine that reads a model file and returns a 2D numpy array, 
        where the first column corresponds to wavelengths, and the second 
        contains model values. See example routine in model_interpolation() 
        description.
    em_lines: dictionary, opt
        Dictionary of emission lines to be added on top of the model spectrum.
        Each dict entry should be the name of the line, assigned to a tuple of
        4 values: 
        1) the wavelength (in mu); 
        2) a string indicating whether line intensity is expressed in flux 
        ('F'), luminosity ('L') or log(L/LSun) ("LogL");
        3) the FWHM of the gaussian (or None if to be set automatically); 
        4) whether the FWHM is expressed in 'nm', 'mu' or 'km/s'. 
        The third and fourth can also be set to None. In that case, the FWHM of 
        the gaussian will automatically be set to the equivalent width of the
        line, calculated from the flux to be injected and the continuum 
        level (measured in the grid model to which the line is injected). 
        Examples: 
        em_lines = {'BrG':(2.1667,'F', None, None)};
        em_lines = {'BrG':(2.1667,'LogL', 100, 'km/s')}
    em_grid: dictionary pointing to lists, opt
        Dictionary where each entry corresponds to an emission line and points
        to a list of values to inject for emission line fluxes. For computation 
        efficiency, interpolation will be performed between the points of this 
        grid during the MCMC sampling. Dict entries should match labels and 
        em_lines.
    dlbda_obs: numpy 1d ndarray or list, optional
        Respective spectral channel width and FWHM of the photometric filters 
        used for the input spectrum. It is used to infer which part(s) of a 
        combined spectro+photometric spectrum should involve 
        convolution+subsampling (model resolution higher than measurements),
        interpolation (the opposite), or convolution by the transmission curve
        of a photometric filter. If not provided, it will be inferred from the
        difference between consecutive lbda_obs points (i.e. inaccurate for a 
        combined spectrum).
    instru_corr : numpy 2d ndarray or list, optional
        Spectral correlation throughout post-processed images in which the 
        spectrum is measured. It is specific to the combination of instrument, 
        algorithm and radial separation of the companion from the central star.
        Can be computed using distances.spectral_correlation(). In case of
        a spectrum obtained with different instruments, build it with
        distances.combine_corrs(). If not provided, it will consider the 
        uncertainties in each spectral channels are independent. See Greco & 
        Brandt (2017) for details.
    instru_res : float or list of floats/strings, optional
        The mean instrumental resolving power(s) OR filter names. This is 
        used to convolve the model spectrum. If several instruments are used, 
        provide a list of resolving power values / filter names, one for 
        each instrument used.
    instru_idx: numpy 1d array, optional
        1d array containing an index representing each instrument used 
        to obtain the spectrum, label them from 0 to n_instru. Zero for points 
        that don't correspond to any instru_res provided above, and i in 
        [1,n_instru] for points associated to instru_res[i-1]. This parameter 
        must be provided if the spectrum consists of points obtained with 
        different instruments.
    use_weights: bool, optional
        For the likelihood calculation, whether to weigh each point of the 
        spectrum based on the resolving power or bandwith of photometric
        filters used. Weights will be proportional to dlbda_obs/lbda_obs if 
        dlbda_obs is provided, or set to 1 for all points otherwise.
    filter_reader: python routine, optional
        External routine that reads a filter file and returns a 2D numpy array, 
        where the first column corresponds to wavelengths, and the second 
        contains transmission values. Important: if not provided, but strings 
        are detected in instru_res, the default file reader will be used. 
        It assumes the following format for the files:
        - first row containing header
        - starting from 2nd row: 1st column: wavelength, 2nd col.: transmission
        - Unit of wavelength can be provided in parentheses of first header \
        key name: e.g. "WL(AA)" for angstrom, "wavelength(mu)" for micrometer \
        or "lambda(nm)" for nanometer. Note: Only what is in parentheses \
        matters.
        Important: filter files should all have the same format and WL units.
    AV_bef_bb: bool, optional
        If both extinction and an extra bb component are free parameters, 
        whether to apply extinction before adding the BB component (e.g. 
        extinction mostly from circumplanetary dust) or after the BB component
        (e.g. mostly insterstellar extinction).
    units_obs : str, opt {'si','cgs','jy'}
        Units of observed spectrum. 'si' for W/m^2/mu; 'cgs' for ergs/s/cm^2/mu 
        or 'jy' for janskys.
    units_mod: str, opt {'si','cgs','jy'}
        Units of the model. 'si' for W/m^2/mu; 'cgs' for ergs/s/cm^2/mu or 'jy'
        for janskys. If different to units_obs, the spectrum units will be 
        converted.
    interp_order: int, opt, {-1,0,1} 
        Interpolation mode for model interpolation.
        -1: log interpolation (i.e. linear interpolatlion on log(Flux))
        0: nearest neighbour model.
        1: Order 1 spline interpolation.
    priors: dictionary, opt
        If not None, sets prior estimates for each parameter of the model. Each 
        entry should be set to either None (no prior) or a tuple of 2 elements 
        containing prior estimate and uncertainty on the estimate.
        Missing entries (i.e. provided in bounds dictionary but not here) will
        be associated no prior.
        e.g. priors = {'Teff':(1600,100), 'logg':(3.5,0.5), 'R':(1.6,0.1), 
        'Av':(1.8,0.2), 'M':(10,3)}
        Important: dictionary entry names should match exactly those of bounds.
    physical: bool, opt
        In case of extra black body component(s) to a photosphere, whether to 
        force lower temperature than the photosphere effective temperature.
    interp_nonexist: bool, opt
        Whether to interpolate non-existing models in the grid. Only used if 
        resamp_before is set to True.
    w : float or tuple
        The relative size of the bounds (around the initial state ``init``) for 
        each parameter. If a float the same relative size is considered for 
        each parameter. E.g. if 0.1, bounds will be set to:
        (0.9*params[0], 1.1*params[0]),
        ...
        (0.9*params[N-1], 1.1*params[N-1]),
        to True), or make it and write it if it does not.
    output_dir: str, optional
        The name of the output directory which contains the output files in the 
        case  ``save`` is True.   
    grid_name: str, optional
        Name of the fits file containing the model grid (numpy array) AFTER
        convolution+resampling as the observed spectrum given as input.
        If provided, will read it if it exists (and resamp_before is set
    method : {"single", "multi", "classic"}, str optional
        Flavor of nested sampling.
    npoints : int optional
        Number of active points. Must be >> ndim+1, otherwise will produce bad 
        results. For UltraNest, this is the minimum number of active points.
    dlogz : float, optional
        Target evidence uncertainty. Iterations will stop when the estimated 
        contribution of the remaining prior volume to the total evidence falls 
        below this threshold. Explicitly, the stopping criterion is 
        log(z + z_est) - log(z) < dlogz
        where z is the current evidence from all saved samples, and z_est is the
        estimated contribution from the remaining volume. This option and
        decline_factor are mutually exclusive. If neither is specified, the
        default is dlogz=0.2.
    decline_factor : float, optional
        If supplied, iteration will stop when the weight (likelihood times prior
        volume) of newly saved samples has been declining for
        decline_factor * nsamples consecutive samples. A value of 1.0 seems to
        work pretty well.
    rstate : random instance, optional
        RandomState instance. If not given, the global random state of the
        numpy.random module will be used.
    kwargs: optional
        Additional optional arguments to either the `nestle.sample` or 
        `ultranest.ReactiveNestedSampler.run` functions.

    Returns
    -------
    res : nestle object
        ``Nestle`` object with the nested sampling results, including the
        posterior samples.

    Notes
    -----
    Nested Sampling is a computational approach for integrating posterior
    probability in order to compare models in Bayesian statistics. It is similar
    to Markov Chain Monte Carlo (MCMC) in that it generates samples that can be
    used to estimate the posterior probability distribution. Unlike MCMC, the
    nature of the sampling also allows one to calculate the integral of the
    distribution. It also happens to be a pretty good method for robustly
    finding global maxima.

    Nestle documentation:
    http://kbarbary.github.io/nestle/

    Convergence:
    http://kbarbary.github.io/nestle/stopping.html
    Nested sampling has no well-defined stopping point. As iterations continue,
    the active points sample a smaller and smaller region of prior space.
    This can continue indefinitely. Unlike typical MCMC methods, we don't gain
    any additional precision on the results by letting the algorithm run longer;
    the precision is determined at the outset by the number of active points.
    So, we want to stop iterations as soon as we think the active points are
    doing a pretty good job sampling the remaining prior volume - once we've
    converged to the highest-likelihood regions such that the likelihood is
    relatively flat within the remaining prior volume.

    Method:
    The trick in nested sampling is to, at each step in the algorithm,
    efficiently choose a new point in parameter space drawn with uniform
    probability from the parameter space with likelihood greater than the
    current likelihood constraint. The different methods all use the
    current set of active points as an indicator of where the target
    parameter space lies, but differ in how they select new points from  it.
    "classic" is close to the method described in Skilling (2004).
    "single", Mukherjee, Parkinson & Liddle (2006), Determines a single
    ellipsoid that bounds all active points,
    enlarges the ellipsoid by a user-settable factor, and selects a new point
    at random from within the ellipsoid.
    "multiple", Shaw, Bridges & Hobson (2007) and Feroz, Hobson & Bridges 2009
    (Multinest). In cases where the posterior is multi-modal,
    the single-ellipsoid method can be extremely inefficient: In such
    situations, there are clusters of active points on separate
    high-likelihood regions separated by regions of lower likelihood.
    Bounding all points in a single ellipsoid means that the ellipsoid
    includes the lower-likelihood regions we wish to avoid
    sampling from.
    The solution is to detect these clusters and bound them in separate
    ellipsoids. For this, we use a recursive process where we perform
    K-means clustering with K=2. If the resulting two ellipsoids have a
    significantly lower total volume than the parent ellipsoid (less than half),
    we accept the split and repeat the clustering and volume test on each of
    the two subset of points. This process continues recursively.
    Alternatively, if the total ellipse volume is significantly greater
    than expected (based on the expected density of points) this indicates
    that there may be more than two clusters and that K=2 was not an
    appropriate cluster division.
    We therefore still try to subdivide the clusters recursively. However,
    we still only accept the final split into N clusters if the total volume
    decrease is significant.

    """
    # ----------------------- Preparation/Formatting --------------------------

    nparams = len(init)
    
    if grid_param_list is not None:
        if model_grid is None and model_reader is None:
            msg = "Either model_grid or model_reader have to be provided"
            raise TypeError(msg)
        n_gparams = len(grid_param_list)
        gp_dims = []
        for nn in range(n_gparams):
            gp_dims.append(len(grid_param_list[nn]))
        gp_dims = tuple(gp_dims)
    else:
        n_gparams = 0
        
    # format emission line dictionary and em grid
    if len(em_grid)>0:
        n_em = len(em_grid)
        em_dims = []
        for lab in labels:
            if lab in em_grid.keys():
                em_dims.append(len(em_grid[lab]))
        em_dims = tuple(em_dims)
        # update the grids depending on input units => make it surface flux
        idx_R = labels.index('R')
        for key, val in em_lines.items():
            if val[1] == 'L':
                idx_line = labels.index(key)
                # adapt grid
                R_si = init[idx_R]*con.R_jup.value
                conv_fac = 4*np.pi*R_si**2
                tmp = np.array(em_grid[key])/conv_fac
                em_grid[key] = tmp.tolist()
                # adapt ini state 
                init = list(init)
                init[idx_line] /= conv_fac
                init = tuple(init)
                #adapt bounds
                bounds_ori = list(bounds[key])
                bounds[key] = (bounds_ori[0]/conv_fac, bounds_ori[1]/conv_fac) 
            elif val[1] == 'LogL':
                idx_line = labels.index(key)
                R_si = init[idx_R]*con.R_jup.value
                conv_fac = con.L_sun.value/(4*np.pi*R_si**2)
                tmp = np.power(10,np.array(em_grid[key]))*conv_fac
                em_grid[key] = tmp.tolist()  
                # adapt ini state 
                init = list(init)
                init[idx_line] = conv_fac*10**init[idx_line]
                init = tuple(init)
                #adapt bounds
                bounds_ori = list(bounds[key])
                bounds[key] = (conv_fac*10**bounds_ori[0], 
                               conv_fac*10**bounds_ori[1]) 
            if em_lines[key][2] is not None:
                if em_lines[key][-1] == 'km/s':
                    v = em_lines[key][2]
                    em_lines_tmp = list(em_lines[key])
                    em_lines_tmp[2] = (1000*v/con.c.value)*em_lines[key][0]
                    em_lines_tmp[3] = 'mu'
                    em_lines[key] = tuple(em_lines_tmp)
                elif em_lines[key][-1] == 'nm':
                    em_lines_tmp = list(em_lines[key])
                    em_lines_tmp[2] = em_lines[key][2]/1000
                    em_lines_tmp[3] = 'mu'
                    em_lines[key] = tuple(em_lines_tmp)
                elif em_lines[key][-1] != 'mu':
                    msg = "Invalid unit of FWHM for line injection"
                    raise ValueError(msg)
                
    if model_grid is not None and grid_param_list is not None:
        if model_grid.ndim-2 != n_gparams:
            msg = "Ndim-2 of model_grid should match len(grid_param_list)"
            raise TypeError(msg)

 
    # Check model grid parameters extend beyond bounds to avoid extrapolation
    if grid_param_list is not None:
        for pp in range(n_gparams):
            if grid_param_list[pp][0]>bounds[labels[pp]][0]:
                msg= "Grid has to extend beyond bounds for {}."
                msg+="\n Consider increasing the lower bound to >{}."
                raise ValueError(msg.format(labels[pp],grid_param_list[pp][0]))
            if grid_param_list[pp][-1]<bounds[labels[pp]][1]:
                msg= "Grid has to extend beyond bounds for {}."
                msg+="\n Consider decreasing the upper bound to <{}."
                raise ValueError(msg.format(labels[pp],grid_param_list[pp][1]))
                
    # Check initial state is within bounds for all params (not only grid)
    for pp in range(nparams):
        if init[pp]<bounds[labels[pp]][0]:
            msg= "Initial state has to be within bounds for {}."
            msg+="\n Consider decreasing the lower bound to <{}."
            raise ValueError(msg.format(labels[pp],init[pp]))            
        if init[pp]>bounds[labels[pp]][1]:
            msg= "Initial state has to be within bounds for {}."
            msg+="\n Consider decreasing the upper bound to >{}."
            raise ValueError(msg.format(labels[pp],init[pp]))
        
    # Prepare model grid: convolve+resample models as observations 
    if resamp_before and grid_param_list is not None:
        if isfile(output_dir+grid_name):
            model_grid = open_fits(output_dir+grid_name)
            # check its shape is consistent with grid_param_list
            if model_grid.shape[:n_gparams] != gp_dims:
                msg="the loaded model grid ({}) doesn't have expected dims ({})"
                raise TypeError(msg.format(model_grid.shape,gp_dims))
            elif model_grid.shape[-2] != len(lbda_obs):
                msg="the loaded model grid doesn't have expected WL dimension"
                raise TypeError(msg)
            elif model_grid.shape[-1] != 2:
                msg="the loaded model grid doesn't have expected last dimension"
                raise TypeError(msg)
            elif len(em_grid) > 0:
                if model_grid.shape[n_gparams:n_gparams+n_em] != em_dims:
                    msg="loaded model grid ({}) doesn't have expected dims ({})"
                    raise TypeError(msg.format(model_grid.shape,em_dims))
        else:
            model_grid = make_resampled_models(lbda_obs, grid_param_list, 
                                               model_grid, model_reader, 
                                               em_lines, em_grid, dlbda_obs, 
                                               instru_res, instru_idx, 
                                               filter_reader, interp_nonexist)
            if output_dir and grid_name:
                write_fits(output_dir+grid_name, model_grid)
        # note: if model_grid is provided, it is still resampled to the 
        # same wavelengths as observed spectrum. However, if a fits name is 
        # provided in grid_name and that file exists, it is assumed the model 
        # grid in it is already resampled to match lbda_obs.    


    # -------------- Definition of utilities for nested sampling --------------

    def prior_transform(x):
        """
        Computes the transformation from the unit distribution `[0, 1]` to 
        parameter space. Uniform distributions are assumed for all parameters.
        
        
        Notes
        -----
        The prior transform function is used to specify the Bayesian prior for 
        the problem, in a round-about way. It is a transformation from a space 
        where variables are independently and uniformly distributed between 0 
        and 1 to the parameter space of interest. For independent parameters, 
        this would be the product of the inverse cumulative distribution 
        function (also known as the percent point function or quantile function) 
        for each parameter.
        http://kbarbary.github.io/nestle/prior.html
        """
        # uniform priors
        pt = []
        for p in range(nparams):            
            pmin = bounds[labels[p]][0]
            pmax = bounds[labels[p]][1]
            pt.append(x[p] * (pmax - pmin) + pmin)
                
        if priors is not None:
            # replace with Gaussian prior where relevant
            for key, prior in priors.items():
                if key == 'M' and 'logg' in labels:
                    msg = "Mass prior only available for MCMC sampling"
                    raise ValueError(msg)
                else:
                    idx_prior = labels.index(key)
                    pt[idx_prior] = prior[0] + prior[1] * ndtri(x)
            
        return np.array(pt)

    def f(param):
        
        return lnlike(param, labels, grid_param_list, lbda_obs, spec_obs, 
                      err_obs, dist, model_grid=model_grid,
                      model_reader=model_reader, em_lines=em_lines, 
                      em_grid=em_grid, dlbda_obs=dlbda_obs, 
                      instru_corr=instru_corr, instru_res=instru_res, 
                      instru_idx=instru_idx, use_weights=use_weights,
                      filter_reader=filter_reader, AV_bef_bb=AV_bef_bb, 
                      units_obs=units_obs, units_mod=units_mod, 
                      interp_order=interp_order)
        
        
    # ------------------------ Actual sampling --------------------------------
    if verbose:  start = time_ini()

    if verbose:
        print('Prior bounds on parameters:')
        for p in range(nparams):            
            pmin = bounds[labels[p]][0]
            pmax = bounds[labels[p]][1]
            print('{} [{},{}]'.format(labels[p], pmin, pmax))
        print('\nUsing {} active points'.format(npoints))


    if sampler == 'nestle':
        final_res = nestle.sample(f, prior_transform, ndim=nparams, 
                                  method=method, npoints=npoints, dlogz=dlogz, 
                                  **kwargs)
    elif sampler == 'ultranest':
        if isinstance(labels, tuple):
            labels = list(labels)
        un_sampler = ultranest.ReactiveNestedSampler(labels, f, prior_transform,
                                                     log_dir=output_dir, 
                                                     resume=True)
        res = un_sampler.run(min_num_live_points=npoints, dlogz=dlogz, **kwargs)
        final_res = (un_sampler, res)

    if verbose:
        print('\nTotal running time:')
        timing(start)
        
    return final_res


def show_nestle_results(ns_object, labels, method, burnin=0.4, bins=None, 
                        cfd=68.27, units=None, ndig=None, labels_plot=None, 
                        save=False, output_dir='/', plot=False,  **kwargs):
    """ 
    Show the results obtained with the Nestle sampler: summary, parameters with 
    errors, walk and corner plots. Returns best-fit values and uncertatinties.
    
    Parameters
    ----------
    ns_object: numpy.array
        The nestle object returned from `nested_spec_sampling`.
    labels: Tuple of strings
        Tuple of labels in the same order as initial_state, that is:
        - first all parameters related to loaded models (e.g. 'Teff', 'logg')
        - then the planet photometric radius 'R', in Jupiter radius
        - (optionally) the flux of emission lines (labels should match those \
        in the em_lines dictionary), in units of the model spectrum (times mu)
        - (optionally) the optical extinction 'Av', in mag
        - (optionally) the ratio of total to selective optical extinction 'Rv'
        - (optionally) 'Tbb1', 'Rbb1', 'Tbb2', 'Rbb2', etc. for each extra bb \
        contribution. 
    method : {"single", "multi", "classic"}, str optional
        Flavor of nested sampling.
    burnin: float, default: 0
        The fraction of a walker we want to discard.
    bins: int, optional
        The number of bins used to sample the posterior distributions.
    cfd: float, optional
        The confidence level given in percentage.
    units: tuple, opt
        Tuple of strings containing units for each parameter. If provided,
        mcmc_res will be printed on top of each 1d posterior distribution along 
        with these units.
    ndig: tuple, opt
        Number of digits precision for each printed parameter.
    labels_plot: tuple, opt
        Labels corresponding to parameter names, used for the plot. If None,
        will use "labels" passed in kwargs.
    save: boolean, default: False
        If True, a pdf file is created.
    output_dir: str, optional
        The name of the output directory which contains the output files in the 
        case  ``save`` is True. 
    plot: bool, optional
        Whether to show the plots (instead of saving them).
    kwargs:
        Additional optional arguments passed to `confidence` (matplotlib 
        optional arguments for histograms).
                    
    Returns
    -------
    final_res: numpy ndarray
         Best-fit parameters and uncertainties (corresponding to 68% confidence
         interval). Dimensions: nparams x 2.
    
    """
    res = ns_object
    nsamples = res.samples.shape[0]
    indburnin = int(np.percentile(np.array(range(nsamples)), burnin * 100))
    nparams = len(labels)

    print(res.summary())

    print(
        '\nNatural log of prior volume and Weight corresponding to each sample')
    if save or plot:
        plt.figure(figsize=(12, 4))
        plt.subplot(1, 2, 1)
        plt.plot(res.logvol, '.', alpha=0.5, color='gray')
        plt.xlabel('samples')
        plt.ylabel('logvol')
        plt.vlines(indburnin, res.logvol.min(), res.logvol.max(),
                   linestyles='dotted')
        plt.subplot(1, 2, 2)
        plt.plot(res.weights, '.', alpha=0.5, color='gray')
        plt.xlabel('samples')
        plt.ylabel('weights')
        plt.vlines(indburnin, res.weights.min(), res.weights.max(),
                   linestyles='dotted')
        if save:
            plt.savefig(output_dir+'Nestle-{}_results.pdf'.format(method))
        if plot:
            plt.show()
            
        print("\nWalk plots before the burnin")
        show_walk_plot(np.expand_dims(res.samples, axis=0), labels)
        if burnin > 0:
            print("\nWalk plots after the burnin")
            show_walk_plot(np.expand_dims(res.samples[indburnin:], axis=0),
                           labels)
        if save:
            plt.savefig(output_dir+'Nestle-{}_walk_plots.pdf'.format(method))
        
    mean, cov = nestle.mean_and_cov(res.samples[indburnin:],
                                    res.weights[indburnin:])
    print("\nWeighted mean +- sqrt(covariance)")
    if ndig is None:
        ndig = [3]*len(labels)   
    for p in range(nparams):
        fmt = "{{:.{0}f}}".format(ndig[p]).format
        print(r"{0} = {1} +/- {2}".format(labels[p], fmt(mean[p]), 
                                          fmt(np.sqrt(cov[p, p]))))
    if save:
        with open(output_dir+'Nestle-{}_sampling.txt'.format(method), "w") as f:
            f.write('#################################\n')
            f.write('####   CONFIDENCE INTERVALS   ###\n')
            f.write('#################################\n')
            f.write(' \n')
            f.write('Results of the NESTLE SAMPLING fit\n')
            f.write('----------------------------------\n ')
            f.write(' \n')
            f.write("\nWeighted mean +- sqrt(covariance)\n")
            for p in range(nparams):
                fmt = "{{:.{0}f}}".format(ndig[p]).format
                f.write(r"{0} = {1} +/- {2}\n".format(labels[p], fmt(mean[p]), 
                                                      fmt(np.sqrt(cov[p, p]))))
                        
    final_res = np.zeros([nparams,2]) 
    for p in range(nparams):     
        final_res[p] = [mean[p], np.sqrt(cov[p, p])]

    if bins is None:
        bins = int(np.sqrt(res.samples[indburnin:].shape[0]))
        print("\nHist bins =", bins)
    
    if save or plot:
        show_corner_plot(res.samples[indburnin:], burnin=burnin, save=save, 
                         output_dir=output_dir, mcmc_res=final_res, units=units, 
                         ndig=ndig, labels_plot=labels_plot, 
                         plot_name='Nestle-{}_corner_plot.pdf'.format(method), 
                         labels=labels)
    if save:
        plt.savefig(output_dir+'Nestle-{}_corner.pdf'.format(method))
    if plot:
        plt.show()
            
    print('\nConfidence intervals')
    if save or plot:
        _ = confidence(res.samples[indburnin:], labels, cfd=cfd, bins=bins,
                       weights=res.weights[indburnin:],
                       gaussian_fit=True, verbose=True, save=False, **kwargs)        
     
    if save:
        fn = output_dir+'Nestle-{}_confi_hist_gaussfit.pdf'
        plt.savefig(fn.format(method))

    if plot:
        plt.show()
        
    return final_res


def show_ultranest_results(un_object, labels, grid_param_list, dist, bins=None, 
                           cfd=68.27, ndig=None, save=False, output_dir='/', 
                           plot=False, col='b', units=None, labels_plot=None,
                           lbda_obs=None, spec_obs=None, spec_obs_err=None, 
                           n_pts=None, title=None, figsize=(8,6), **kwargs):
    """ 
    Shows the results obtained with the Ultranest sampler: summary, parameters 
    with errors, walk and corner plots. Returns best-fit values and 
    uncertatinties.
    
    Parameters
    ----------
    un_object: object
        The UltraNest Sampler object returned by nested_spec_sampling.
    labels: Tuple of strings
        Tuple of labels in the same order as initial_state, that is:
        - first all parameters related to loaded models (e.g. 'Teff', 'logg')
        - then the planet photometric radius 'R', in Jupiter radius
        - (optionally) the flux of emission lines (labels should match those \
        in the em_lines dictionary), in units of the model spectrum (times mu)
        - (optionally) the optical extinction 'Av', in mag
        - (optionally) the ratio of total to selective optical extinction 'Rv'
        - (optionally) 'Tbb1', 'Rbb1', 'Tbb2', 'Rbb2', etc. for each extra bb \
        contribution. 
    grid_param_list : list of 1d numpy arrays/lists OR None
        - If list, should contain list/numpy 1d arrays with available grid of \
        model parameters. 
        - Set to None for a pure n-blackbody fit, n=1,2,...
        - Note1: model grids should not contain grids on radius and Av, but \
        these should still be passed in initial_state (Av optional).
        - Note2: for a combined grid model + black body, just provide \
        the grid parameter list here, and provide values for 'Tbbn' and 'Rbbn' \
        in initial_state, labels and bounds.
    dist :  float
        Distance in parsec, used for flux scaling of the models.
    bins: int, optional
        The number of bins used to sample the posterior distributions.
    cfd: float, optional
        The confidence level given in percentage.
    ndig: tuple, opt
        Number of digits precision for each printed parameter.
    save: boolean, default: False
        If True, a pdf file is created.
    output_dir: str, optional
        The name of the output directory which contains the output files in the 
        case  ``save`` is True. 
    plot: bool, optional
        Whether to show the best-fit model against data.
    col: str, optional
        Color used for data points.
    lbda_obs, spec_obs, spec_obs_err: 1d ndarrays, optional
        Must be provided if plot is set to True
    n_pts: None or int, optional
        If None, models will be sampled at the same resolution as measured
        spectrum for plot. If an integer, corresponds to the number of sampled 
        points (uniform sampling) between the first and last point of the
        measured spectrum.
    title: str, optional
        If plot is set to True, title of the plot.
    kwargs:
        Additional optional arguments passed to `nested_spec_sampling` - only
        used if plot is set to True and model parameters are evaluated.
                    
    Returns
    -------
    final_res: numpy ndarray
         Best-fit parameters and uncertainties (corresponding to 68% confidence
         interval). Dimensions: nparams x 2.
    
    """
    un_sampler, res = un_object
    
    # print results
    un_sampler.print_results()
    
    if save or plot:
        # show run
        un_sampler.plot_run()
        if save:
            plt.savefig(output_dir+'UltraNest_run.pdf')
        if plot:
            plt.show()
                    
        # show trace
        cornerplot
        un_sampler.plot_trace()
        if save:
            plt.savefig(output_dir+'UltraNest_trace.pdf')
        if plot:
            plt.show()
    
    burned_res = un_burning(res)        
    
    # show results as table
    df = pd.DataFrame(data=res['samples'], columns=res['paramnames'])
    df.describe()
    
    # plot the best fit
    if save or plot:
        if lbda_obs is None or spec_obs is None or spec_obs_err is None:
            msg = "For the plot to be made, lbda_obs, spec_obs and spec_obs_err"
            msg += " must be provided"
            raise ValueError(msg)
        plt.figure(figsize=figsize)
        if title is not None:
            plt.title(title)
        plt.xlabel(r'Wavelength ($\mu$m)')
        plt.ylabel(r'$\lambda F_{\lambda}$ (W m-2)')
        plt.errorbar(x=lbda_obs, y=lbda_obs*spec_obs, yerr=spec_obs_err,
                     marker='o', ls=' ', color=col)
        
        if n_pts is None:
            lbda_grid = lbda_obs.copy()
        else:
            lbda_grid = np.linspace(lbda_obs[0], lbda_obs[-1], n_pts)
        band = PredictionBand(lbda_grid)
        
        # go through the solutions
        for params in un_sampler.results['samples']:
            # compute for each time the model
            mod = make_model_from_params(params, labels, grid_param_list, dist, 
                                         lbda_obs=lbda_grid, **kwargs)
            band.add(lbda_grid*mod[1])
        
        band.line(color='k')
        # add 1 sigma quantile
        band.shade(color='k', alpha=0.3)
        # add wider quantile (0.01 .. 0.99)
        band.shade(q=0.49, color='gray', alpha=0.2)
        if save:
            plt.savefig(output_dir+"Best-fit_UltraNest.pdf",  
                        bbox_inches='tight')

    print('\nConfidence intervals')
    if bins is None:
        bins = int(np.sqrt(burned_res[0].shape[0]))
        print("\nHist bins =", bins)
    if save or plot:
        val_max, confidenceInterval = confidence(burned_res[0], labels, 
                                                 cfd=cfd, bins=bins, 
                                                 weights=burned_res[1],
                                                 gaussian_fit=False, 
                                                 verbose=True, save=save)  
        if save:
            plt.savefig(output_dir+'UltraNest_confi_hist.pdf')
        if plot:
            plt.show()
        mu, sig = confidence(burned_res[0], labels, cfd=cfd, bins=bins,
                             weights=burned_res[1], gaussian_fit=True, 
                             verbose=True, save=save)       
        if save:
            plt.savefig(output_dir+'UltraNest_confi_hist_gaussfit.pdf')     
        if plot:
            plt.show()
            
    nparams = len(labels)
    if ndig is None:
        ndig = [3]*len(labels)   
    for p in range(nparams):
        fmt = "{{:.{0}f}}".format(ndig[p]).format
        print(r"{0} = {1} +/- {2}".format(labels[p], fmt(mu[p]), 
                                          fmt(np.sqrt(sig[p]))))
    if save:
        with open(output_dir+'UltraNest_sampling.txt', "w") as f:
            f.write('#################################\n')
            f.write('####   CONFIDENCE INTERVALS   ###\n')
            f.write('#################################\n')
            f.write(' \n')
            f.write('Results of the NESTED SAMPLING (UltraNest)\n')
            f.write('------------------------------------------\n ')
            f.write(' \n')
            f.write("\nMost likely value [{}% confidence interval]\n".format(cfd))
            for p in range(nparams):
                fmt = "{{:.{0}f}}".format(ndig[p]).format
                line = r"{0} = {1} [{2}, {3}] \n"
                f.write(line.format(labels[p], fmt(val_max[labels[p]]), 
                                    fmt(confidenceInterval[labels[p]][0]),
                                    fmt(confidenceInterval[labels[p]][1])))
            f.write("\nGaussian fit results: Weighted mean +- sqrt(covariance) \n")
            for p in range(nparams):
                fmt = "{{:.{0}f}}".format(ndig[p]).format
                f.write(r"{0} = {1} +/- {2}\n".format(labels[p], fmt(mu[p]), 
                                                      fmt(sig[p])))                


    if plot or save:
        # show corner plot
        show_corner_plot(burned_res[0], save=save, units=units, ndig=ndig, 
                         mcmc_res=(val_max, confidenceInterval), 
                         labels_plot=labels_plot, output_dir=output_dir, 
                         plot_name='UltraNest-corner_plot.pdf', 
                         labels=labels, weights=burned_res[1])
        if plot:
            plt.show()
        
    final_res = val_max, confidenceInterval
        
    return final_res