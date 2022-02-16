#! /usr/bin/env python

"""
Module with functions for posterior sampling of the model spectra parameters 
using nested sampling (``nestle``).
"""



__author__ = 'V. Christiaens, C. A. Gomez Gonzalez',
__all__ = ['nested_spec_sampling',
           'nested_sampling_results']

import nestle
import corner
import numpy as np
from matplotlib import pyplot as plt
from ..config import time_ini, timing
from .mcmc_sampling import lnlike, confidence, show_walk_plot


def nested_spec_sampling(init, lbda_obs, spec_obs, err_obs, dist, 
                         grid_param_list, labels, bounds, resamp_before=True, 
                         model_grid=None, model_reader=None, em_lines={}, 
                         em_grid={}, dlbda_obs=None, instru_corr=None, 
                         instru_fwhm=None, instru_idx=None, filter_reader=None, 
                         AV_bef_bb=False, units_obs='si', units_mod='si', 
                         interp_order=1, priors=None, physical=True, 
                         interp_nonexist=True, w=0.1, method='single', 
                         npoints=100, dlogz=0.1, decline_factor=None, 
                         rstate=None, verbose=True):
                            
    
    """ Runs a nested sampling algorithm in order to determine the position and
    the flux of the planet using the 'Negative Fake Companion' technique. The
    result of this procedure is a a ``nestle`` object containing the samples
    from the posterior distributions of each of the 3 parameters. It provides
    pretty good results (value plus error bars) compared to a more CPU intensive
    Monte Carlo approach with the affine invariant sampler (``emcee``).

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
        Spectral channel width for the observed spectrum. It should be provided 
        IF one wants to weigh each point based on the spectral 
        resolution of the respective instruments (as in Olofsson et al. 2016).
    instru_corr : numpy 2d ndarray or list, optional
        Spectral correlation throughout post-processed images in which the 
        spectrum is measured. It is specific to the combination of instrument, 
        algorithm and radial separation of the companion from the central star.
        Can be computed using distances.spectral_correlation(). In case of
        a spectrum obtained with different instruments, build it with
        distances.combine_corrs(). If not provided, it will consider the 
        uncertainties in each spectral channels are independent. See Greco & 
        Brandt (2017) for details.
    instru_fwhm : float OR list of either floats or strings, optional
        The instrumental spectral fwhm provided in nm. This is used to convolve
        the model spectrum. If several instruments are used, provide a list of 
        instru_fwhm values, one for each instrument whose spectral resolution
        is coarser than the model - including broad band filter FWHM if 
        relevant.
        If strings are provided, they should correspond to filenames (including 
        full paths) of text files containing the filter information for each 
        observed wavelength. Strict format: 
    instru_idx: numpy 1d array, optional
        1d array containing an index representing each instrument used 
        to obtain the spectrum, label them from 0 to n_instru. Zero for points 
        that don't correspond to any instru_fwhm provided above, and i in 
        [1,n_instru] for points associated to instru_fwhm[i-1]. This parameter 
        must be provided if the spectrum consists of points obtained with 
        different instruments.
    filter_reader: python routine, optional
        External routine that reads a filter file and returns a 2D numpy array, 
        where the first column corresponds to wavelengths, and the second 
        contains transmission values. Important: if not provided, but strings 
        are detected in instru_fwhm, the default file reader will be used. 
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
    method : {"single", "multi", "classic"}, str optional
        Flavor of nested sampling. Single ellipsoid works well for the NEGFC and
        is the default.
    npoints : int optional
        Number of active points. At least ndim+1 (4 will produce bad results).
        For problems with just a few parameters (<=5) like the NEGFC, good
        results are obtained with 100 points (default).
    dlogz : Estimated remaining evidence
        Iterations will stop when the estimated contribution of the remaining
        prior volume to the total evidence falls below this threshold.
        Explicitly, the stopping criterion is log(z + z_est) - log(z) < dlogz
        where z is the current evidence from all saved samples, and z_est is the
        estimated contribution from the remaining volume. This option and
        decline_factor are mutually exclusive. If neither is specified, the
        default is dlogz=0.5.
    decline_factor : float, optional
        If supplied, iteration will stop when the weight (likelihood times prior
        volume) of newly saved samples has been declining for
        decline_factor * nsamples consecutive samples. A value of 1.0 seems to
        work pretty well.
    rstate : random instance, optional
        RandomState instance. If not given, the global random state of the
        numpy.random module will be used.

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

    # If companion flux is too low MCMC will not converge. Solution: scale up 
    # the intensities in the cube after injecting the negfc.
    if init[2] < 100:
        scale_fac = 100./init[2]
    else:
        scale_fac = 1

    def prior_transform(x):
        """
        Computes the transformation from the unit distribution `[0, 1]` to 
        parameter space. Uniform distributions are assumed for all parameters.
        
        The default distributions used are
        radius: Uniform distribution transformed into polar coordinates
            This distribution assumes uniform distribution for the (x,y) 
            coordinates transformed to polar coordinates.
        theta: Uniform distribution
            This distribution is derived the same as the radial distribution, 
            but there is no change on the prior for theta after the 
            change-of-variables transformation.
        flux: Poisson-invariant scale distribution
            This distribution is the Jeffrey's prior for Poisson data
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
        nparams = len(init)
        
        if np.isscalar(w):
            w_list = [w]*nparams
        elif nparams != len(w):
            msg = "length of w should be equal to the number of parameters"
            raise ValueError(msg)
        else:
            w_list = w
            
        pt = []
        for p in range(nparams):            
            tmin = init[p] * (1-w_list[p])
            tmax = init[p] * (1+w_list[p])
            pt.append(x[p] * (tmax - tmin) + tmin)
        
        return np.array(pt)

    def f(param):
        return lnlike(params, labels, grid_param_list, lbda_obs, spec_obs, err_obs, dist, 
                   model_grid=None, model_reader=None, em_lines={}, em_grid={}, 
                   dlbda_obs=None, instru_corr=None, instru_fwhm=None, instru_idx=None, 
                   filter_reader=None, AV_bef_bb=False, units_obs='si', units_mod='si', 
                   interp_order=1
            
            param=param, cube=cube, angs=angs, plsc=plsc,
                      psf_norm=psf, fwhm=fwhm, annulus_width=annulus_width,
                      aperture_radius=aperture_radius, initial_state=init,
                      cube_ref=cube_ref, svd_mode=svd_mode, scaling=scaling,
                      algo=algo, delta_rot=delta_rot, fmerit='sum', ncomp=ncomp, 
                      collapse=collapse, algo_options=algo_options, 
                      weights=weights, scale_fac=scale_fac)

    # -------------------------------------------------------------------------
    if verbose:  start = time_ini()

    if verbose:
        print('Prior bounds on parameters:')
        print('Radius [{},{}]'.format(init[0] - w[0], init[0] + w[0], ))
        print('Theta [{},{}]'.format(init[1] - w[1], init[1] + w[1]))
        print('Flux [{},{}]'.format(init[2] - w[2], init[2] + w[2]))
        print('\nUsing {} active points'.format(npoints))

    res = nestle.sample(f, prior_transform, ndim=3, method=method,
                        npoints=npoints, rstate=rstate, dlogz=dlogz,
                        decline_factor=decline_factor)

    # if verbose:  print; timing(start)
    if verbose:
        print('\nTotal running time:')
        timing(start)
    return res


def nested_sampling_results(ns_object, burnin=0.4, bins=None, save=False,
                            output_dir='/', plot=False):
    """ Shows the results of the Nested Sampling, summary, parameters with 
    errors, walk and corner plots.
    """
    res = ns_object
    nsamples = res.samples.shape[0]
    indburnin = int(np.percentile(np.array(range(nsamples)), burnin * 100))

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
        if plot:
            plt.show()
    
        plt.savefig(output_dir+'Nested_results.pdf')
            
        print("\nWalk plots before the burnin")
        show_walk_plot(np.expand_dims(res.samples, axis=0))
        if burnin > 0:
            print("\nWalk plots after the burnin")
            show_walk_plot(np.expand_dims(res.samples[indburnin:], axis=0))
        plt.savefig(output_dir+'Nested_walk_plots.pdf')
        
    mean, cov = nestle.mean_and_cov(res.samples[indburnin:],
                                    res.weights[indburnin:])
    print("\nWeighted mean +- sqrt(covariance)")
    print("Radius = {:.3f} +/- {:.3f}".format(mean[0], np.sqrt(cov[0, 0])))
    print("Theta = {:.3f} +/- {:.3f}".format(mean[1], np.sqrt(cov[1, 1])))
    print("Flux = {:.3f} +/- {:.3f}".format(mean[2], np.sqrt(cov[2, 2])))

    if save:
        with open(output_dir+'Nested_sampling.txt', "w") as f:
            f.write('#################################\n')
            f.write('####   CONFIDENCE INTERVALS   ###\n')
            f.write('#################################\n')
            f.write(' \n')
            f.write('Results of the NESTED SAMPLING fit\n')
            f.write('----------------------------------\n ')
            f.write(' \n')
            f.write("\nWeighted mean +- sqrt(covariance)\n")
            f.write("Radius = {:.3f} +/- {:.3f}\n".format(mean[0], np.sqrt(cov[0, 0])))
            f.write("Theta = {:.3f} +/- {:.3f}\n".format(mean[1], np.sqrt(cov[1, 1])))
            f.write("Flux = {:.3f} +/- {:.3f}\n".format(mean[2], np.sqrt(cov[2, 2])))
                        
    if bins is None:
        bins = int(np.sqrt(res.samples[indburnin:].shape[0]))
        print("\nHist bins =", bins)
    
    if save or plot:
        ranges = None
        fig = corner.corner(res.samples[indburnin:], bins=bins,
                            labels=["$r$", r"$\theta$", "$f$"],
                            weights=res.weights[indburnin:], range=ranges,
                            plot_contours=True)
        fig.set_size_inches(8, 8)
    if save:
        plt.savefig(output_dir+'Nested_corner.pdf')
            
    print('\nConfidence intervals')
    if save or plot:
        _ = confidence(res.samples[indburnin:], cfd=68, bins=bins,
                       weights=res.weights[indburnin:],
                       gaussian_fit=True, verbose=True, save=False)
                   
    if save:
        plt.savefig(output_dir+'Nested_confi_hist_flux_r_theta_gaussfit.pdf')

    final_res = np.array([[mean[0], np.sqrt(cov[0, 0])],
                          [mean[1], np.sqrt(cov[1, 1])],
                          [mean[2], np.sqrt(cov[2, 2])]])
    return final_res