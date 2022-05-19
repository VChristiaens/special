#! /usr/bin/env python

"""
Module to estimate the spectral correlation between channels of an IFS datacube.

.. [GRE16]
   | Greco & Brandt 2016
   | **The Measurement, Treatment, and Impact of Spectral Covariance and 
     Bayesian Priors in Integral-field Spectroscopy of Exoplanets**
   | *The Astrophysical Journal, Volume 833, Issue 1, p. 134*
   | `https://arxiv.org/abs/1602.00691
     <https://arxiv.org/abs/1602.00691>`_
"""

__author__ = 'V. Christiaens'
__all__ = ['spectral_correlation',
           'combine_spec_corrs']

from astropy.stats import gaussian_fwhm_to_sigma
import numpy as np
from scipy.optimize import curve_fit


def spectral_correlation(array, awidth=2, r_in=1, r_out=None, pl_xy=None,
                         mask_r=4, fwhm=4, sp_fwhm_guess=3, full_output=False):
    """ Computes the spectral correlation between (post-processed) IFS frames, 
    as a function of radius, implemented as Eq. 7 of [GRE16]. This is a crucial 
    step for an unbiased fit of a measured IFS spectrum to either synthetic or 
    template spectra.
    
    Parameters
    ----------
    array : numpy ndarray
        Input cube or 3d array, of dimensions n_ch x n_y x n_x; where n_y and 
        n_x should be odd values (star should be centered on central pixel).
    awidth : int, optional
        Width in pixels of the concentric annuli used to compute the spectral 
        correlation as a function of radial separation. Greco & Brandt 2017 
        noted no significant differences for annuli between 1 and 3 pixels 
        width on GPI data.
    r_in: int, optional
        Innermost radius where the spectral correlation starts to be computed.
    r_out: int, optional
        Outermost radius where the spectral correlation is computed. If left as 
        None, it will automatically be computed up to the edge of the frame. 
    pl_xy: tuple of tuples of 2 floats, optional
        x,y coordinates of all companions present in the images. If provided, 
        a circle centered on the location of each companion will be masked out 
        for the spectral correlation computation.
    mask_r: float, optional
        if pl_xy is provided, this should also be provided. Size of the 
        aperture around each companion (in terms of fwhm) that is discarded to 
        not bias the spectral correlation computation.
    fwhm: float, optional
        if pl_xy is provided, this should also be provided. By default we  
        consider a 2FWHM aperture mask around each companion to not bias the 
        spectral correlation computation.
    sp_fwhm_guess: float, optional
        Initial guess on the spectral FWHM of all channels.
    full_output: bool, opt
        Whether to also output the fitted spectral FWHM for each channel, and 
        the vector of radial separation at which each spectral correlation
        matrix is calculated.
    Note: radii that are skipped will be filled with zeros in the output cube.

    Returns
    -------
    sp_corr : numpy ndarray
        3d array of spectral correlation, as a function of radius with 
        dimensions: n_rad x n_ch x n_ch, where n_rad = int((r_out-r_in)/2)
    sp_fwhm: numpy ndarray
        (if full_output is True) 2d array containing the spectral fwhm at each 
        radius, for each spectral channel. Dims: n_rad x n_ch
    sp_rad: numpy ndarray
        (if full_output is True) 1d array containing the radial separation of
        each measured spectral correlation matrix. Dims: n_rad    
    """

    if not isinstance(awidth,int) or not isinstance(r_in,int):
        raise TypeError("Inputs should be integers")

    if array.ndim != 3:
        raise TypeError("Input array should be 3D.")
        
    n_ch, n_y, n_x = array.shape
    n_r = min((n_y-1)/2.,(n_x-1)/2.)
    if n_r%1:
        raise TypeError("Input array y and x dimensions should be odd")
    
    if r_out is None:
        r_out = n_r

    test_rads = np.arange(r_in-1,r_out-1)
    n_rad = max(1,int(np.floor(test_rads.shape[0]/awidth)))
    
    #n_rad = int(np.ceil(n_r/ann_width)) # effective number of annuli probed
    
    sp_corr = np.zeros([n_rad,n_ch,n_ch])
    sp_rad= np.zeros([n_rad])
    if full_output:
        sp_fwhm = np.zeros([n_rad,n_ch])
        def gauss_1fp(x, *p):
            sigma = p[0]*gaussian_fwhm_to_sigma
            return np.exp(-x**2/(2.*sigma**2))
    mask_f = np.zeros_like(array[0])

    if pl_xy is not None:
        mask = np.ones_like(array[0])
        for i in range(len(pl_xy)):
            if not isinstance(pl_xy[i], tuple):
                raise TypeError("Format of companions coordinates incorrect")
            mask_i = get_circle(mask, radius=mask_r*fwhm, cy=pl_xy[i][1], 
                                cx=pl_xy[i][0], mode="mask")
            mask_f[np.where(mask_i)] = 1

    for ann in range(n_rad):
        inner_radius = r_in+ (ann * awidth)
        ind = get_annulus_segments(array[0], inner_radius, awidth)
        yy = ind[0][0]
        xx = ind[0][1]
        yy_final = [yy[i] for i in range(len(ind[0][0])) if not mask_f[yy[i],
                                                                       xx[i]]]
        xx_final = [xx[i] for i in range(len(ind[0][0])) if not mask_f[yy[i],
                                                                       xx[i]]]
        matrix = array[:, yy_final, xx_final]  # shape (z, npx_annsegm)
        sp_rad[ann*awidth:(ann+1)*awidth] = r_in+(ann+0.5)*awidth
        
        for zi in range(n_ch):
            for zj in range(n_ch):
                num = np.nanmean(matrix[zi]*matrix[zj])
                denom = np.sqrt(np.nanmean(matrix[zi]*matrix[zi])* \
                                np.nanmean(matrix[zj]*matrix[zj]))
                sp_corr[ann*awidth:(ann+1)*awidth,zi,zj] = num/denom
            if full_output:
                p0 = (sp_fwhm_guess,)
                x = np.arange(n_ch)-zi
                y = sp_corr[ann*awidth,zi]# norm y
                y = y-np.amin(y)
                y = y/np.amax(y)
                coeff, var_matrix = curve_fit(gauss_1fp, x, y, p0=p0)
                sp_fwhm[ann*awidth:(ann+1)*awidth,zi] = coeff[0]
                
    # Zero is adopted for uncorrelated channels
    sp_corr[np.where(sp_corr<0)] = 0
                
    if full_output:
        return sp_corr, sp_fwhm, sp_rad
    else:
        return sp_corr
    

def combine_spec_corrs(arr_list):
    """ Combines the spectral correlation matrices of different instruments 
    into a single square matrix (required for input of spectral fits).
    
    Parameters
    ----------
    arr_list : list or tuple of numpy ndarrays
        List/tuple containing the distinct square spectral correlation matrices 
        OR ones (for independent photometric measurements). 

    Returns
    -------
    combi_corr : numpy 2d ndarray
        2d square ndarray representing the combined spectral correlation.
        
    """
    n_arr = len(arr_list)
    
    size = 0
    for nn in range(n_arr):
        if isinstance(arr_list[nn],np.ndarray):
            if arr_list[nn].ndim != 2:
                raise TypeError("Arrays of the tuple should be 2d")
            elif arr_list[nn].shape[0] != arr_list[nn].shape[1]:
                raise TypeError("Arrays of the tuple should be square")
            size+=arr_list[nn].shape[0]
        elif arr_list[nn] == 1:
            size+=1
        else:
            raise TypeError("Tuple can only have square 2d arrays or ones")
            
    combi_corr = np.zeros([size,size])
    
    size_tmp = 0
    for nn in range(n_arr):
        if isinstance(arr_list[nn],np.ndarray):
            mm = arr_list[nn].shape[0]
            combi_corr[size_tmp:size_tmp+mm,size_tmp:size_tmp+mm]=arr_list[nn]
            size_tmp+=mm
        elif arr_list[nn] == 1:
            combi_corr[size_tmp,size_tmp]=1
            size_tmp+=1      
        
    return combi_corr

    
    
def get_annulus_segments(data, inner_radius, width, nsegm=1, theta_init=0,
                         optim_scale_fact=1, mode="ind"):
    """
    Return indices or values in segments of a centerered annulus (as in VIP).

    The annulus is defined by ``inner_radius <= annulus < inner_radius+width``.

    Parameters
    ----------
    data : 2d numpy ndarray or tuple
        Input 2d array (image) or tuple with its shape.
    inner_radius : float
        The inner radius of the donut region.
    width : float
        The size of the annulus.
    nsegm : int
        Number of segments of annulus to be extracted.
    theta_init : int
        Initial azimuth [degrees] of the first segment, counting from the
        positive x-axis counterclockwise.
    optim_scale_fact : float
        To enlarge the width of the segments, which can then be used as
        optimization segments (e.g. in LOCI).
    mode : {'ind', 'val', 'mask'}, optional
        Controls what is returned: indices of selected pixels, values of
        selected pixels, or a boolean mask.

    Returns
    -------
    indices : list of ndarrays
        [mode='ind'] Coordinates of pixels for each annulus segment.
    values : list of ndarrays
        [mode='val'] Pixel values.
    masked : list of ndarrays
        [mode='mask'] Copy of ``data`` with masked out regions.

    Notes
    -----
    Moving from ``get_annulus`` to ``get_annulus_segments``:

    .. code::python
        # get_annulus handles one single segment only, so note the ``[0]`` 
        after the call to get_annulus_segments if you want to work with one 
        single segment only.

        get_annulus(arr, 2, 3, output_indices=True)
        # is the same as
        get_annulus_segments(arr, 2, 3)[0]

        get_annulus(arr, inr, w, output_values=True)
        # is the same as
        get_annulus_segments(arr, inr, w, mode="val")[0]

        get_annulus(arr, inr, w)
        # is the same as
        get_annulus_segments(arr, inr, w, mode="mask")[0]

        # the only difference is the handling of the border values:
        # get_annulus_segments is `in <= ann < out`, while get_annulus is
        # `in <= ann <= out`. But that should make no difference in practice.

    """
    array = frame_or_shape(data)

    if not isinstance(nsegm, int):
        raise TypeError('`nsegm` must be an integer')

    cy, cx = frame_center(array)
    azimuth_coverage = np.deg2rad(int(np.ceil(360 / nsegm)))
    twopi = 2 * np.pi

    yy, xx = np.mgrid[:array.shape[0], :array.shape[1]]
    rad = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    phi = np.arctan2(yy - cy, xx - cx)
    phirot = phi % twopi
    outer_radius = inner_radius + (width*optim_scale_fact)
    masks = []

    for i in range(nsegm):
        phi_start = np.deg2rad(theta_init) + (i * azimuth_coverage)
        phi_end = phi_start + azimuth_coverage

        if phi_start < twopi and phi_end > twopi:
            masks.append((rad >= inner_radius) & (rad < outer_radius) &
                         (phirot >= phi_start) & (phirot <= twopi) |
                         (rad >= inner_radius) & (rad < outer_radius) &
                         (phirot >= 0) & (phirot < phi_end - twopi))
        elif phi_start >= twopi and phi_end > twopi:
            masks.append((rad >= inner_radius) & (rad < outer_radius) &
                         (phirot >= phi_start - twopi) &
                         (phirot < phi_end - twopi))
        else:
            masks.append((rad >= inner_radius) & (rad < outer_radius) &
                         (phirot >= phi_start) & (phirot < phi_end))

    if mode == "ind":
        return [np.where(mask) for mask in masks]
    elif mode == "val":
        return [array[mask] for mask in masks]
    elif mode == "mask":
        return [array*mask for mask in masks]
    else:
        raise ValueError("mode '{}' unknown!".format(mode))
        
        
def get_circle(array, radius, cy=None, cx=None, mode="mask"):
    """
    Return a centered circular region from a 2d ndarray (as in VIP).

    Parameters
    ----------
    array : numpy ndarray
        Input 2d array or image.
    radius : int
        The radius of the circular region.
    cy, cx : int, optional
        Coordinates of the circle center. If one of them is ``None``, the 
        center of ``array`` is used.
    mode : {'mask', 'val'}, optional
        Controls what is returned: array with circular mask applied, or values
        of the pixels in the circular region.

    Returns
    -------
    masked : numpy ndarray
        [mode="mask"] Input array with the circular mask applied.
    values : numpy ndarray
        [mode="val"] 1d array with the values of the pixels in the circular
        region.

    Notes
    -----
    An alternative implementation would use ``skimage.draw.disk``. ``disk``
    performs better on large ``array``s (e.g. 1000px, 10.000px), while the
    current implementation is faster for small ``array``s (e.g. 100px). See
    `test_shapes.py` for benchmark details.

    """
    if array.ndim != 2:
        raise TypeError('Input array is not a frame or 2d array.')
    sy, sx = array.shape
    if cy is None or cx is None:
        cy, cx = frame_center(array, verbose=False)

    # ogrid is a multidim mesh creator (faster than mgrid):
    yy, xx = np.ogrid[:sy, :sx]
    circle = (yy - cy) ** 2 + (xx - cx) ** 2  # eq of circle. sq dist to center
    circle_mask = circle < radius ** 2  # boolean mask
    if mode == "mask":
        return array * circle_mask
    elif mode == "val":
        return array[circle_mask]
    else:
        raise ValueError("mode '{}' unknown!".format(mode))
        
        
def frame_center(array, verbose=False):
    """
    Return the coordinates y,x of the frame(s) center.
    If odd: dim/2-0.5
    If even: dim/2

    Parameters
    ----------
    array : 2d/3d/4d numpy ndarray
        Frame or cube.
    verbose : bool optional
        If True the center coordinates are printed out.

    Returns
    -------
    cy, cx : int
        Coordinates of the center.

    """
    if array.ndim == 2:
        shape = array.shape
    elif array.ndim == 3:
        shape = array[0].shape
    elif array.ndim == 4:
        shape = array[0, 0].shape
    else:
        raise ValueError('`array` is not a 2d, 3d or 4d array')

    cy = shape[0] / 2
    cx = shape[1] / 2

    if shape[0]%2:
        cy-=0.5
    if shape[1]%2:
        cx-=0.5        

    if verbose:
        print('Center px coordinates at x,y = ({}, {})'.format(cx, cy))  
    
    return int(cy), int(cx)



def frame_or_shape(data):
    """
    Sanitize ``data``, always return a 2d frame.

    If ``data`` is a 2d frame, it is returned unchanged. If it is a shaped,
    return an empty array of that shape.

    Parameters
    ----------
    data : 2d ndarray or shape tuple

    Returns
    -------
    array : 2d ndarray

    """
    if isinstance(data, np.ndarray):
        array = data
        if array.ndim != 2:
            raise TypeError('`data` is not a frame or 2d array')
    elif isinstance(data, tuple):
        array = np.zeros(data)
    else:
        raise TypeError('`data` must be a tuple (shape) or a 2d array')

    return array