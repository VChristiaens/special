#! /usr/bin/env python

"""
Module with utilities to estimate the spectral type and gravity of an object 
based on spectral indices.
"""

__author__ = 'V. Christiaens'
__all__ = ['spectral_idx',
           'sp_idx_to_spt',
           'idx_to_spt',
           'spt_to_idx']

import numpy as np
from .utils_spec import find_nearest

def spectral_idx(lbda, spec, line='H2O-1.5', spec_err=None, verbose=False):
    """
    Computes a spectral index. Implemented so far:
        - the H2O 1.5 mu index (Allers et al. 2007)
        - the Na 1.1 mu index (Allers et al. 2007)
        
    Parameters
    ----------
    lbda : numpy ndarray
        1d numpy array containing the wavelengths of the spectrum in microns.
    spec : numpy ndarray
        1d numpy array containing the measured flux (arbitrary units accepted).
    line: str, optional {'H2O-1.5', 'Na-1.1'}
        Which line to use for the calculation of the spectral index.
    spec_err: numpy ndarray, optional
        1d numpy array containing the uncertainties on the measured flux 
        (arbitrary units accepted). If provided the uncertainty on the spectral
        index will also be returned.
    verbose: bool, optional
        Whther to print more information.

    Returns
    -------
    index : float
        Value of the spectral index
    index_err: float
       [if spec_err is provided] Uncertainty on the spectral index.
    """
    
    if line == 'H2O-1.5':
        lbda_num0 = 1.55
        lbda_num1 = 1.56
        lbda_den0 = 1.492
        lbda_den1 = 1.502
    elif line == 'Na-1.1':
        lbda_num0 = 1.15
        lbda_num1 = 1.16
        lbda_den0 = 1.134
        lbda_den1 = 1.144
    else:
        raise ValueError("line not recognised")
        
    idx_num0 = find_nearest(lbda, lbda_num0)
    idx_num1 = find_nearest(lbda, lbda_num1)
    idx_den0 = find_nearest(lbda, lbda_den0)
    idx_den1 = find_nearest(lbda, lbda_den1)
    
    if verbose:
        msg = "Closest indices to {}, {}, {} and {} mu".format(lbda_num0,lbda_num1,lbda_den0,lbda_den1)
        msg+= "correspond to {}, {}, {} and {} resp.".format(lbda[idx_num0],lbda[idx_num1],lbda[idx_den0],lbda[idx_den1])
        print(msg)
    
    if idx_den0 == idx_num1:
        raise ValueError("Same index for num1 and den0: not enough spec res")
    
    num = np.mean(spec[idx_num0:idx_num1+1])
    den = np.mean(spec[idx_den0:idx_den1+1])
    index = num/den
    
    if spec_err is None:
        return index
    else:
        num_err = np.sum(np.power(spec_err[idx_num0:idx_num1+1]/((idx_num1+1-idx_num0)*den),2))
        den_err = np.sum(np.power(spec_err[idx_num0:idx_num1+1]*(idx_den1+1-idx_den0)*num/(den**2),2))
        index_err = np.sqrt(num_err**2+den_err**2)
        return index, index_err


def sp_idx_to_spt(idx, line='H2O-1.5', idx_err=None):
    """
    Estimates a spectral type from a spectral index. Implemented so far:
        - the H2O 1.5 mu index (Allers et al. 2007)
        
    Note on scale of SpT:
        5 = M5
        15 = L5
        
    Parameters
    ----------
    idx : float
        Value of spectral index
    line: str, optional {'H2O-1.5'}
        The line corresponding to the spectral index.
    idx_err: float, optional
        Uncertainty on the spectral index value

    Returns
    -------
    spt : float
        Value of the spectral type
    spt_err: float
       [if idx_err is provided] Uncertainty on the spectral type.

    """

    if line == 'H2O-1.5':
        a = 0.75
        da = 0.03
        b = 0.044
        db = 0.004
    else:
        raise ValueError("not implemented for this line")
        
    spt = (idx-a)/b
    tot_err = np.power(da/b,2) # err related to a
    tot_err += np.power(db*(idx-a)/b**2,2)  # err related to b
    if idx_err is not None:
        tot_err += np.power(idx_err/b, 2)
    spt_err = np.sqrt(tot_err)
    
    return spt, spt_err


def spt_to_idx(spt, verbose=False):
    """
    Converts a string representing spectral type into an integer index.
    Convention (from splat): K0 = 0, M0=10, L0=20, T0=30, Y9 = 49
    
    Parameters
    ----------
    spt : str
        String representing the spectral index
    verbose: bool, optional
        Whther to print more information.

    Returns
    -------
    idx : float or int
        Index value of the spectral type
    """

    spt = spt.split()[0]
    
    if 'K' in spt:
        idx = 0
        loc = spt.find('K')
    elif 'M' in spt:
        idx = 10
        loc = spt.find('M')
    elif 'L' in spt:
        idx = 20
        loc = spt.find('L')
    elif 'T' in spt:
        idx = 30
        loc = spt.find('T')
    elif 'Y' in spt:
        idx = 40
        loc = spt.find('Y')
    else:
        idx = np.nan
        if verbose:
            print("Spectral type not between K0 and Y9. Idx set to nan.")

    if len(spt)>1: # ie. not just "M" or "L"
        if len(spt[loc+1:]) > 2:
            if spt[loc+2] == ".":
                idx+=float(spt[loc+1:loc+4])
        elif spt[loc+1].isdigit():
            idx+=float(spt[loc+1])
    else:
        idx+=5 # in case just "M" or "L" is given, let's assume it's M5 or L5
    
    return idx

    
def idx_to_spt(idx, convention='splat'):
    """
    Converts an integer index into spectral type.
    Convention (from splat): K0 = 0, M0=10, L0=20, T0=30, Y9 = 49
    Convention (from Allers+07): M5 = 5, L5 = 15
    
    Parameters
    ----------
    idx : float or int
        Index value of the spectral type 
    verbose: bool, optional
        Whther to print more information.

    Returns
    -------
    spt : str
        String representing the spectral index
    """
    
    if convention == 'splat':
        if idx < 10:
            spt = 'K'
        elif idx < 20:
            spt = 'M'
        elif idx < 30:
            spt = 'L'
        elif idx < 40:
            spt = 'T'
        elif idx < 50:
            spt = 'Y'
        else:
            raise ValueError("Index not between 0 (K0) and 49 (Y9)")
    elif convention == 'Allers07':
        if idx < 10:
            spt = 'M'
        elif idx < 20:
            spt = 'L'
        else:
            raise ValueError("Index not between 0 (M0) and 20 (L9)")
        
    if not (idx%10)%1:
       spt+="{:.0f}".format(idx%10)
    else:
        spt+="{:.1f}".format(idx%10)
    
    return spt
    