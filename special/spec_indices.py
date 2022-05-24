#! /usr/bin/env python

"""
Module with utilities to estimate the spectral type and gravity of an object 
based on spectral indices.

.. [GOR03]
   | Gorlova et al. 2003
   | **Gravity Indicators in the Near-Infrared Spectra of Brown Dwarfs**
   | *The Astrophysical Journal, Volume 593, Issue 1, pp. 1074-1092*
   | `https://arxiv.org/abs/astro-ph/0305147
     <https://arxiv.org/abs/astro-ph/0305147>`_
     
.. [SLE04]
   | Slesnick et al. 2004
   | **The Spectroscopically Determined Substellar Mass Function of the Orion 
     Nebula Cluster**
   | *The Astrophysical Journal, Volume 610, Issue 1, pp. 1045-1063*
   | `https://arxiv.org/abs/astro-ph/0404292
     <https://arxiv.org/abs/astro-ph/0404292>`_
          
.. [ALL07]
   | Allers et al. 2007
   | **Characterizing Young Brown Dwarfs Using Low-Resolution Near-Infrared 
     Spectra**
   | *The Astrophysical Journal, Volume 657, Issue 1, pp. 511-520*
   | `https://arxiv.org/abs/astro-ph/0611408
     <https://arxiv.org/abs/astro-ph/0611408>`_          
"""

__author__ = 'V. Christiaens'
__all__ = ['spectral_idx',
           'sp_idx_to_spt',
           'sp_idx_to_gravity',
           'digit_to_spt',
           'spt_to_digit']

import numpy as np
from .utils_spec import find_nearest

def spectral_idx(lbda, spec, band='H2O-1.5', spec_err=None, verbose=False):
    """
    Computes a spectral index. Implemented so far:
        - the Na 1.1 mu index [ALL07]
        - the H2O 1.3 mu index [GOR03]
        - the H2O 1.5 mu index [ALL07]
        - the H2O 2 index [SLE04]
        - the CO 2.3 index [GOR03].
        
    Parameters
    ----------
    lbda : numpy ndarray
        1d numpy array containing the wavelengths of the spectrum in microns.
    spec : numpy ndarray
        1d numpy array containing the measured flux (arbitrary units accepted).
    band: str, optional {'H2O-1.5', 'H2O-1.3', 'H2O-2', 'Na-1.1', 'CO-2.3'}
        Name of the band where the spectral index is defined (spectral feature 
        + wavelength in mu)
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
    
    if band == 'Na-1.1':
        lbda_num0 = 1.15
        lbda_num1 = 1.16
        lbda_den0 = 1.134
        lbda_den1 = 1.144
    elif band == 'H2O-1.3':
        lbda_num0 = 1.334
        lbda_num1 = 1.338
        lbda_den0 = 1.32
        lbda_den1 = 1.324    
    elif band == 'H2O-1.5':
        lbda_num0 = 1.55
        lbda_num1 = 1.56
        lbda_den0 = 1.492
        lbda_den1 = 1.502
    elif band == 'H2O-2':
        lbda_num0 = 2.035
        lbda_num1 = 2.045
        lbda_den0 = 2.145
        lbda_den1 = 2.155
    elif band == 'CO-2.3':
        lbda_num0 = 2.296
        lbda_num1 = 2.298
        lbda_den0 = 2.286
        lbda_den1 = 2.288
    else:
        raise ValueError("Band not recognised")
        
    idx_num0 = find_nearest(lbda, lbda_num0)
    idx_num1 = find_nearest(lbda, lbda_num1)
    # set floor/ceil constraints if the indices are the same
    if idx_num0 == idx_num1:
        idx_num0 = find_nearest(lbda, lbda_num0, constraint='floor')
        idx_num1 = find_nearest(lbda, lbda_num0, constraint='ceil')
    idx_den0 = find_nearest(lbda, lbda_den0)
    idx_den1 = find_nearest(lbda, lbda_den1)
    # set floor/ceil constraints if the indices are the same
    if idx_den0 == idx_den1:
        idx_den0 = find_nearest(lbda, lbda_den0, constraint='floor')
        idx_den1 = find_nearest(lbda, lbda_den1, constraint='ceil')
        
    raise_flag = False
    # Check if enough spectral resolution
    err_msg = "Indices overlap for num and den: not enough spectral resolution"
    if idx_den0 == idx_num1:
        raise_flag = True
        if abs(lbda[idx_den0]-lbda_den0)<abs(lbda[idx_num1]-lbda_num1):
            idx_num1 -= 1
            if idx_num0 == idx_num1+1:
                raise ValueError(err_msg)
        else:
            idx_den0 += 1
            if idx_den0 == idx_den1+1:
                raise ValueError(err_msg)  
    if idx_den1 == idx_num0:
        raise_flag = True
        if abs(lbda[idx_den1]-lbda_den1)<abs(lbda[idx_num0]-lbda_num0):
            idx_num0 += 1
            if idx_num0 == idx_num1+1:
                raise ValueError(err_msg)
        else:
            idx_den1 -= 1
            if idx_den0 == idx_den1+1:
                raise ValueError(err_msg)
        
    if verbose:
        if raise_flag:
            msg = "*WARNING*: spectral resolution is barely sufficient for "
            msg += "non-overlapping indices - results to be taken with caution!"
            print(msg)
        msg = "Closest (non-overlapping) wavelengths to {:.3f}, {:.3f}, {:.3f} and"
        msg+= " {:.3f} mu correspond to {:.3f}, {:.3f}, {:.3f} and {:.3f} mu resp."
        print(msg.format(lbda_num0, lbda_num1, lbda_den0, lbda_den1, 
                         lbda[idx_num0], lbda[idx_num1], lbda[idx_den0],
                         lbda[idx_den1]))
    
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


def sp_idx_to_spt(idx, name='H2O-1.5', idx_err=None, young=True):
    """
    Estimates a spectral type from a spectral index. Implemented so far:
        - the H2O 1.3 mu index [GOR03]
        - the H2O 1.5 mu index [ALL07]
        - the H2O-2 index [SLE04]
        
    Note on scale of SpT: 0 = M0, 10 = L0.
        
    Parameters
    ----------
    idx : float
        Value of spectral index
    name: str, optional {'H2O-1.3', 'H2O-1.5, 'H2O-2'}
        The name of the spectral index.
    idx_err: float, optional
        Uncertainty on the spectral index value
    young: bool, opt
        Whether the object is likely young (only used for 'H2O-1.3' index, which
        is gravity dependent)

    Returns
    -------
    spt : float
        Value of the spectral type
    spt_err: float
       [if idx_err is provided] Uncertainty on the spectral type.

    """
    
    if name == 'H2O-1.3':
        if young:
            a = 10.99
            da = 6.06
            b = 17.45
            db = 7.68
        else:
            a = 21.23
            da = 7.82
            b = 29.11
            db = 9.56
        
        spt = a - b*idx
        tot_err = da**2 # err related to a
        tot_err += np.power(db*idx,2)  # err related to b
        if idx_err is not None:
            tot_err += np.power(idx_err*b, 2)
        spt_err = np.sqrt(tot_err)

    elif name == 'H2O-1.5':
        a = 0.75
        da = 0.03
        b = 0.044
        db = 0.004
        
        spt = (idx-a)/b
        tot_err = np.power(da/b,2) # err related to a
        tot_err += np.power(db*(idx-a)/b**2,2)  # err related to b
        if idx_err is not None:
            tot_err += np.power(idx_err/b, 2)
        spt_err = np.sqrt(tot_err)

    elif name == 'H2O-2':
        a = 34.13
        da = 1.19
        b = 27.10
        db = 1.20
        spt = a - b*idx
        tot_err = da**2 # err related to a
        tot_err += np.power(db*idx,2)  # err related to b
        if idx_err is not None:
            tot_err += np.power(idx_err*b, 2)
        spt_err = np.sqrt(tot_err)

    else:
        raise ValueError("not implemented for this line")

    
    return spt, spt_err


def sp_idx_to_gravity(idx, name='Na-1.1'):
    """
    Provides a qualitative estimate of the gravity/youth based on a 
    gravity-sensitive spectral index. Implemented so far:
        - the Na-1.1 index [ALL07]
        - the CO-2.3 index [GOR03]
        
    Parameters
    ----------
    idx : float
        Value of spectral index
    name: str, optional {'Na-1.1', 'CO-2.3'}
        The name of the spectral index.

    """
    
    if name == 'Na-1.1':
        if idx[0] < 1.4:
            msg = "The object's gravity is consistent with being very young "
            msg+= "(Na index less than 1.4; Allers et al. 2007))"
        else:
            msg= "The Na index does not suggest a very low gravity (i.e. it is"
            msg+= " unlikely to be very young; Allers et al. 2007))"
    elif name == 'CO-2.3':
        if idx[0] < 0.81:
            msg = "The CO-2.3 index is consistent with either a red giant or "
            msg += " a field L-dwarf with log(g) > 4.5 (CO-2.3 index less than"
            msg += " 0.88; Gorlova et al. 2003)"
        elif idx[0] < 0.88:
            msg = "The object's gravity is likely higher than 4.5, with a SpT "
            msg += "earlier than L0 (CO-2.3 index less than 0.88; Gorlova et "
            msg += "al. 2003)"
        else:
            msg= "The object is likely young, with a log(g) lower than 4.5 "
            msg+= "(CO-2.3 index larger than 0.88; Gorlova et al. 2003)"

    print(msg)


def spt_to_digit(spt, convention='splat'):
    """
    Converts a string representing spectral type into an integer index.
    
    Parameters
    ----------
    spt : str
        String representing the spectral index
    convention: str, optional {'splat', 'Allers+07'}
        Which convention to use to convert digit into spectral type.
        Convention from splat: K0 = 0, M0=10, L0=20, T0=30, Y9 = 49
        Convention from [ALL07]: M0 = 0, L0 = 10, ...

    Returns
    -------
    idx : float or int
        Index value of the spectral type
    """

    spt_str = spt.split()[0]
    sub = 0
    if convention == 'Allers+07':
        sub=10
        
    if 'K' in spt_str:
        idx = 0-sub
        loc = spt_str.find('K')
        if convention == 'Allers+07':
            msg="WARNING: Allers+07 convention doesn't handle SpT earlier than M"
            print(msg)        
    elif 'M' in spt_str:
        idx = 10-sub
        loc = spt_str.find('M')
    elif 'L' in spt_str:
        idx = 20-sub
        loc = spt_str.find('L')
    elif 'T' in spt_str:
        idx = 30-sub
        loc = spt_str.find('T')
    elif 'Y' in spt_str:
        idx = 40-sub
        loc = spt_str.find('Y')
    else:
        idx = np.nan
        print("WARNING: Spectral type not between K0 and Y9. Idx set to NaN.")

    if len(spt)>1:
        if len(spt[loc+1:]) > 2:
            if spt[loc+2] == ".":
                idx+=float(spt[loc+1:loc+4])
        elif spt[loc+1].isdigit():
            idx+=float(spt[loc+1])
    else:
        idx+=5 # in case just "M" or "L" is given, let's assume it's M5 or L5
    
    return idx

    
def digit_to_spt(idx, convention='splat'):
    """
    Converts an integer index into spectral type.

    
    Parameters
    ----------
    idx : float or int
        Index value of the spectral type 
    convention: str, optional {'splat', 'Allers+07'}
        Which convention to use to convert digit into spectral type.
        Convention from splat: K0 = 0, M0=10, L0=20, T0=30, Y9 = 49
        Convention from Allers+07: M0 = 0, L0 = 10, ...

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
    