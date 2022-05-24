#! /usr/bin/env python

"""
Module with utility functions to the nested sampling for parameter estimation.
"""

__author__ = 'V. Christiaens'
__all__ = ['un_burning']

import numpy as np

def un_burning(res, logger=None):
    """
    Automatic burning of UltraNest chain based on cumulated sum of weights 
    (as implemented in UltraNest's cornerplot).
     
    Note: this function is necessary to be able to make corner plots showing
    units after best estimates, as ultranest's cornerplots does not feature 
    that option and does burning+corner plot together.
    
    Parameters
    ----------
    res: UltraNest result object
        The UltraNest result.
        
    Returns
    -------
    burned_res: tuple of 2 numpy nd array
        The burned UltraNest chain and associated weights
    """
    
    paramnames = res['paramnames']
    data = np.array(res['weighted_samples']['points'])
    weights = np.array(res['weighted_samples']['weights'])
    cumsumweights = np.cumsum(weights)

    mask = cumsumweights > 1e-4

    if mask.sum() == 1:
        if logger is not None:
            warn = 'Posterior is still concentrated in a single point:'
            for i, p in enumerate(paramnames):
                v = res['samples'][mask,i]
                warn += "\n" + '    %-20s: %s' % (p, v)

            logger.warning(warn)
            logger.info('Try running longer.')
        return
    
    burned_res = (data[mask,:], weights[mask])
    
    return burned_res