#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module with configuration parameters and timing functions (inspired from VIP).
"""

__author__ = 'Valentin Christiaens'

__all__ = ['time_ini',
           'timing',
           'time_fin',
           'sep',
           'figsize',
           'figdpi']

from datetime import datetime

sep = '―' * 80
figsize = (8, 5)
figdpi = 100


def time_ini(verbose=True):
    """
    Set and print the time at which the script started.

    Returns
    -------
    start_time : string
        Starting time.

    """
    start_time = datetime.now()
    if verbose:
        print(sep)
        print("Starting time: " + start_time.strftime("%Y-%m-%d %H:%M:%S"))
        print(sep)
    return start_time


def timing(start_time):
    """
    Print the execution time of a script.

    It requires the initialization  with the function time_ini().
    """
    print("Running time:  " + str(datetime.now()-start_time))
    print(sep)


def time_fin(start_time):
    """
    Return the execution time of a script.

    It requires the initialization  with the function time_ini().
    """
    return str(datetime.now()-start_time)

