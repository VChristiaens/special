#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module with configuration parameters, timing functions and multiprocessing 
utilities (inspired from VIP).
"""

__author__ = 'Valentin Christiaens, C. A. Gomez Gonzalez'

__all__ = ['time_ini',
           'timing',
           'time_fin',
           'sep',
           'figsize',
           'figdpi']


from datetime import datetime
import itertools as itt
import multiprocessing
import sys

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


class Progressbar(object):
    """ Show progress bars. Supports multiple backends.

    Examples
    --------
    .. code:: python

        from vip_hci.var import Progressbar
        Progressbar.backend = "tqdm"

        from time import sleep

        for i in Progressbar(range(50)):
            sleep(0.02)

        # or:

        bar = Progressbar(total=50):
        for i in range(50):
            sleep(0.02)
            bar.update()

        # Progressbar can be disabled globally using
        Progressbar.backend = "hide"

        # or locally using the ``verbose`` keyword:
        Progressbar(iterable, verbose=False)

    Notes
    -----
    - `leave` keyword is natively supported by tqdm, support could be added to
      other backends too?

    """
    backend = "pyprind"

    def __new__(cls, iterable=None, desc=None, total=None, leave=True,
                backend=None, verbose=True):
        if backend is None:
            backend = Progressbar.backend

        if not verbose:
            backend = "hide"

        if backend == "tqdm":
            from tqdm import tqdm
            return tqdm(iterable=iterable, desc=desc, total=total, leave=leave,
                        ascii=True, ncols=80, file=sys.stdout,
                        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed"
                                   "}<{remaining}{postfix}]") # remove rate_fmt
        elif backend == "tqdm_notebook":
            from tqdm import tqdm_notebook
            return tqdm_notebook(iterable=iterable, desc=desc, total=total,
                                 leave=leave)
        elif backend == "pyprind":
            from pyprind import ProgBar, prog_bar
            ProgBar._adjust_width = lambda self: None  # keep constant width
            if iterable is None:
                return ProgBar(total, title=desc, stream=1)
            else:
                return prog_bar(iterable, title=desc, stream=1,
                                iterations=total)
        elif backend == "hide":
            return NoProgressbar(iterable=iterable)
        else:
            raise NotImplementedError("unknown backend")

    def set(b):
        Progressbar.backend = b
        
        
class NoProgressbar(object):
    """ Wraps an ``iterable`` to behave like ``Progressbar``, but without
    producing output.
    """
    def __init__(self, iterable=None):
        self.iterable = iterable

    def __iter__(self):
        return self.iterable.__iter__()

    def __next__(self):
        return self.iterable.__next__()

    def __getattr__(self, key):
        return self.iterable.key

    def update(self):
        pass
    
    
def eval_func_tuple(f_args):
    """ Takes a tuple of a function and args, evaluates and returns result"""
    return f_args[0](*f_args[1:])


class FixedObj(object):
    def __init__(self, v):
        self.v = v


def iterable(v):
    """ Helper function for ``pool_map``: prevents the argument from being
    wrapped in ``itertools.repeat()``.

    Examples
    --------
    .. code-block:: python

        # we have a worker function whic processes a word:

        def worker(word, method):
            # ...

        # we want to process these words in parallel fasion:
        words = ["lorem", "ipsum", "esse", "ea", "eiusmod"]

        # but all with
        method = 1

        # we then would use
        pool_map(3, worker, iterable(words), method)

        # this results in calling
        #
        # worker(words[0], 1)
        # worker(words[1], 1)
        # worker(words[2], 1)
        # ...
    """
    return FixedObj(v)


def pool_map(nproc, fkt, *args, **kwargs):
    """
    Abstraction layer for multiprocessing. When ``nproc=1``, the builtin
    ``map()`` is used. For ``nproc>1`` a ``multiprocessing.Pool`` is created.

    Parameters
    ----------
    nproc : int
        Number of processes to use.
    fkt : callable
        The function to be called with each ``*args``
    *args : function arguments
        Arguments passed to ``fkt`` By default, ``itertools.repeat`` is applied
        on all the arguments, except when you wrap the argument in
        ``iterable()``.
    msg : str or None, optional
        Description to be displayed.
    progressbar_single : bool, optional
        Display a progress bar when single-processing is used. Defaults to
        ``False``.
    verbose : bool, optional
        Show more output. Also disables the progress bar when set to ``False``.

    Returns
    -------
    res : list
        A list with the results.

    """
    multiprocessing.set_start_method('fork', force=True)
    from multiprocessing import Pool
    
    msg = kwargs.get("msg", None)
    verbose = kwargs.get("verbose", True)
    progressbar_single = kwargs.get("progressbar_single", False)
    _generator = kwargs.get("_generator", False)  # not exposed in docstring

    args_r = [a.v if isinstance(a, FixedObj) else itt.repeat(a) for a in args]
    z = zip(itt.repeat(fkt), *args_r)

    if nproc == 1:
        if progressbar_single:
            total = len([a.v for a in args if isinstance(a, FixedObj)][0])
            z = Progressbar(z, desc=msg, verbose=verbose, total=total)
        res = map(eval_func_tuple, z)
        if not _generator:
            res = list(res)
    else:
        if verbose and msg is not None:
            print("{} with {} processes".format(msg, nproc))
        pool = Pool(processes=nproc)
        if _generator:
            res = pool.imap(eval_func_tuple, z)
        else:
            res = pool.map(eval_func_tuple, z)
        pool.close()
        pool.join()

    return res


def pool_imap(nproc, fkt, *args, **kwargs):
    """
    Generator version of ``pool_map``. Useful when showing a progress bar for
    multiprocessing (see examples).

    Parameters
    ----------
    nproc : int
        Number of processes to use.
    fkt : callable
        The function to be called with each ``*args``
    *args : function arguments
        Arguments passed to ``fkt``
    msg : str or None, optional
        Description to be displayed.
    progressbar_single : bool, optional
        Display a progress bar when single-processing is used. Defaults to
        ``True``.
    verbose : bool, optional
        Show more output. Also disables the progress bar when set to ``False``.

    Examples
    --------

    .. code-block:: python

        # using pool_map

        res = pool_map(2, my_worker_function, *args)

        # using pool_imap with a progessbar:

        res = list(Progressbar(pool_imap(2, my_worker_function, *args)))

    """
    kwargs["_generator"] = True
    return pool_map(nproc, fkt, *args, **kwargs)


def repeat(*args):
    """
    Applies ``itertools.repeat`` to every ``args``.


    Examples
    --------

    # instead of using

    import itertools as itt

    my_fkt(itt.repeat(a), itt.repeat(b), itt.repeat(c), d, itt.repeat(e))

    # you could use `repeat`:

    my_fkt(*repeat(a, b, c), d, *repeat(e))

    """
    return [itt.repeat(a) for a in args]
