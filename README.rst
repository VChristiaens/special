spec_fit
=======

.. image:: https://badge.fury.io/py/spec_fit.svg
    :target: https://pypi.python.org/pypi/spec_fit

.. image:: https://img.shields.io/badge/Python-3.5%2C%203.6%2C%203.7%2C%203.8%2C%203.9-brightgreen.svg
    :target: https://pypi.python.org/pypi/spec_fit

.. image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat
    :target: https://github.com/VChristiaens/spec_fit/blob/master/LICENSE


``spec_fit`` is a package for the spectral characterization of (sub-)stellar objects.

This package provides the tools for the analysis of measured spectra (e.g. low-res spectra of directly imaged companions, but not exclusively). It enables:

- fitting of input spectra to different grids of models (the code is grid-independent, download of the grid is left to the user);  
- MCMC sampling of the model grid parameter space;
- searching for the best-fit template spectrum within a given library.

In particular, the MCMC sampler routine makes use of the `emcee` package (Foreman-Mackey et al. 2013):

- is flexible, as it is usable on any grid of models downloaded by the user (provided a snippet function specifying the format of the input);
- can fit (additional) blackbody components;
- can sample the effect of extinction; 
- can sample different extinction laws than ISM (parametrised using RV);
- accepts either uniform or Gaussian priors for each model parameter;
- accepts a prior on the mass of the object (if surface gravity is one of the model parameters);
- considers convolution with the spectral PSF and/or resampling of the model for consistency with the input spectrum.

Furthermore, the log-likelihood expression has the option to include:

- spectral correlation between measurements of adjacent channels of a given instrument;
- weights that are proportional to the relative spectral span of each measurement, in case these are obtained from different instruments.

More details are available in `Christiaens et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021MNRAS.502.6117C/abstract>`_.
Please cite this publication if you use `specfit` for your research, along with `Foreman-Mackey et al. (2013) <https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract>`_ if you use the MCMC sampler.


Documentation
-------------
The documentation for ``spec_fit`` can be found here: TBD


Jupyter notebook tutorial
-------------------------
TBD


TL;DR setup guide
-----------------
.. code-block:: bash

    $ pip install spec_fit


Installation and dependencies
-----------------------------
The benefits of using a Python package manager (distribution), such as
(ana)conda or Canopy, are many. Mainly, it brings easy and robust package
management and avoids messing up with your system's default python. An
alternative is to use package managers like apt-get for Ubuntu or
Homebrew/MacPorts/Fink for macOS. We recommend using 
`Miniconda <https://conda.io/miniconda>`_.

``spec_fit`` depends on existing packages from the Python ecosystem, such as
``numpy``, ``scipy``, ``matplotlib``, ``pandas`` and ``astropy``. There are different ways of
installing ``spec_fit`` suitable for different scenarios.


Using pip
^^^^^^^^^
The easiest way to install ``spec_fit`` is through the Python Package Index, aka
`PyPI <https://pypi.org/>`_, with the ``pip`` package manager. Simply run:

.. code-block:: bash

  $ pip install spec_fit

With ``pip`` you can easily uninstall, upgrade or install a specific version of
``spec_fit ``. For upgrading the package run:

.. code-block:: bash

  $ pip install --upgrade spec_fit

Alternatively, you can use ``pip install`` and point to the GitHub repo:

.. code-block:: bash

  $ pip install git+https://github.com/VChristiaens/spec_fit.git

Using the setup.py file
^^^^^^^^^^^^^^^^^^^^^^^
You can download ``spec_fit`` from its GitHub repository as a zip file. A ``setup.py``
file (setuptools) is included in the root folder of ``spec_fit``. Enter the package's
root folder and run:

.. code-block:: bash

  $ python setup.py install


Using Git
^^^^^^^^^
If you plan to contribute or experiment with the code you need to make a 
fork of the repository (click on the fork button in the top right corner) and 
clone it:

.. code-block:: bash

  $ git clone https://github.com/<replace-by-your-username>/spec_fit.git

If you do not create a fork, you can still benefit from the ``git`` syncing
functionalities by cloning the repository (but will not be able to contribute):

.. code-block:: bash

  $ git clone https://github.com/VChristiaens/spec_fit.git

Before installing the package, it is highly recommended to create a dedicated
conda environment to not mess up with the package versions in your base 
environment. This can be done easily with (replace spec_env by the name you want
for your environment):

.. code-block:: bash

  $ conda create -n spec_env python=3.9 ipython

Note: installing ipython while creating the environment with the above line will
avoid a commonly reported issue which stems from trying to import VIP from 
within a base python2.7 ipython console.

To install spec_fit, simply cd into the spec_fit directory and run the setup file 
in 'develop' mode:

.. code-block:: bash

  $ cd VIP
  $ python setup.py develop

If cloned from your fork, make sure to link your spec_fit directory to the upstream 
source, to be able to easily update your local copy when a new version comes 
out or a bug is fixed:

.. code-block:: bash

  $ git add remote upstream https://github.com/VChristiaenss/spec_fit.git


Loading specfit
^^^^^^^^^^^^^^^
Finally, start Python or IPython and check
that you are able to import ``spec_fit``:

.. code-block:: python

  import spec_fit as specfit

If everything went fine with the installation, you will see a welcome message.
Now you can start characterizing exoplanets and other (sub)stellar objects!
