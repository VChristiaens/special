.. image:: docs/_static/Special_logo.png
   :align: center
   :width: 825px


.. image:: https://badge.fury.io/py/special.svg?branch=main&service=github
    :target: https://badge.fury.io/py/special

.. image:: https://img.shields.io/badge/Python-3.6%2C%203.7%2C%203.8%2C%203.9-brightgreen.svg
    :target: https://pypi.python.org/pypi/special.svg

.. image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat
    :target: https://github.com/VChristiaens/special/blob/master/LICENSE

.. image:: https://joss.theoj.org/papers/10.21105/joss.04456/status.svg
   :target: https://doi.org/10.21105/joss.04456

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.7181950.svg
   :target: https://doi.org/10.5281/zenodo.7181950

``special`` is a Python package for the SPEctral Characterization of directly ImAged Low-mass companions. While some tools are specific to the characterisation of low-mass (M, L, T) dwarfs down to giant planets at optical/IR wavelengths, the main routines of ``special`` (MCMC and nested samplers) can also be used in a more general way for the characterisation of any type of object with a measured spectrum, provided a relevant input model grid, regardless of the observational method used to obtain the spectrum (direct imaging or not) and regardless of the format of the spectra (multi-band photometry, low-resolution or medium-resolution spectrum, or a combination thereof).

This package provides the following tools for the analysis of measured spectra:

* calculation of the spectral correlation between channels of an IFS datacube (relevant to directly imaged companions with an IFS, where the uncertainty reflects spectrally correlated residual speckle noise);
* calculation of empirical spectral indices for MLT-dwarfs;
* fitting of input spectra to either photo-/atmospheric model grids or a blackbody model, including additional parameters such as (extra) black body component(s), extinction and total-to-selective extinction ratio;
* using either MCMC (`emcee <https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract>`_) or nested (`nestle <http://github.com/kbarbary/nestle>`_ or `UltraNest <https://johannesbuchner.github.io/UltraNest/>`_) samplers to infer posterior distributions on spectral model parameters in a Bayesian framework;
* searching for the best-fit template spectrum within a given template library, with up to two free parameters (relative flux scaling and extinction).


The MCMC and nested sampler routines have been adapted to:

* be flexible, as they are usable on any grid of models provided by the user (along with a snippet function specifying how to read the format of the input files);
* sample the effect of (additional) blackbody components;
* sample the effect of extinction (AV); 
* sample different extinction laws than ISM (parametrised using the total-to-selective extinction ratio RV);
* sample a list of potential emission lines;
* accept either uniform or Gaussian priors for each model parameter;
* accept a prior on the mass of the object (if surface gravity is one of the model parameters, and for the MCMC sampler only);
* consider convolution with the line spread function, photometric filters transmission and/or resampling of the model for consistency with the input spectrum - in particular convolution and resampling are done in two consecutive steps, and multiple resolving powers can be provided as input;
* use a log-likelihood expression that can include i) spectral correlation between measurements of adjacent channels of a given instrument, and ii) additional weights that are proportional to the relative spectral bandwidth of each measurement, in case these are obtained from different instruments (e.g. photometry+spectroscopy).


Documentation
-------------
The documentation for ``special`` can be found `here <https://special.readthedocs.io/en/latest/>`_.
``special`` was originally implemented as ``specfit``, a former module of the ``VIP`` package, before undergoing significant expansion. It was first presented in `Christiaens et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021MNRAS.502.6117C/abstract>`_ . More details will be available in an upcoming publication (Christiaens et al., subm. to JOSS).


Jupyter notebook tutorial
-------------------------
A Jupyter notebook tutorial examplifying most possibilities within ``special`` is available in the 
``special-extras`` `repository <https://github.com/VChristiaens/special_extras>`_. 
Alternatively, you can execute this repository on 
`Binder <https://mybinder.org/v2/gh/VChristiaens/special_extras/main>`_ (in the tutorials directory), or go through it in `the documentation <https://special.readthedocs.io/en/latest/tutorials/walkthrough.html>`_.


TL;DR setup guide
-----------------
.. code-block:: bash

    $ pip install special


Installation and dependencies
-----------------------------
The benefits of using a Python package manager (distribution), such as
(ana)conda or Canopy, are many. Mainly, it brings easy and robust package
management and avoids messing up with your system's default python. An
alternative is to use package managers like apt-get for Ubuntu or
Homebrew/MacPorts/Fink for macOS. We recommend using 
`Miniconda <https://conda.io/miniconda>`_.

``special`` depends on existing packages from the Python ecosystem, such as
``numpy``, ``scipy``, ``matplotlib``, ``pandas`` and ``astropy``. There are different ways of
installing ``special`` suitable for different scenarios.


Using pip
^^^^^^^^^
The easiest way to install ``special`` is through the Python Package Index, aka
`PyPI <https://pypi.org/>`_, with the ``pip`` package manager. Simply run:

.. code-block:: bash

  $ pip install special

With ``pip`` you can easily uninstall, upgrade or install a specific version of
``special``. For upgrading the package run:

.. code-block:: bash

  $ pip install --upgrade special

Alternatively, you can use ``pip install`` and point to the GitHub repo:

.. code-block:: bash

  $ pip install git+https://github.com/VChristiaens/special.git

Using the setup.py file
^^^^^^^^^^^^^^^^^^^^^^^
You can download ``special`` from its GitHub repository as a zip file. A ``setup.py``
file (setuptools) is included in the root folder of ``special``. Enter the package's
root folder and run:

.. code-block:: bash

  $ python setup.py install


Using Git
^^^^^^^^^
If you plan to contribute or experiment with the code you need to make a 
fork of the repository (click on the fork button in the top right corner) and 
clone it:

.. code-block:: bash

  $ git clone https://github.com/<replace-by-your-username>/special.git

If you do not create a fork, you can still benefit from the ``git`` syncing
functionalities by cloning the repository (but will not be able to contribute):

.. code-block:: bash

  $ git clone https://github.com/VChristiaens/special.git

Before installing the package, it is highly recommended to create a dedicated
conda environment to not mess up with the package versions in your base 
environment. This can be done easily with (replace spec_env by the name you want
for your environment):

.. code-block:: bash

  $ conda create -n spec_env python=3.9 ipython

Note: installing ipython while creating the environment with the above line will
avoid a commonly reported issue which stems from trying to import ``special`` from 
within a base python2.7 ipython console.

To install ``special``, simply cd into the special directory and run the setup file 
in 'develop' mode:

.. code-block:: bash

  $ cd special
  $ python setup.py develop

If cloned from your fork, make sure to link your special directory to the upstream 
source, to be able to easily update your local copy when a new version comes 
out or a bug is fixed:

.. code-block:: bash

  $ git add remote upstream https://github.com/VChristiaens/special.git


Loading `special`
^^^^^^^^^^^^^^^^^
Finally, start Python or IPython and check that you are able to import ``special``:

.. code-block:: python

  import special

Now you can start characterizing exoplanets and other (sub)stellar objects!



About `special`
---------------

Contributions
^^^^^^^^^^^^^
Feel free to fork the repository and submit a pull request with either new features or bug fixes. External contributions are very welcome. In particular, please check the expected future `areas for development <https://github.com/VChristiaens/special/projects/1>`_.


Questions and suggestions
^^^^^^^^^^^^^^^^^^^^^^^^^
``special`` was developed by Valentin Christiaens. Feel free to contact me at valentin.christiaens@uliege.be if you have any question or suggestion.


Acknowledgements
^^^^^^^^^^^^^^^^
Please cite `Christiaens et al. (2022) <https://joss.theoj.org/papers/10.21105/joss.04456>`_ if you use ``special`` for your research, along with (where relevant):

- `Foreman-Mackey et al. (2013) <https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract>`_ if you use the ``emcee`` MCMC sampler;
- `Skilling (2004) <https://ui.adsabs.harvard.edu/abs/2004AIPC..735..395S/abstract>`_, `Mukherjee et al. (2006) <https://ui.adsabs.harvard.edu/abs/2006ApJ...638L..51M/abstract>`_, or `Feroz et al. (2009) <https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.1601F/abstract>`_ if you use the nested sampler `nestle` in 'classic', 'single' or 'multi' mode, respectively. Please also mention the ``nestle`` `GitHub repository <http://github.com/kbarbary/nestle>`_;
- `Buchner (2021) <https://ui.adsabs.harvard.edu/abs/2021JOSS....6.3001B/abstract>`_ if you use the `UltraNest` nested sampler.
