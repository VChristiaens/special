.. image:: docs/_static/Special_logo.png
   :align: center
   :width: 825px


.. image:: https://badge.fury.io/py/special.svg?branch=main&service=github
    :target: https://badge.fury.io/py/special

.. image:: https://img.shields.io/badge/Python-3.6%2C%203.7%2C%203.8%2C%203.9-brightgreen.svg
    :target: https://pypi.python.org/pypi/special.svg

.. image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat
    :target: https://github.com/VChristiaens/special/blob/master/LICENSE


``special`` is a package for the SPEctral Characterization of directly ImAged Low-mass companions. While some tools are specific to the characterisation of low-mass (M, L, T) dwarfs down to giant planets, ``special`` can also be used in a more general way for the characterisation of any object with a measured spectrum, provided an input model grid.

This package provides the following tools for the analysis of measured spectra:

* calculation of the spectral correlation between channels of an IFS datacube;
* calculation of empirical spectral indices for MLT-dwarfs;
* fitting of input spectra to different grids of models, including additional parameters such as (extra) black body component(s) and extinction;
* using either the MCMC (`emcee <https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract>`_) or nested (`nestle <http://github.com/kbarbary/nestle>`_) sampler to infer posterior distributions on spectral model parameters;
* searching for the best-fit template spectrum within a given template library, with up to two free parameters (flux scaling and relative extinction).

More details are available in `Christiaens et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021MNRAS.502.6117C/abstract>`_ (note it was originally implemented as ``specfit``, a former module of the ``VIP`` package).


Documentation
-------------
The documentation for ``special`` can be found `here <https://special.readthedocs.io/en/latest/>`_.


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
Please cite `Christiaens et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021MNRAS.502.6117C/abstract>`_ if you use ``special`` for your research, along with:

- `Foreman-Mackey et al. (2013) <https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract>`_ if you use the ``emcee`` MCMC sampler;
- `Skilling (2004) <https://ui.adsabs.harvard.edu/abs/2004AIPC..735..395S/abstract>`_, `Mukherjee et al. (2006) <https://ui.adsabs.harvard.edu/abs/2006ApJ...638L..51M/abstract>`_, or `Feroz et al. (2009) <https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.1601F/abstract>`_ if you use the nested sampler in 'classic', 'single' or 'multi' mode, respectively. Please also mention the ``nestle`` `GitHub repository <http://github.com/kbarbary/nestle>`_.
