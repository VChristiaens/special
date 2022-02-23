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



