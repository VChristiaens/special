"""
``special`` has helping functions for the analysis of (low-res) 
spectra, including:

- fitting of input spectra to models and templates;
- mcmc sampling of model parameter space;
- nested sampling of model parameter space;
- best fit search within a template library;
- utility functions for the spectral fit.
"""

__version__ = "0.1.2"

from .config import *
from .utils_mcmc import *
from .utils_spec import *
from .fits import *
from .spec_corr import *
from .spec_indices import *
from .chi import *
from .model_resampling import *
from .mcmc_sampling import *
from .nested_sampling import *
from .template_fit import *
