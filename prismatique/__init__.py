"""``prismatique`` is a Python library that functions essentially as a wrapper
to the Python library ``pyprismatic``, which itself is a thin wrapper to
``prismatic``, a CUDA/C++ package for fast image simulations in scanning
transmission electron microscopy and high-resolution transmission electron
microscopy. You can find more information about ``pyprismatic`` and
``prismatic`` `here <https://prism-em.com/>`_.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# Import child modules and packages of current package.
# import prismatique.sim
import prismatique.worker
import prismatique.thermal
import prismatique.discretization
import prismatique.sample
import prismatique.scan
import prismatique.tilt
import prismatique.aperture
import prismatique.cbed
import prismatique.stem
import prismatique.hrtem
import prismatique.load
import prismatique.version



############################
## Authorship information ##
############################

__author__       = "Matthew Fitzpatrick"
__copyright__    = "Copyright 2023"
__credits__      = ["Matthew Fitzpatrick"]
__version__      = prismatique.version.version
__full_version__ = prismatique.version.full_version
__maintainer__   = "Matthew Fitzpatrick"
__email__        = "mrfitzpa@uvic.ca"
__status__       = "Development"



###################################
## Useful background information ##
###################################

# See e.g. ``https://docs.python.org/3/reference/import.html#regular-packages``
# for a brief discussion of ``__init__.py`` files.



##################################
## Define classes and functions ##
##################################

# List of public objects in package.
__all__ = ["show_config"]



def show_config():
    r"""Print information about the version of ``prismatique`` and libraries it 
    uses.

    """
    print(version.version_summary)

    return None