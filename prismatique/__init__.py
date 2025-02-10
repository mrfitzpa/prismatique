# -*- coding: utf-8 -*-
# Copyright 2024 Matthew Fitzpatrick.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/gpl-3.0.html>.
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

# Get version of current package.
from prismatique.version import __version__



##################################
## Define classes and functions ##
##################################

# List of public objects in package.
__all__ = []



###########################
## Define error messages ##
###########################
