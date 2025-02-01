r"""For specifying HRTEM systems, HRTEM simulation parameters, and running HRTEM
simulations.

The algorithm used to perform HRTEM simulations is similar to the PRISM
algorithm used to perform STEM simulations. See the documentation for the
subpackage :mod:`prismatique.stem` for a discussion on the PRISM algorithm.

The PRISM algorithm can be used to perform STEM simulations as an alternative to
the multislice algorithm. Here we discuss very briefly aspects of both the PRISM
and multislice algorithms, although it is recommended that the user see
Ref. [Ophus1]_ for more details on the former and Ref. [Kirkland1]_ for an
exposition on the latter. 

As discussed in the documentation for the class
:class:`prismatique.thermal.Params`, the intensity pattern for a given beam in a
HRTEM experiment depends on the state operator :math:`\hat{\rho}_{t}` of a
transmitted beam electron, which is assumed to be a weighted sum of pure states

.. math ::
    &\hat{\rho}_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
    \boldsymbol{\delta}_{\beta}\right)\\&\quad=\left|\psi_{t}\left(\delta_{f};
    \mathbf{u}_{1},\ldots,\mathbf{u}_{N};
    \boldsymbol{\delta}_{\beta}\right)\right\rangle 
    \left\langle \psi_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
    \boldsymbol{\delta}_{\beta}\right)\right|,
    :label: pure_state_operator_for_transmitted_electron_3

where
:math:`\left|\psi_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
\boldsymbol{\delta}_{\beta}\right)\right\rangle` being the state vector of a
transmitted beam electron for a perfectly coherent beam and a sample in a frozen
atomic configuration :math:`\left\{ \mathbf{u}_{j}\right\} _{j=1}^{N}`, with the
:math:`\delta_{f}` and :math:`\boldsymbol{\delta}_{\beta}` implicitly specifying
the defocus and beam tilt respectively. Note that the weighted sum is over the
defocus, beam tilt, and :math:`\left\{ \mathbf{u}_{j}\right\} _{j=1}^{N}`. 

Each :math:`\left|\psi_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
\boldsymbol{\delta}_{\beta}\right)\right\rangle` is calculated by

.. math ::
    \left|\psi_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
    \boldsymbol{\delta}_{\beta}\right)\right\rangle =
    \hat{Q}\left\{ S_{m_{x},m_{y}}\right\} \left(x,y\right),
    :label: calculating_psi_t_in_hrtem

.. math ::
    \hat{Q}\left\{ g\right\} \left(x,y\right)=
    \left(\hat{A}\circ\hat{D}\right)\left\{ g\right\} \left(x,y\right),
    :label: Q_operator

where the notation :math:`\hat{O}\left\{ g\right\} \left(q_{x},q_{y}\right)`
denotes an operator :math:`\hat{O}` transforming a function :math:`g` into a
space with coordinates :math:`q_{x}` and :math:`q_{y}` [e.g. :math:`\hat{O}`
could be a two-dimensional Fourier transform]; the notation
:math:`\left(\hat{O}_{1}\circ\hat{O}_{2}\right)\left\{ g\right\}
\left(q_{x},q_{y}\right)` denotes the composition of two operators
:math:`\hat{O}_{1}` and :math:`\hat{O}_{2}` transforming a function :math:`g`
into a space with coordinates :math:`q_{1}` and :math:`q_{2}` with
:math:`\hat{O}_{2}` being applied first; :math:`\hat{A}` is given by
Eq. :eq:`A_operator`;

.. math ::
    \hat{D}\left\{ g\right\} \left(k_{x},k_{y}\right)\equiv
    \hat{\mathcal{F}}_{\text{2D}}\left\{ g\right\} \left(k_{x},k_{y}\right)
    e^{-i\chi\left(k_{x},k_{y}\right)}\xi\left(k_{x},k_{y}\right),
    :label: D_operator

with :math:`\hat{\mathcal{F}}_{\text{2D}}` being given by
Eq. :eq:`fourier_transform_2D`, :math:`\xi\left(k_{x},k_{y}\right)` being an
aperture function modelling the objective aperture, and
:math:`\chi\left(k_{x},k_{y}\right)` being the phase deviation due to coherent
lens aberrations of the objective lens; and :math:`S_{m_{x},m_{y}}` is the
element of the :math:`S`-matrix corresponding to a perfectly coherent plane wave
incident at the sample with the beam tilt implicitly specified by
:math:`\boldsymbol{\delta}_{\beta}`. See the documentation for the subpackage
:mod:`prismatique.stem` for a discussion on :math:`S`-matrices and the
documentation for the class :class:`embeam.coherent.Aberration` for a definition
of :math:`\chi\left(k_{x},k_{y}\right)`. Note that
:math:`\xi\left(k_{x},k_{y}\right)=0` for all scattering directions
:math:`\left(k_{x},k_{y}\right)` blocked by the objective aperture, and
:math:`\xi\left(k_{x},k_{y}\right)=1` otherwise.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# Import child modules and packages of current package.
import prismatique.hrtem.system
import prismatique.hrtem.image
import prismatique.hrtem.output
import prismatique.hrtem.sim



############################
## Authorship information ##
############################

__author__     = "Matthew Fitzpatrick"
__copyright__  = "Copyright 2023"
__credits__    = ["Matthew Fitzpatrick"]
__maintainer__ = "Matthew Fitzpatrick"
__email__      = "mrfitzpa@uvic.ca"
__status__     = "Development"



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = [""]



###########################
## Define error messages ##
###########################