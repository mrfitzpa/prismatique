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
r"""For specifying STEM systems, STEM simulation parameters, and running STEM
simulations.

The PRISM algorithm can be used to perform STEM simulations as an alternative to
the multislice algorithm. Here we discuss very briefly aspects of both the PRISM
and multislice algorithms, although it is recommended that the user see
Ref. [Ophus1]_ for more details on the former and Ref. [Kirkland1]_ for an
exposition on the latter. 

As discussed in the documentation for the class
:class:`prismatique.thermal.Params`, the intensity pattern for a given probe in
a STEM experiment or a given beam in a HRTEM experiment depends on the state
operator :math:`\hat{\rho}_{t}` of a transmitted beam electron, which is assumed
to be a weighted sum of pure states

.. math ::
    &\hat{\rho}_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
    \boldsymbol{\delta}_{\beta}\right)\\&\quad=\left|\psi_{t}\left(\delta_{f};
    \mathbf{u}_{1},\ldots,\mathbf{u}_{N};
    \boldsymbol{\delta}_{\beta}\right)\right\rangle 
    \left\langle \psi_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
    \boldsymbol{\delta}_{\beta}\right)\right|,
    :label: pure_state_operator_for_transmitted_electron_2

where
:math:`\left|\psi_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
\boldsymbol{\delta}_{\beta}\right)\right\rangle` being the state vector of a
transmitted beam electron for a perfectly coherent beam and a sample in a frozen
atomic configuration :math:`\left\{ \mathbf{u}_{j}\right\} _{j=1}^{N}`, with the
:math:`\delta_{f}` and :math:`\boldsymbol{\delta}_{\beta}` implicitly specifying
the defocus and beam tilt respectively. Note that the weighted sum is over the
defocus, beam tilt, and :math:`\left\{ \mathbf{u}_{j}\right\} _{j=1}^{N}`. Each
:math:`\left|\psi_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
\boldsymbol{\delta}_{\beta}\right)\right\rangle` is assumed to have been evolved
from an initial pure state
:math:`\left|\psi_{i}\left(\delta_{f};\boldsymbol{\delta}_{\beta}\right)
\right\rangle` that is determined by the conditions of the corresponding
perfectly coherent beam at incidence with the sample. Each
:math:`\left|\psi_{i}\left(\delta_{f};\boldsymbol{\delta}_{\beta}\right)
\right\rangle` is evolved using the fast electron approximation. In this
approximation, the full beam electron wavefunction :math:`\Psi\left(\mathbf{r}
\right)`, for a given perfectly coherent beam and frozen atomic configuration,
is assumed to be of the following factored form [Kirkland1]_:

.. math ::
    \Psi\left(\mathbf{r}\right)=\psi\left(\mathbf{r}\right)e^{2\pi ikz},
    :label: beam_electron_wavefunction

where

.. math ::
    \mathbf{r}=x\hat{\mathbf{x}}+y\hat{\mathbf{y}}+z\hat{\mathbf{z}};
    :label: r_vector_in_prism_params

:math:`k=\left|\mathbf{k}\right|=1/\lambda` is the magnitude of the electron's
total wavevector with :math:`\lambda` being the beam electron wavelength;
:math:`e^{2\pi ikz}` is a plane wave travelling along the :math:`z`-axis
[i.e. the optic axis]; :math:`\psi\left(\mathbf{r}\right)` satisfies

.. math ::
    \frac{\partial\psi\left(\mathbf{r}\right)}{\partial z}= 
    \left[\frac{i\lambda}{4\pi}
    \left\{ \frac{\partial^{2}}{\partial x^{2}}
    +\frac{\partial^{2}}{\partial y^{2}}\right\} 
    +i\sigma V\left(\mathbf{r}\right)\right]\psi\left(\mathbf{r}\right),
    :label: fast_electron_equation

with

.. math ::
    \sigma=\frac{2\pi m_{e}e\lambda}{h^{2}},
    :label: sigma

:math:`\hbar` being Planck's constant, :math:`e` being the elementary
charge, :math:`m_{e}` being the rest mass of an electron, and
:math:`V\left(\mathbf{r}\right)` being the total potential due to the
sample, which is treated classically; and

.. math ::
    \psi\left(x\hat{\mathbf{x}}+y\hat{\mathbf{y}}\right)
    =\left\langle x,y\left|\psi_{i}\left(\delta_{f};
    \boldsymbol{\delta}_{\beta}\right)\right\rangle \right.;
    :label: beam_electron_wavefunction_at_z_eq_0

with :math:`\left|x,y\right\rangle` being the electron transverse position
eigenvector. As discussed in the documentation for the
:class:`prismatique.discretization.Params`, the total potential
:math:`V\left(\mathbf{r}\right)` is replaced by an effective potential
:math:`\tilde{V}\left(\mathbf{r}\right)`, which is obtained by coarse-graining
:math:`V\left(\mathbf{r}\right)` along the :math:`z`-axis:

.. math ::
    \tilde{V}\left(\mathbf{r}\right)=\sum_{n=0}^{N_{\text{slices}}-1}
    \delta\left(z-z_{n}\right)\tilde{V}_{n}\left(x,y\right),
    :label: coarse_graining_V

where :math:`N_{\text{slices}}` is a positive integer;

.. math ::
    z_{n}=n\delta z,
    :label: z_n_in_prism_params

with

.. math ::
    \delta z = \frac{\Delta Z}{N_{\text{slices}}},
    :label: slice_thickness_in_prism_params

:math:`\Delta Z` is the :math:`z`-dimension of the sample's supercell in units
of length. See the documentation for the class
:class:`prismatique.discretization.Params` for a discussion on sample supercells
and the calculation of :math:`\tilde{V}_{n}\left(x,y\right)`. Note that
:math:`z=z_{0}` and :math:`z=z_{N_{\text{slices}}}` are the
:math:`z`-coordinates of the entrance and exit surfaces of the sample's
supercell respectively. This coarse-graining approximation introduces a global
error that is linear in the slice thickness :math:`\delta z`. Upon replacing
:math:`V\left(\mathbf{r}\right)` with :math:`\tilde{V}\left(\mathbf{r}\right)`,
Eq. :eq:`fast_electron_equation` can be solved using the following iterative
scheme [Kirkland1]_:

.. math ::
    \zeta_{n+1}\left(x,y\right)=
    \left(\hat{K}_{\delta z}\circ\hat{C}_{n;\varphi=\sigma}\right)
    \left\{ \zeta_{n}\right\} \left(x,y\right),
    :label: zeta_nP1

where 

.. math ::
    \zeta_{n=0}\left(x,y\right)
    =\psi\left(x\hat{\mathbf{x}}+y\hat{\mathbf{y}}\right)
    \equiv\psi^{\left(\text{incident}\right)}\left(x,y\right),
    :label: zeta_0

with the entrance surface of the sample being located at :math:`z=0` and
:math:`\psi^{\left(\text{incident}\right)}\left(x,y\right)` being the incident
beam wavefunction;;

.. math ::
    \psi\left(x\hat{\mathbf{x}}+y\hat{\mathbf{y}}
    +\Delta Z\hat{\mathbf{z}}\right)&=\left\langle x,y\left|\delta_{f};
    \mathbf{u}_{1},\ldots,\mathbf{u}_{N};
    \boldsymbol{\delta}_{\beta}\right\rangle \right.\\
    &\equiv\psi^{\left(\text{exit}\right)}\left(x,y\right)\\
    &=\zeta_{n=N_{\text{slices}}}\left(x,y\right),
    :label: introduce_exit_wave

with the exit surface of the sample being located at
:math:`z=N_{\text{slices}}\delta z=\Delta Z` and
:math:`\psi^{\left(\text{exit}\right)}\left(x,y\right)` being the so-called exit
wavefunction; the notation :math:`\hat{O}\left\{ g\right\}
\left(q_{x},q_{y}\right)` denotes an operator :math:`\hat{O}` transforming a
function :math:`g` into a space with coordinates :math:`q_{x}` and :math:`q_{y}`
[e.g. :math:`\hat{O}` could be a two-dimensional Fourier transform]; the
notation :math:`\left(\hat{O}_{1}\circ\hat{O}_{2}\right)\left\{ g\right\}
\left(q_{x},q_{y}\right)` denotes the composition of two operators
:math:`\hat{O}_{1}` and :math:`\hat{O}_{2}` transforming a function :math:`g`
into a space with coordinates :math:`q_{1}` and :math:`q_{2}` with
:math:`\hat{O}_{2}` being applied first;

.. math ::
    \hat{K}_{\delta z^{\prime}}\left\{ g\right\} \left(x,y\right)=
    \left(\hat{A}\circ\hat{B}_{\delta z^{\prime}}\right)\left\{ g\right\} 
    \left(x,y\right),
    :label: K_operator

with

.. math ::
    \hat{A}\left\{ \tilde{g}\right\} \left(x,y\right)\equiv
    \hat{\mathcal{F}}_{\text{2D}}^{-1}\left\{ \tilde{g}\right\} 
    \left(x,y\right),
    :label: A_operator

.. math ::
    \hat{B}_{\delta z^{\prime}}\left\{ g\right\} \left(k_{x},k_{y}\right)
    \equiv\hat{\mathcal{F}}_{\text{2D}}\left\{ g\right\} 
    \left(k_{x},k_{y}\right)e^{-i\pi\lambda\left(k_{x}^{2}+k_{y}^{2}\right)
    \delta z^{\prime}},
    :label: B_operator

.. math ::
    \hat{C}_{n;\varphi}\left\{ g\right\} \left(x,y\right)\equiv 
    g\left(x,y\right)e^{i\varphi\tilde{V}_{n}\left(x,y\right)},
    :label: C_operator

.. math ::
    \hat{\mathcal{F}}_{\text{2D}}\left\{ g\right\} \left(k_{x},k_{y}\right)
    =\int_{-\infty}^{\infty}\int_{-\infty}^{\infty}dx\,dy\,g\left(x,y\right)
    e^{-2\pi i\left(k_{x}x+k_{y}y\right)},
    :label: fourier_transform_2D

and

.. math ::
    \hat{\mathcal{F}}_{\text{2D}}^{-1}\left\{ \tilde{g}\right\} 
    \left(x,y\right)=\int_{-\infty}^{\infty}\int_{-\infty}^{\infty}
    dk_{x}\,dk_{y}\,\tilde{g}\left(k_{x},k_{y}\right)
    e^{2\pi i\left(k_{x}x+k_{y}y\right)}.
    :label: inverse_fourier_transform_2D

The operator :math:`\hat{C}_{n;\varphi}` expresses the transformation that the
wavefunction :math:`\psi^{\left(\text{incident}\right)}\left(x,y\right)`
undergoes when a beam/probe electron passes through a 2D potential
:math:`\tilde{V}_{n}\left(x,y\right)` perdendicularly, i.e. the wavefunction
acquires a phase shift :math:`\varphi\tilde{V}_{n}\left(x,y\right)` while
keeping the amplitude the same. The operator :math:`\hat{K}_{\delta z^{\prime}}`
expresses the transformation that the wavefunction
:math:`\psi_{x_{p},y_{p}}\left(x,y,z\right)` undergoes when a beam/probe
electron propagates through free space in the :math:`z`-direction by a distance
:math:`\delta z^{\prime}`. In this discussion, we will refer to
:math:`\hat{K}_{\delta z^{\prime}}` as the free-space propagator. If
:math:`\psi^{\left(\text{incident}\right)}\left(x,y\right)` is the incident
wavefunction of a perfect coherent probe in a STEM experiment, then
Eqs. :eq:`zeta_nP1`-:eq:`inverse_fourier_transform_2D` essentially yield the
multislice algorithm for STEM simulations. 

The key observation involved in the PRISM algorithm is that the incident
wavefunction of a perfect coherent probe in a STEM experiment [Eq. :eq:`zeta_0`]
can be expanded in terms of plane waves as [Ophus1]_:

.. math ::
    \psi^{\left(\text{incident}\right)}\left(x,y\right)=
    \int_{-\infty}^{\infty}\int_{-\infty}^{\infty}dk_{x}\,dk_{y}\,
    \alpha_{k_{x},k_{y}}\left(x_{p},y_{p};\theta_{x},\theta_{y}\right)
    \phi_{k_{x},k_{y}}\left(x,y\right),
    :label: probe_expansion_1

where the :math:`\phi_{k_{x},k_{y}}\left(x^{\prime},y^{\prime}\right)` are the
plane wave basis functions:

.. math ::
    \phi_{k_{x},k_{y}}\left(x^{\prime},y^{\prime}\right)=
    e^{2\pi i\left(k_{x}x^{\prime}+k_{y}y^{\prime}\right)},
    :label: plane_wave

:math:`\alpha_{k_{x},k_{y}}\left(x_{p},y_{p};\theta_{x},\theta_{y}\right)` are
the expansion coefficients:

.. math ::
    \alpha_{k_{x},k_{y}}\left(x_{p},y_{p};\theta_{x},\theta_{y}\right)
    &=\phi_{k_{x},k_{y}}\left(-x_{p}+h\tan\left[\theta_{x}\right],
    -y_{p}+h\tan\left[\theta_{y}\right]\right)\\
    &\mathrel{\phantom{=}}\mathop{\times}
    e^{-i\chi\left(k_{x},k_{y}\right)}\zeta\left(k_{x},k_{y}\right),
    :label: alpha_expansion_coefficients

with :math:`\left(x_{p},y_{p}\right)` denoting the coordinates in the
:math:`xy`-plane of the probe position, :math:`\theta_{x}` and
:math:`\theta_{y}` being the probe tilt angles along the :math:`x`- and
:math:`y`-axes respectively, :math:`h` being the height [i.e.  the length along
the :math:`z`-axis] of the sample's supercell,
:math:`\chi\left(k_{x},k_{y}\right)` being the phase deviation due to coherent
lens aberrations, and :math:`\zeta\left(k_{x},k_{y}\right)` being an aperture
function. See the documentation for the class
:class:`embeam.coherent.Aberration` for a definition of
:math:`\chi\left(k_{x},k_{y}\right)`, and Ref. [DaCosta1]_ for the expression of
:math:`\zeta\left(k_{x},k_{y}\right)` implemented in ``prismatic``. One
important property of the aperture function, as defined in Ref. [DaCosta1]_, is
that for some :math:`k_{\text{max}}>0` we have

.. math ::
    \zeta\left(k_{x},k_{y}\right)=0,\quad\text{for }k_{x}^{2}+k_{y}^{2}
    \ge k_{\text{max}}^{2}.
    :label: aperture_function

At this stage, it is important to point out that both real-space and
momentum/Fourier-space need to be discretized in order to handle wavefunctions
numerically. The smallest possible Fourier-space pixel sizes in the
:math:`k_{x}`- and :math:`k_{y}`-directions, which we denote here as
:math:`\Delta k_{x}` and :math:`\Delta k_{y}` respectively, are determined by
the :math:`x`- and :math:`y`-dimensions of the sample's supercell:

.. math ::
    \Delta k_{x}=\frac{1}{\Delta X},\quad\Delta k_{y}=\frac{1}{\Delta Y},
    :label: Delta_ks_in_prism_params

where :math:`\Delta X` and :math:`\Delta Y` are the :math:`x`- and
:math:`y`-dimensions of the sample's supercell in units of length. It is worth
noting at this point that our system is subject to the following periodic
boundary conditions:

.. math ::
    &\Psi\left(x+l_{x}\Delta X\hat{\mathbf{x}}+y+l_{y}\Delta Y\hat{\mathbf{y}}
    +z\hat{\mathbf{z}}\right)\\&\quad=\Psi\left(x\hat{\mathbf{x}}
    +y\hat{\mathbf{y}}+z\hat{\mathbf{z}}\right),
    \quad\text{for }l_{x},l_{y}\in\mathbb{Z},
    :label: pbc_1

.. math ::
    &V\left(x+l_{x}\Delta X\hat{\mathbf{x}}+y+l_{y}\Delta Y\hat{\mathbf{y}}
    +z\hat{\mathbf{z}}\right)\\&\quad=V\left(x\hat{\mathbf{x}}
    +y\hat{\mathbf{y}}+z\hat{\mathbf{z}}\right),
    \quad\text{for }l_{x},l_{y}\in\mathbb{Z}.
    :label: pbc_2

In ``prismatique/prismatic``, for pixels of size :math:`\left(\Delta k_{x},
\Delta k_{y}\right)`, the :math:`k_{x}`- and :math:`k_{y}`-dimensions of
Fourier-space in units of said pixels is equivalent to half the :math:`x`- and
:math:`y`-dimensions of the sample's supercell in units of pixels. The reason
that the former has half the dimensions of the latter is because ``prismatic``
removes all momenta beyond the Nyquist limit, which is half of the original
Fourier-space pixel set, in order to prevent aliasing.

In discretized Fourier-space with pixels of size :math:`\left(\Delta k_{x},
\Delta k_{y}\right)`, we rewrite
:math:`\psi^{\left(\text{incident}\right)}\left(x,y\right)` as

.. math ::
    \psi^{\left(\text{incident}\right)}\left(x,y\right)\approx
    \sum_{m_{x}=-N_x/4}^{N_x/4-1}\sum_{m_{y}=-N_y/4}^{N_y/4-1}
    \alpha_{m_{x}\Delta k_{x},m_{y}\Delta k_{y}}
    \left(x_{p},y_{p};\theta_{x},\theta_{y}\right)\phi_{m_{x}
    \Delta k_{x},m_{y}\Delta k_{y}}\left(x,y\right),
    :label: probe_expansion_2

where :math:`N_x` and :math:`N_y` are the :math:`x`- and :math:`y`-dimensions of
the sample's supercell in units of pixels.

The PRISM algorithm involves expanding
:math:`\psi^{\left(\text{incident}\right)}\left(x,y\right)` in terms of plane
waves in a manner similar to that in Eq. :eq:`probe_expansion_2`:

.. math ::
    \psi^{\left(\text{incident}\right)}\left(x,y\right)\approx
    \sum_{m_{x}=-\tilde{N}_x}^{\tilde{N}_x-1}
    \sum_{m_{y}=-\tilde{N}_y}^{\tilde{N}_y-1}
    \alpha_{m_{x}f_{x}\Delta k_{x},m_{y}f_{x}
    \Delta k_{y}}\left(x_{p},y_{p};\theta_{x},\theta_{y}\right)
    \phi_{m_{x}f_{x}\Delta k_{x},m_{y}f_{y}\Delta k_{y}}\left(x,y\right),
    :label: probe_expansion_3

where :math:`f_{x}` and :math:`f_{y}` are positive integers called the :math:`x`
and :math:`y` interpolation factors introduced in the documentation for the
class :class:`prismatique.discretization.Params`, and

.. math ::
    \left(\tilde{N}_{x},\tilde{N}_{y}\right)=\left(\frac{N_{x}}{4f_{x}},
    \frac{N_{y}}{4f_{y}}\right),
    :label: reduced_supercell_xy_dims_in_pixels_in_prism_params

are what we refer to as the "sample supercell's reduced xy-dimensions in units
of pixels", also introduced in the documentation for the class
:class:`prismatique.discretization.Params`. Note that for :math:`f_{x}+f_{y}>2`,
Eq. :eq:`probe_expansion_3` approximates
:math:`\psi^{\left(\text{incident}\right)}\left(x,y\right)` further as compared
to Eq. :eq:`probe_expansion_2`. Hence, choosing interpolation factors that
satisfy :math:`f_{x}+f_{y}>2` will result in a loss of simulation accuracy since
we are discarding some plane wave basis functions, however the upside is an
improvement in simulation wall time. Further below we discuss how wall time
depends on :math:`f_{x}` and :math:`f_{y}`.

By discarding plane wave basis functions according to :math:`f_x` and
:math:`f_y`, and implementing the anti-alias procedure described above, we are
effectively working in a different discretized Fourier-space that is a subset of
the original discretized Fourier-space: the new Fourier-space pixel sizes in the
:math:`k_{x}`- and :math:`k_{y}`-directions are respectively

.. math ::
    \Delta \tilde{k}_{x}=f_x \Delta k_x,
    \quad\Delta \tilde{k}_{y}=f_y \Delta k_y,
    :label: Delta_tilde_ks_in_prism_params

and the :math:`k_{x}`- and :math:`k_{y}`-dimensions of this new discretized
Fourier-space in units of pixels are respectively

.. math ::
    n_x=2\tilde{N}_x,\quad n_y=2\tilde{N}_y.
    :label: stem_detector_dims_in_pixels

The pixelated STEM detector that is simulated in ``prismatic`` shares the same
pixel sizes and dimensions as this new discretized Fourier-space. Noting that a
beam electron's scattering angles with the :math:`x`- and :math:`y`-axes,
:math:`\theta_x` and :math:`\theta_y`, are related to the electron's momentum
by:
    
.. math ::
    \theta_x = \lambda k_x,\quad\theta_y = \lambda k_y,
    :label: k_to_theta_in_prism_params

it follows from Eqs. :eq:`Delta_ks_in_prism_params`, and
:eq:`reduced_supercell_xy_dims_in_pixels_in_prism_params`-:eq:`stem_detector_dims_in_pixels`
that the angular field of view of the aforementioned pixelated STEM detector in
the :math:`x`- and :math:`y`-directions are

.. math ::
    \Delta \boldsymbol{\Theta}_x 
    = \frac{2 \lambda f_x \tilde{N}_x}{\Delta X},\quad
    \Delta \boldsymbol{\Theta}_y
    = \frac{2 \lambda f_y \tilde{N}_y}{\Delta Y}.
    :label: angular_FOV_in_prism_params

Equations :eq:`Delta_tilde_ks_in_prism_params`-:eq:`angular_FOV_in_prism_params`
are also applicable to the case where the multislice algorithm is used, i.e. we
are still operating in this new discretized Fourier-space, however the
interpolation factors are trivially set to :math:`f_x=f_y=1`.

In the PRISM algorithm, we use Eq. :eq:`zeta_nP1` again but this time with

.. math ::
    \zeta_{n=0}\left(x,y\right)=\phi_{m_{x}f_{x}
    \Delta k_{x},m_{y}f_{y}\Delta k_{y}}\left(x,y\right),
    :label: prism_alg

and

.. math ::
    S_{m_{x},m_{y}}\left(x,y\right)=
    \zeta_{n=N_{\text{slices}}}\left(x,y\right),
    :label: S_matrix

for each :math:`\left(m_{x},m_{y}\right)` pair in the double sum of
Eq. :eq:`probe_expansion_3`. Once we have calculated all the necessary elements
of :math:`S_{m_{x},m_{y}}\left(x,y\right)`, we can then use it to calculate the
exit wavefunctions:

.. math ::
    \psi^{\left(\text{exit}\right)}\left(x,y\right)\approx
    \sum_{m_{x}=-\tilde{N}_x}^{\tilde{N}_x-1}
    \sum_{m_{y}=-\tilde{N}_y}^{\tilde{N}_y-1}
    \alpha_{m_{x}f_{x}\Delta k_{x},m_{y}f_{x}
    \Delta k_{y}}\left(x_{p},y_{p};\theta_{x},\theta_{y}\right)
    S_{m_{x},m_{y}}\left(x,y\right).
    :label: calc_exit_wave_from_S

with :math:`S_{m_{x},m_{y}}\left(x,y\right)` being the so-called
:math:`S`-matrix or scattering matrix. In other words, the plane wave basis
functions are transformed individually according to Eq. :eq:`zeta_nP1`, and then
added together using the same expansion coefficients
:math:`\alpha_{m_{x}f_{x}\Delta k_{x},m_{y}f_{x}\Delta
k_{y}}\left(x_{p},y_{p};\theta_{x},\theta_{y}\right)` in accordance with
Eq. :eq:`probe_expansion_3` to get the exit wavefunction
:math:`\psi^{\left(\text{exit}\right)}\left(x,y\right)`. In the PRISM algorithm,
a slight modification to Eq. :eq:`calc_exit_wave_from_S` where instead we
calculate

.. math ::
    \psi^{\left(\text{exit}\right)}\left(x,y\right)
    \approx
    \sum_{m_{x}=-\tilde{N}_x}^{\tilde{N}_x-1}
    \sum_{m_{y}=-\tilde{N}_y}^{\tilde{N}_y-1}\alpha_{m_{x}f_{x}
    \Delta k_{x},m_{y}f_{x}\Delta k_{y}}\left(x_{p},y_{p};
    \theta_{x},\theta_{y}\right)\tilde{S}_{m_{x},m_{y};x_{p},y_{p}}
    \left(x,y\right),
    :label: calc_exit_wave_from_S_tilde

with

.. math ::
    & \tilde{S}_{m_{x},m_{y};x_{p},y_{p}}\left(x,y\right)\\
    & \quad=\begin{cases}
    S_{m_{x},m_{y}}\left(x,y\right), & \text{if }\frac{\Delta X}{4f_{x}}
    \le x-x_{p}<\frac{\Delta X}{4f_{x}}\text{ and }\frac{\Delta Y}{4f_{y}}
    \le y-y_{p}<\frac{\Delta Y}{4f_{y}},\\
    0, & \text{otherwise},
    \end{cases}
    :label: S_tilde_matrix

i.e. :math:`\tilde{S}_{m_{x},m_{y};x_{p},y_{p}}\left(x,y\right)` is obtained by
cropping :math:`S_{m_{x},m_{y}}\left(x,y\right)` to a rectangular region that is
centered about :math:`\left(x_{p},y_{p}\right)`, and has :math:`x`- and
:math:`y`-dimensions in units of length smaller than those of the sample's
supercell by factors of :math:`4f_x` and :math:`4f_y` respectively.  Similarly
to discarding plane wave basis functions, the cropping of
:math:`S_{m_{x},m_{y}}\left(x,y\right)` will result in a loss of simulation
accuracy, especially for significant beam spread, however this loss comes with
faster wall times. To minimize the loss of simulation accuracy, ``prismatic``
offers the option to refocus the exit wavefunction upon calculating it as
described above.  To refocus the exit wavefunction, the free-space propagator
:math:`\hat{K}_{\delta z^{\prime}}` [Eq. :eq:`K_operator`] is applied to the
:math:`S`-matrix [DaCosta1]_:

.. math ::
    S_{m_{x},m_{y}}\left(x,y\right)\to\hat{K}_{\delta z^{\prime}}
    \left\{ S_{m_{x},m_{y}}\right\} \left(x,y\right),
    :label: refocusing_S

where the propagation displacement :math:`\delta z^{\prime}` is determined
automatically by ``prismatic`` such that the spread of the refocused exit
wavefunction has been reduced, thus reducing the loss of information due to
cropping. Since the free-space propagator imparts only a phase shift to the
Fourier components of the exit wavefunction, the corresponding diffraction
pattern — where only intensity measurements are recorded — is unaffected
[DaCosta1]_.

The main advantage of the PRISM algorithm is that
:math:`S_{m_{x},m_{y}}\left(x,y\right)` can be recycled to calculate the exit
wavefunction for each probe position :math:`\left(x_{p},y_{p}\right)`.  If we
assume that :math:`N_x=N_y\gg1`, then we have [Ophus1]_:

.. math ::
    \frac{t_{\text{Multi}}}{t_{\text{PRISM}}}\approx
    \frac{8\left(f_{x}f_{y}\right)^{2}N_{\text{slices}}N_{\text{probes}}}
    {N_{\text{PW}}
    \left(8f_{x}f_{y}N_{\text{slices}}+N_{\text{probes}}\right)},
    :label: comparing_algs_1

where :math:`t_{\text{Multi}}` and :math:`t_{\text{PRISM}}` are the wall times
for the multislice and PRISM algorithms respectively, :math:`N_{\text{probes}}`
is the number of probe positions considered, and :math:`N_{\text{PW}}` is the
number of plane waves considered in our expansion of
Eq. :eq:`probe_expansion_3`.  If
:math:`N_{\text{probes}}\gg8f_{x}f_{y}N_{\text{slices}}`, which is often the
case, then the performance of the PRISM algorithm will scale like:

.. math ::
    \frac{t_{\text{Multi}}}{t_{\text{PRISM}}}\approx
    \frac{8\left(f_{x}f_{y}\right)^{2}N_{\text{slices}}}{N_{\text{PW}}}.
    :label: comparing_algs_2

As mentioned above, the state operator :math:`\hat{\rho}_{t}` of a transmitted
beam electron is assumed to be a weighted sum of pure states
:math:`\hat{\rho}_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
\boldsymbol{\delta}_{\beta}\right)` over defocus, beam tilt, and frozen atomic
configurations. In ``prismatique``, the weighted sum is actually only over
defocus and frozen atomic configurations, not beam tilt, for STEM
simulations.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# Import child modules and packages of current package.
import prismatique.stem.system
import prismatique.stem.output
import prismatique.stem.sim



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = [""]



###########################
## Define error messages ##
###########################
