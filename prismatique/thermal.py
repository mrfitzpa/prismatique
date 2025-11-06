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
r"""For specifying simulation parameters related to the thermal properties of 
the sample and its environment.

Note that the documentation in this module draws from Ref. [Loane1]_, directly
copying certain passages when convenient.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For validating and converting objects.
import czekitout.check
import czekitout.convert

# For defining classes that support enforced validation, updatability,
# pre-serialization, and de-serialization.
import fancytypes



# For recycling helper functions and/or constants.
import prismatique.worker.cpu



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params"]



def _check_and_convert_enable_thermal_effects(params):
    obj_name = "enable_thermal_effects"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    enable_thermal_effects = czekitout.convert.to_bool(**kwargs)

    return enable_thermal_effects



def _pre_serialize_enable_thermal_effects(enable_thermal_effects):
    obj_to_pre_serialize = enable_thermal_effects
    serializable_rep = obj_to_pre_serialize
    
    return serializable_rep



def _de_pre_serialize_enable_thermal_effects(serializable_rep):
    enable_thermal_effects = serializable_rep

    return enable_thermal_effects



def _check_and_convert_num_frozen_phonon_configs_per_subset(params):
    obj_name = \
        "num_frozen_phonon_configs_per_subset"
    kwargs = \
        {"obj": params[obj_name], "obj_name": obj_name}
    num_frozen_phonon_configs_per_subset = \
        czekitout.convert.to_positive_int(**kwargs)

    return num_frozen_phonon_configs_per_subset



def _pre_serialize_num_frozen_phonon_configs_per_subset(
        num_frozen_phonon_configs_per_subset):
    obj_to_pre_serialize = num_frozen_phonon_configs_per_subset
    serializable_rep = obj_to_pre_serialize
    
    return serializable_rep



def _de_pre_serialize_num_frozen_phonon_configs_per_subset(serializable_rep):
    num_frozen_phonon_configs_per_subset = serializable_rep

    return num_frozen_phonon_configs_per_subset



def _check_and_convert_num_subsets(params):
    obj_name = "num_subsets"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    num_subsets = czekitout.convert.to_positive_int(**kwargs)

    return num_subsets



def _pre_serialize_num_subsets(num_subsets):
    obj_to_pre_serialize = num_subsets
    serializable_rep = obj_to_pre_serialize
    
    return serializable_rep



def _de_pre_serialize_num_subsets(serializable_rep):
    num_subsets = serializable_rep

    return num_subsets



def _check_and_convert_rng_seed(params):
    obj_name = "rng_seed"
    obj = params[obj_name]

    current_func_name = "_check_and_convert_rng_seed"
    
    if obj is not None:
        kwargs = {"obj": obj, "obj_name": obj_name}
        try:
            rng_seed = czekitout.convert.to_nonnegative_int(**kwargs)
        except:
            err_msg = globals()[current_func_name+"_err_msg_1"]
            raise TypeError(err_msg)
    else:
        rng_seed = obj

    return rng_seed



def _pre_serialize_rng_seed(rng_seed):
    obj_to_pre_serialize = rng_seed
    serializable_rep = obj_to_pre_serialize
    
    return serializable_rep



def _de_pre_serialize_rng_seed(serializable_rep):
    rng_seed = serializable_rep

    return rng_seed



_module_alias = \
    prismatique.worker.cpu
_default_enable_thermal_effects = \
    False
_default_num_frozen_phonon_configs_per_subset = \
    1
_default_num_subsets = \
    1
_default_rng_seed = \
    None
_default_skip_validation_and_conversion = \
    _module_alias._default_skip_validation_and_conversion



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related to the thermal properties of the 
    sample and its environment.

    In STEM experiments, the intensity pattern for a given probe collected by
    the STEM detector measures the diagonal elements of the state operator
    :math:`\hat{\rho}_{t}` of a transmitted beam electron, in the electron's
    transverse momentum basis:

    .. math ::
        I_{\text{STEM}}\left(k_{x},k_{y}\right)\Delta k_{x}\Delta k_{y}
        \approx N_{e}\left\langle k_{x},k_{y}\right|
        \hat{\rho}_{t}\left|k_{x},k_{y}\right\rangle \Delta k_{x}\Delta k_{y},
        :label: I_STEM

    where :math:`N_{e}` is the number of electrons collected by the detector for
    said intensity pattern; :math:`\left|k_{x},k_{y}\right\rangle` is the
    electron transverse momentum eigenvector; :math:`\Delta k_{x}` is the
    resolution of the STEM detector in the :math:`k_{x}`-direction;
    :math:`\Delta k_{y}` is the resolution of the STEM detector in the
    :math:`k_{y}`-direction;
    :math:`I_{\text{STEM}}\left(k_{x},k_{y}\right)\Delta k_{x}\Delta k_{y}` is
    the measured intensity, i.e. electron count, over the STEM detector pixel
    with the effective transverse momentum coordinates
    :math:`\left(k_{x},k_{y}\right)`. Similarly, in HRTEM experiments, the
    intensity image for a given beam collected by the TEM detector measures the
    diagonal elements of the state operator :math:`\hat{\rho}_{t}` of a
    transmitted beam electron, in the transverse (i.e. :math:`xy`) position
    basis:

    .. math ::
        I_{\text{HRTEM}}\left(x,y\right)\Delta x\Delta y\approx 
        N_{e}\left\langle x,y\right|\hat{\rho}_{t}\left|x,y\right\rangle 
        \Delta x\Delta y,
        :label: I_HRTEM

    where :math:`\left|x,y\right\rangle` is the electron transverse position
    eigenvector; :math:`\Delta x` is the resolution of the TEM detector in the
    :math:`x`-direction; :math:`\Delta y` is the resolution of the TEM detector
    in the :math:`y`-direction; :math:`I_{\text{HRTEM}}\left(x,y\right)\Delta
    x\Delta y` is the measured intensity, i.e. electron count, over the TEM
    detector pixel with the effective transverse position coordinates
    :math:`\left(x,y\right)`.

    Due to various sources of decoherence, a given transmitted beam electron
    will generally be in a mixed state. In ``prismatique``, we consider thermal
    fluctuations in the sample and chromatic aberrations as the dominant sources
    of decoherence for STEM and HRTEM simulations, where for the latter we also
    consider finite source size effects. At finite temperature, the atoms in the
    sample collectively vibrate to form phonons. Typically, the phonons have
    vibration periods on the order of :math:`10^{-13}\ \text{s}` [Loane1]_. A
    beam electron will typically pass through the sample in a time on the order
    of :math:`10^{-16}\ \text{s}` [Loane1]_, which is several orders of magnitude faster
    than the typical phonon vibration period. Hence, at a semiclassical level,
    any beam electron passing through the sample will interact with the
    vibrating atoms as if they were effectively stationary, i.e. the atoms are
    effectively in a frozen atomic or phonon configuration during the time
    interval it takes for a beam electron to pass through the sample. For
    typical beam currents, the average time between successive electrons passing
    through the specimen is about :math:`10^3` phonon vibration periods,
    therefore the effective frozen phonon configurations that successive beam
    electrons interact with are essentially uncorrelated. In scenarios where
    chromatic aberrations are present, and there are small fluctuations in the
    electron beam energy over time, there will be correspondingly small
    fluctuations in the defocus of the electron beam over time. In HRTEM
    experiments, the electron beam will generally illuminate the sample with a
    small distribution of angles, due to finite source size effects from the
    electron gun.

    To model the decoherence effects, we consider the following general state
    operator for a transmitted beam electron:

    .. math ::
        \hat{\rho}_{t}&=\int_{-\infty}^{\infty}d\delta_{f}\,
        \prod_{j=1}^{N}\left\{ \int d\mathbf{u}_{j}\right\} \,
        \int d\boldsymbol{\delta}_{\beta}\\
        &\phantom{=}\quad\mathop{\times}p_{\sigma_{f}}\left(\delta_{f}\right)
        p_{a}\left(\mathbf{u}_{1},\ldots,\mathbf{u}_{N}\right)
        p_{\sigma_{\beta}}\left(\boldsymbol{\delta}_{\beta}\right)\\
        &\phantom{=}\quad\mathop{\times}
        \hat{\rho}_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
        \boldsymbol{\delta}_{\beta}\right),
        :label: mixed_state_operator_for_transmitted_electron

    where :math:`\delta_{f}` is the defocal offset, i.e. a deviation from the
    target operating defocus :math:`\Delta f`;
    :math:`p_{\sigma_{f}}\left(\delta_{f}\right)` is the distribution of
    :math:`\delta_{f}`:

    .. math ::
        p_{\sigma_{f}}\left(\delta_{f}\right)=\frac{1}{\sigma_{f}\sqrt{2\pi}}
        \exp\left(-\frac{1}{2}\frac{\delta_{f}^{2}}{\sigma_{f}^{2}}\right),
        :label: p_sigma_f

    with :math:`\sigma_{f}` being the defocal spread;
    :math:`\boldsymbol{\delta}_{\beta}` is the deviation from the target
    operating beam tilt :math:`\mathbf{k}_{xy,\beta}`;
    :math:`p_{\sigma_{\beta}}\left(\boldsymbol{\delta}_{\beta}\right)` is the
    distribution of :math:`\boldsymbol{\delta}_{\beta}`:

    .. math ::
        p_{\sigma_{\beta}}\left(\boldsymbol{\delta}_{\beta}\right)=
        \frac{1}{\sigma_{\beta}\sqrt{2\pi}}
        \exp\left(-\frac{1}{2}
        \frac{\left|\boldsymbol{\delta}_{\beta}\right|^{2}}{\sigma_{\beta}^{2}}
        \right),
        :label: p_sigma_beta

    with :math:`\sigma_{\beta}` being the beam tilt spread; :math:`N` is the
    number of atoms in the sample's supercell [see the documentation for the
    class :class:`prismatique.discretization.Params` for a discussion on sample
    supercells]; :math:`\mathbf{u}_{j}` is the displacement of the
    :math:`j^{\text{th}}` atom of the sample's supercell from the expectation
    value of its position at zero-temperature;
    :math:`p_{a}\left(\mathbf{u}_{1},\ldots,\mathbf{u}_{N}\right)` is the
    distribution of the sample's supercell atomic configuration :math:`\left\{
    \mathbf{u}_{j}\right\} _{j=1}^{N}`; and

    .. math ::
        &\hat{\rho}_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
        \boldsymbol{\delta}_{\beta}\right)\\
        &\quad=\left|\psi_{t}\left(\delta_{f};
        \mathbf{u}_{1},\ldots,\mathbf{u}_{N};
        \boldsymbol{\delta}_{\beta}\right)\right\rangle \left\langle 
        \psi_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
        \boldsymbol{\delta}_{\beta}\right)\right|,
        :label: pure_state_operator_for_transmitted_electron_1

    with
    :math:`\left|\psi_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
    \boldsymbol{\delta}_{\beta}\right)\right\rangle` being the state vector of a
    transmitted beam electron for a perfectly coherent beam operating at a
    defocus of :math:`\Delta f+\delta_{f}` and a beam tilt of
    :math:`\mathbf{k}_{xy,\beta}+\boldsymbol{\delta}_{\beta}`, and a sample with
    a frozen supercell atomic configuration :math:`\left\{
    \mathbf{u}_{j}\right\} _{j=1}^{N}`. The defocal spread :math:`\sigma_{f}` is
    calculated by:

    .. math ::
        \sigma_{f}=C_{c}\sqrt{\left(\frac{\sigma_{E}/e}{V}\right)^{2}
        +\left(2\frac{\sigma_{I}}{I}\right)^{2}
        +\left(\frac{\sigma_{V}}{V}\right)^{2}},
        :label: sigma_f_in_stem_probe_model_params__1

    where :math:`C_c` is the chromatic aberration coefficient; :math:`I` is the
    mean current of the lens [the probe forming lens in STEM simulations and the
    objective lens in HRTEM simulations]; :math:`\sigma_I` is the standard
    deviation of the current of the lens; :math:`V` is the mean accelerating
    voltage; :math:`e` is the elementary charge; :math:`\sigma_V` is the
    standard deviation of the accelerating voltage; and :math:`\sigma_E` is the
    standard deviation of the electrons in the gun when operating a voltage
    supply that does not fluctuate [i.e. :math:`\sigma_E` is the intrinsic
    energy spread of the gun]. As mentioned above, we do not explicitly consider
    finite source size effects for STEM simulations in ``prismatique``, hence
    for such scenarios we can set

    .. math ::
        \sigma_{\beta}\to0,\quad\left(\text{for STEM}\right).
        :label: sigma_beta_to_zero_for_STEM

    In ``prismatic``, :math:`p_{a}\left(
    \mathbf{u}_{1},\ldots,\mathbf{u}_{N}\right)` is approximated using a simple
    Einstein model [Loane1]_:

    .. math ::
        p_{a}\left(\mathbf{u}_{1},\ldots,\mathbf{u}_{N}\right)=\prod_{i=1}^{N}
        \left\{ p_{a}\left(\left|\mathbf{u}_{i}\right|\right)\right\},
        :label: p_a_for_atomic_config

    .. math ::
        p_{a}\left(\mathbf{u}_{i}\right)=
        \left(\frac{3}{2\pi u_{i,\text{rms}}^{2}}\right)^{3/2}
        e^{-\frac{3}{2}\left(
        \frac{\mathbf{u}_{i}}{u_{i,\text{rms}}}\right)^{2}},
        :label: p_a_for_single_atom

    .. math ::
        u_{i,\text{rms}}=\sqrt{\left\langle \mathbf{u}_{i}^{2}\right\rangle }
        = \frac{A\left(T\right)}{m_i^{\frac{1}{2}}},
        :label: u_i_rms

    where :math:`A\left(T\right)` is some function of the temperature :math:`T`
    which would need to be determined experimentally, and :math:`m_i` is the
    mass of the :math:`i^{\text{th}}` atom. The incoherent average over the
    phonon configurations in
    Eq. :eq:`mixed_state_operator_for_transmitted_electron` is approximated by
    sampling a random set of frozen phonon configurations from the distribution
    :math:`p_{a}\left(\mathbf{u}_{1},\ldots,\mathbf{u}_{N}\right)`. The larger
    the sampling size of frozen phonon configurations, the more accurate the
    estimate. For a discussion on the convergence of this sampling approach, see
    Ref. [DaCosta1]_, in particular Sec. 4.4. Note that throughout the
    documentation of the ``prismatique`` library, we refer to
    :math:`u_{i,\text{rms}}` as the effective root-mean-squared displacement of
    the :math:`i^{\text{th}}` atom.

    Parameters
    ----------
    enable_thermal_effects : `bool`, optional
        If ``enable_thermal_effects`` is set to ``True``, then for each atom in
        the sample supercell, :math:`u_{i,\text{rms}}` is set according to the
        specifications given in the "atomic coordinates file", which among other
        things specifies the zero-temperature expectation value of the atomic
        coordinates for each atom in the unit cell of the sample. If
        ``enable_thermal_effects`` is set to ``False``, then for each atom,
        :math:`u_{i,\text{rms}}` is set to zero. See the documentation for the
        class :class:`prismatique.sample.ModelParams`, in particular the
        documentation for the construction parameter ``atomic_coords_filename``
        for a description of the file format of the atomic coordinates file.
    num_frozen_phonon_configs_per_subset : `int`, optional
        Sometimes one might be implementing a simulation that requires
        significant computational resources, namely RAM. In such cases, rather
        than call a single session of ``prismatic`` to calculate
        :math:`p_{\text{STEM}}\left(k_{x},k_{y}\left|\mathbf{u}_{1},\ldots,
        \mathbf{u}_{N}\right.\right)` for all sampled frozen phonon
        configurations, it is advantageous to partition the set of sampled
        frozen phonon configurations into disjoint subsets and subsequently call
        multiple sessions of ``prismatic`` sequentially, wherein for each
        session
        :math:`p_{\text{STEM}}\left(k_{x},k_{y}\left|\mathbf{u}_{1},\ldots,
        \mathbf{u}_{N}\right.\right)` is calculated for a different subset of
        the sampled frozen phonon configurations. This way, one need only store
        the potential slice data or the S-matrix data [depending on the
        algorithm used] for a single subset of frozen phonon configurations at a
        time, thus reducing memory requirements. See the documentation for the
        class :class:`prismatique.sample.ModelParams` for a discussion on
        potential slices, and the subpackage :mod:`prismatique.stem` for a
        discussion on S-matrices.

        ``num_frozen_phonon_configs_per_subset`` is the number of frozen phonon
        configurations to sample per subset.
    num_subsets : `int`, optional
        Continuing from above, ``num_subsets`` is the number of subsets of 
        frozen phonon configurations. The total number of frozen phonon 
        configurations to sample is 
        ``num_frozen_phonon_configs_per_subset*num_subsets``.
    rng_seed : `int` | `None`, optional
        A seed to pass to the random number generator used to sample the frozen
        phonon configurations. If set, ``rng_seed`` must be a non-negative 
        integer or `None`.
    skip_validation_and_conversion : `bool`, optional
        Let ``validation_and_conversion_funcs`` and ``core_attrs`` denote the
        attributes :attr:`~fancytypes.Checkable.validation_and_conversion_funcs`
        and :attr:`~fancytypes.Checkable.core_attrs` respectively, both of which
        being `dict` objects.

        Let ``params_to_be_mapped_to_core_attrs`` denote the `dict`
        representation of the constructor parameters excluding the parameter
        ``skip_validation_and_conversion``, where each `dict` key ``key`` is a
        different constructor parameter name, excluding the name
        ``"skip_validation_and_conversion"``, and
        ``params_to_be_mapped_to_core_attrs[key]`` would yield the value of the
        constructor parameter with the name given by ``key``.

        If ``skip_validation_and_conversion`` is set to ``False``, then for each
        key ``key`` in ``params_to_be_mapped_to_core_attrs``,
        ``core_attrs[key]`` is set to ``validation_and_conversion_funcs[key]
        (params_to_be_mapped_to_core_attrs)``.

        Otherwise, if ``skip_validation_and_conversion`` is set to ``True``,
        then ``core_attrs`` is set to
        ``params_to_be_mapped_to_core_attrs.copy()``. This option is desired
        primarily when the user wants to avoid potentially expensive deep copies
        and/or conversions of the `dict` values of
        ``params_to_be_mapped_to_core_attrs``, as it is guaranteed that no
        copies or conversions are made in this case.

    """
    ctor_param_names = ("enable_thermal_effects",
                        "num_frozen_phonon_configs_per_subset",
                        "num_subsets",
                        "rng_seed")
    kwargs = {"namespace_as_dict": globals(),
              "ctor_param_names": ctor_param_names}
    
    _validation_and_conversion_funcs_ = \
        fancytypes.return_validation_and_conversion_funcs(**kwargs)
    _pre_serialization_funcs_ = \
        fancytypes.return_pre_serialization_funcs(**kwargs)
    _de_pre_serialization_funcs_ = \
        fancytypes.return_de_pre_serialization_funcs(**kwargs)

    del ctor_param_names, kwargs

    
    
    def __init__(self,
                 enable_thermal_effects=\
                 _default_enable_thermal_effects,
                 num_frozen_phonon_configs_per_subset=\
                 _default_num_frozen_phonon_configs_per_subset,
                 num_subsets=\
                 _default_num_subsets,
                 rng_seed=\
                 _default_rng_seed,
                 skip_validation_and_conversion=\
                 _default_skip_validation_and_conversion):
        ctor_params = {key: val
                       for key, val in locals().items()
                       if (key not in ("self", "__class__"))}
        kwargs = ctor_params
        kwargs["skip_cls_tests"] = True
        fancytypes.PreSerializableAndUpdatable.__init__(self, **kwargs)

        return None



    @classmethod
    def get_validation_and_conversion_funcs(cls):
        validation_and_conversion_funcs = \
            cls._validation_and_conversion_funcs_.copy()

        return validation_and_conversion_funcs


    
    @classmethod
    def get_pre_serialization_funcs(cls):
        pre_serialization_funcs = \
            cls._pre_serialization_funcs_.copy()

        return pre_serialization_funcs


    
    @classmethod
    def get_de_pre_serialization_funcs(cls):
        de_pre_serialization_funcs = \
            cls._de_pre_serialization_funcs_.copy()

        return de_pre_serialization_funcs



def _check_and_convert_thermal_params(params):
    obj_name = "thermal_params"
    obj = params[obj_name]

    accepted_types = (Params, type(None))
    
    if isinstance(obj, accepted_types[-1]):
        thermal_params = accepted_types[0]()
    else:
        kwargs = {"obj": obj,
                  "obj_name": obj_name,
                  "accepted_types": accepted_types}
        czekitout.check.if_instance_of_any_accepted_types(**kwargs)

        kwargs = obj.get_core_attrs(deep_copy=False)
        thermal_params = accepted_types[0](**kwargs)

    return thermal_params



def _pre_serialize_thermal_params(thermal_params):
    obj_to_pre_serialize = thermal_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_thermal_params(serializable_rep):
    thermal_params = Params.de_pre_serialize(serializable_rep)
    
    return thermal_params



_default_thermal_params = None
    


###########################
## Define error messages ##
###########################

_check_and_convert_rng_seed_err_msg_1 = \
    ("The object ``rng_seed`` must be either a nonnegative integer or of the "
     "type `NoneType`.")
