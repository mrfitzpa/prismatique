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
r"""For specifying simulation parameters and calculating quantities related to
beam tilt series in HRTEM.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For general array handling.
import numpy as np

# For validating and converting objects.
import czekitout.check
import czekitout.convert

# For defining classes that support enforced validation, updatability,
# pre-serialization, and de-serialization.
import fancytypes

# For calculating the electron beam wavelength and validating attributes of the
# class :class:`embeam.gun.ModelParams`.
import embeam



# For validating instances of the classes
# :class:`prismatique.sample.ModelParams`, and
# :class:`prismatique.sample.PotentialSliceSubsetIDs``; and for calculating
# quantities related to the modelling of the sample.
import prismatique.sample



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params",
           "step_size",
           "series"]



def _check_and_convert_offset(params):
    obj_name = "offset"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    offset = czekitout.convert.to_pair_of_floats(**kwargs)

    return offset



def _pre_serialize_offset(offset):
    obj_to_pre_serialize = offset
    serializable_rep = obj_to_pre_serialize
    
    return serializable_rep



def _de_pre_serialize_offset(serializable_rep):
    offset = serializable_rep

    return offset



def _check_and_convert_window(params):
    obj_name = "window"

    current_func_name = "_check_and_convert_window"

    try:
        kwargs = {"obj": params[obj_name], "obj_name": obj_name}
        window = czekitout.convert.to_tuple_of_nonnegative_floats(**kwargs)

        if (len(window) != 2) and (len(window) != 4):
            raise
        
        if ((not (0 <= window[0] <= window[1]))
            or (not (0 <= window[-2] <= window[-1]))):
            raise
    except:
        err_msg = globals()[current_func_name+"_err_msg_1"]
        raise TypeError(err_msg)

    return window



def _pre_serialize_window(window):
    obj_to_pre_serialize = window
    serializable_rep = obj_to_pre_serialize
    
    return serializable_rep



def _de_pre_serialize_window(serializable_rep):
    window = serializable_rep

    return window



def _check_and_convert_spread(params):
    obj_name = "spread"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    spread = czekitout.convert.to_nonnegative_float(**kwargs)

    return spread



def _pre_serialize_spread(spread):
    obj_to_pre_serialize = spread
    serializable_rep = obj_to_pre_serialize
    
    return serializable_rep



def _de_pre_serialize_spread(serializable_rep):
    spread = serializable_rep

    return spread



_module_alias = \
    prismatique.sample
_default_offset = \
    (0, 0)
_default_window = \
    (0, 0)
_default_spread = \
    0
_default_skip_validation_and_conversion = \
    _module_alias._default_skip_validation_and_conversion



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related to the beam tilt series in the HRTEM
    simulation.

    As discussed in the documentation for the class
    :class:`prismatique.thermal.Params`, due to finite source effects from the
    electron gun, the electron beam in HRTEM experiments will generally
    illuminate the sample with a small distribution of angles. The resulting
    decoherence effects are encoded in the expression we use for the state
    operator for a transmitted beam electron,
    Eq. :eq:`mixed_state_operator_for_transmitted_electron`, wherein an
    incoherent averaging over a set of beam tilt angles is performed.

    The incident beam's tilt angles with the :math:`x`- and :math:`y`-axes,
    :math:`\theta_x` and :math:`\theta_y`, are related to the incident beam
    electron momentum by:

    .. math ::
        \theta_x = \lambda k_x,\quad\theta_y = \lambda k_y,
        :label: k_to_theta_in_tilt_params
    
    where :math:`\lambda` is the incident beam electron wavelength, and
    :math:`k_x` and :math:`k_y` are the :math:`x`- and :math:`y`-components of
    the incident beam electron momentum.

    As discussed in the documentation for the class
    :class:`prismatique.discretization.Params`, real-space [and thus
    momentum/Fourier/angular-space] need to be discretized in order to handle
    wavefunctions and probabilities numerically. The smallest possible
    Fourier-space pixel sizes in the :math:`k_{x}`- and
    :math:`k_{y}`-directions, which we denote here as :math:`\Delta k_{x}` and
    :math:`\Delta k_{y}` respectively, are determined by the :math:`x`- and
    :math:`y`-dimensions of the sample's supercell:

    .. math ::
        \Delta k_{x}=\frac{1}{\Delta X},\quad\Delta k_{y}=\frac{1}{\Delta Y},
        :label: Delta_ks_in_tilt_params

    where :math:`\Delta X` and :math:`\Delta Y` are the :math:`x`- and
    :math:`y`-dimensions of the sample's supercell in units of length [see the
    documentation for the class :class:`prismatique.discretization.Params` for a
    discussion on supercells]. An angular space
    :math:`\boldsymbol{\Theta}_{f_{x},f_{y}}` that is useful to our current
    discussion is the set of angles:

    .. math ::
        \boldsymbol{\Theta}_{f_{x},f_{y}}=\left\{ \left.\left(l_{x}f_{x}\lambda
        \Delta k_{x},l_{y}f_{y}\lambda
        \Delta k_{y}\right)\right|l_{x},l_{y}\in\mathbb{Z}\right\},
        :label: discretized_angular_space

    where :math:`f_{x}` and :math:`f_{y}` are positive integers called the
    :math:`x` and :math:`y` interpolation factors, which are introduced in the
    documentation for the class :class:`prismatique.discretization.Params`.

    In ``prismatique``, users can select all the angles from either a
    rectangular or radial window/region of the discretized angular-space
    :math:`\boldsymbol{\Theta}_{f_{x},f_{y}}` as a set of beam tilt angles for
    which to perform HRTEM simulations. Users can generate beam tilt series to
    conveniently model spatially coherent HRTEM experiments at different beam
    tilts, or to model a single spatially incoherent HRTEM beam, for which to
    use in simulations.

    Parameters
    ----------
    offset : array_like` (`float`, shape=(``2``,)), optional
        ``offset`` specifies the offset of the window for the beam tilt series
        discussed above: ``offset[0]`` specifies the :math:`x`-coordinate of the
        offset in mrads; ``offset[1]`` specifies the :math:`y`-coordinate of the
        offset in mrads.
    window : array_like` (`float`, shape=(``2``,)) | array_like` (`float`, shape=(``4``,)), optional
        If ``window`` is an array of length 2, then ``window`` specifies a
        radial window for the beam tilt series: ``window[0]`` and ``window[1]``
        specify the minimum and maximum possible radial tilts with respect to 
        the tilt offset in mrads.

        If ``window`` is an array of length 4, then ``window`` specifies a
        rectangular window for the beam tilt series: ``window[0]`` and
        ``window[1]`` specify the minimum and maximum possible absolute values
        of the :math:`x`-tilt with respect to the tilt offset in mrads;
        ``window[2]`` and ``window[3]`` specify the minimum and maximum possible
        absolute values of the :math:`y`-tilt with respect to the tilt offset in
        mrads.
    spread : `float`, optional
        The beam tilt spread :math:`\sigma_{\beta}` in units of mrads,
        introduced in Eq. :eq:`p_sigma_beta`. Must be nonnegative.
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
    ctor_param_names = ("offset",
                        "window",
                        "spread")
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
                 offset=\
                 _default_offset,
                 window=\
                 _default_window,
                 spread=\
                 _default_spread,
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



def _check_and_convert_tilt_params(params):
    obj_name = "tilt_params"
    obj = params[obj_name]

    accepted_types = (Params, type(None))
    
    if isinstance(obj, accepted_types[-1]):
        tilt_params = accepted_types[0]()
    else:
        kwargs = {"obj": obj,
                  "obj_name": obj_name,
                  "accepted_types": accepted_types}
        czekitout.check.if_instance_of_any_accepted_types(**kwargs)

        kwargs = obj.get_core_attrs(deep_copy=False)
        tilt_params = accepted_types[0](**kwargs)

    return tilt_params



def _pre_serialize_tilt_params(tilt_params):
    obj_to_pre_serialize = tilt_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_tilt_params(serializable_rep):
    tilt_params = Params.de_pre_serialize(serializable_rep)
    
    return tilt_params



_default_tilt_params = None



def _check_and_convert_sample_specification(params):
    params["accepted_types"] = (prismatique.sample.ModelParams,
                                prismatique.sample.PotentialSliceSubsetIDs)

    module_alias = prismatique.sample
    func_alias = module_alias._check_and_convert_sample_specification
    sample_specification = func_alias(params)

    del params["accepted_types"]

    return sample_specification



def _check_and_convert_mean_beam_energy(params):
    obj_name = "mean_beam_energy"

    module_alias = embeam.gun
    cls_alias = module_alias.ModelParams
    func_alias = cls_alias.get_validation_and_conversion_funcs()[obj_name]
    mean_beam_energy = func_alias(params)

    return mean_beam_energy



def _check_and_convert_skip_validation_and_conversion(params):
    module_alias = prismatique.sample
    func_alias = module_alias._check_and_convert_skip_validation_and_conversion
    skip_validation_and_conversion = func_alias(params)

    return skip_validation_and_conversion



_module_alias = \
    prismatique.sample
_default_mean_beam_energy = \
    80
_default_skip_validation_and_conversion = \
    _module_alias._default_skip_validation_and_conversion



def step_size(sample_specification,
              mean_beam_energy=\
              _default_mean_beam_energy,
              skip_validation_and_conversion=\
               _default_skip_validation_and_conversion):
    r"""Determine beam tilt series step size in the HRTEM simulation from a 
    subset of the simulation parameters.

    For additional context on the tilt step size, see the documentation for the
    class :class:`prismatique.tilt.Params`.

    Parameters
    ----------
    sample_specification : :class:`prismatique.sample.ModelParams` | :class:`prismatique.sample.PotentialSliceSubsetIDs`
        The simulation parameters specifying the sample model. 

        If ``sample_specification`` is of the type
        :class:`prismatique.sample.ModelParams`, then ``sample_specifications``
        specifies sample model parameters that are used to construct the model
        from scratch, i.e. the potential slices for each frozen phonon
        configuration subset are calculated from said model parameters. See the
        documentation for the classes :class:`prismatique.discretization.Params`
        and :class:`prismatique.thermal.Params` for discussions on potential
        slices and frozen phonon configuration subsets respectively. Note that
        of parameters stored in ``sample_specification``, only the following are
        used:

        - sample_specification

          * atomic_coords_filename

          * unit_cell_tiling

          * discretization_params

            + interpolation_factors

        Otherwise, if ``sample_specification`` is an instance of the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs` then
        ``sample_specification`` specifies a set of files, where each file
        stores the pre-calculated potential slices for a frozen phonon
        configuration subset. See the documentation for the aforementioned
        class for a further discussion on specifying pre-calculated
        potential slices. 
    mean_beam_energy : `float`, optional
        The mean electron beam energy :math:`E` in units of keV. Must be
        positive.
    skip_validation_and_conversion : `bool`, optional
        If ``skip_validation_and_conversion`` is set to ``False``, then
        validations and conversions are performed on the above
        parameters. 

        Otherwise, if ``skip_validation_and_conversion`` is set to ``True``, no
        validations and conversions are performed on the above parameters. This
        option is desired primarily when the user wants to avoid potentially
        expensive validation and/or conversion operations.

    Returns
    -------
    tilt_step_size : `array_like` (`float`, shape=(``2``,))
        ``tilt_step_size[0]`` and ``tilt_step_size[1]`` are the tilt step sizes 
        in the :math:`x`- and :math:`y`-directions in units of mrads 
        respectively.

    """
    params = locals()

    func_alias = _check_and_convert_skip_validation_and_conversion
    skip_validation_and_conversion = func_alias(params)

    if (skip_validation_and_conversion == False):
        global_symbol_table = globals()
        for param_name in params:
            if param_name == "skip_validation_and_conversion":
                continue
            func_name = "_check_and_convert_" + param_name
            func_alias = global_symbol_table[func_name]
            params[param_name] = func_alias(params)

    kwargs = params
    del kwargs["skip_validation_and_conversion"]
    tilt_step_size = _step_size(**kwargs)

    return tilt_step_size



def _step_size(sample_specification, mean_beam_energy):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    if "thermal_params" in sample_specification_core_attrs:
        discretization_params = \
            sample_specification_core_attrs["discretization_params"]
        discretization_params_core_attrs = \
            discretization_params.get_core_attrs(deep_copy=False)
        interpolation_factors = \
            np.array(discretization_params_core_attrs["interpolation_factors"])
    else:
        interpolation_factors = \
            np.array(sample_specification_core_attrs["interpolation_factors"])

    kwargs = \
        {"sample_specification": sample_specification,
         "mean_beam_energy": mean_beam_energy}
    smallest_possible_angular_space_pixel_size = \
        _calc_smallest_possible_angular_space_pixel_size(**kwargs)
    
    tilt_step_size = \
        smallest_possible_angular_space_pixel_size * interpolation_factors
    tilt_step_size = \
        tuple(float(elem) for elem in tilt_step_size)

    return tilt_step_size



def _calc_smallest_possible_angular_space_pixel_size(sample_specification,
                                                     mean_beam_energy):
    kwargs = {"sample_specification": sample_specification}
    sample_supercell_dims = prismatique.sample._supercell_dims(**kwargs)

    kwargs = {"beam_energy": mean_beam_energy,
              "skip_validation_and_conversion": True}
    electron_beam_wavelength = embeam.wavelength(**kwargs)
    
    result = electron_beam_wavelength / sample_supercell_dims[:2] * 1000
    result = tuple(float(elem) for elem in result)

    return result



def series(sample_specification,
           mean_beam_energy=\
           _default_mean_beam_energy,
           tilt_params=\
           _default_tilt_params,
           skip_validation_and_conversion=\
           _default_skip_validation_and_conversion):
    r"""Determine the beam tilt series in the HRTEM simulation from a subset of 
    the simulation parameters.

    For additional context on beam tilts, see the documentation for the class
    :class:`prismatique.tilt.Params`.

    Parameters
    ----------
    sample_specification : :class:`prismatique.sample.ModelParams` | :class:`prismatique.sample.PotentialSliceSubsetIDs`
        The simulation parameters specifying the sample model. 

        If ``sample_specification`` is of the type
        :class:`prismatique.sample.ModelParams`, then ``sample_specifications``
        specifies sample model parameters that are used to construct the model
        from scratch, i.e. the potential slices for each frozen phonon
        configuration subset are calculated from said model parameters. See the
        documentation for the classes :class:`prismatique.discretization.Params`
        and :class:`prismatique.thermal.Params` for discussions on potential
        slices and frozen phonon configuration subsets respectively. Note that
        of parameters stored in ``sample_specification``, only the following are
        used:

        - sample_specification

          * atomic_coords_filename

          * unit_cell_tiling

          * discretization_params

            + sample_supercell_reduced_xy_dims_in_pixels

            + interpolation_factors

        Otherwise, if ``sample_specification`` is an instance of the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs` then
        ``sample_specification`` specifies a set of files, where each file
        stores the pre-calculated potential slices for a frozen phonon
        configuration subset. See the documentation for the aforementioned
        class for a further discussion on specifying pre-calculated
        potential slices. 
    mean_beam_energy : `float`, optional
        The mean electron beam energy :math:`E` in units of keV. Must be
        positive.
    tilt_params : :class:`prismatique.tilt.Params` | `None`, optional
        The simulation parameters related to the beam tilt series in the HRTEM
        simulation to model a set of spatially coherent HRTEM experiments at
        different beam tilts, or to model a single spatially incoherent HRTEM
        beam. See the documentation for the class
        :class:`prismatique.tilt.Params` for a discussion on said parameters.
        If ``tilt_params`` is set to `None` [i.e. the default value], then the
        aforementioned simulation parameters are set to default values.
    skip_validation_and_conversion : `bool`, optional
        If ``skip_validation_and_conversion`` is set to ``False``, then
        validations and conversions are performed on the above
        parameters. 

        Otherwise, if ``skip_validation_and_conversion`` is set to ``True``, no
        validations and conversions are performed on the above parameters. This
        option is desired primarily when the user wants to avoid potentially
        expensive validation and/or conversion operations.

    Returns
    -------
    tilt_series : `array_like` (`float`, shape=(``num_tilts``, ``2``))
        If we let ``num_tilts`` be the number of beam tilts, then
        ``tilt_series[i][0]`` and ``tilt_series[i][1]`` are the beam tilt angles
        along the :math:`x`- and :math:`y`-axes respectively for the ``i`` th
        beam tilt in units of mrad, where ``0<=i<num_tilts``.

    """
    func_alias = _check_and_convert_skip_validation_and_conversion
    skip_validation_and_conversion = func_alias(params)

    if (skip_validation_and_conversion == False):
        global_symbol_table = globals()
        for param_name in params:
            if param_name == "skip_validation_and_conversion":
                continue
            func_name = "_check_and_convert_" + param_name
            func_alias = global_symbol_table[func_name]
            params[param_name] = func_alias(params)

    kwargs = params
    del kwargs["skip_validation_and_conversion"]
    tilt_series = _series(**kwargs)

    return tilt_series



def _series(sample_specification, mean_beam_energy, tilt_params):
    tilt_params_core_attrs = tilt_params.get_core_attrs(deep_copy=False)

    kwargs = \
        {"sample_specification": sample_specification,
         "mean_beam_energy": mean_beam_energy,
         "tilt_offset": tilt_params_core_attrs["offset"],
         "tilt_window": tilt_params_core_attrs["window"],
         "for_a_hrtem_calculation": True}
    angular_mesh, beam_mask = \
        prismatique.sample._calc_angular_mesh_and_beam_mask(**kwargs)

    tilt_series = []
    for k_x_idx in range(beam_mask.shape[0]):
        for k_y_idx in range(beam_mask.shape[1]):
            if beam_mask[k_x_idx][k_y_idx]:
                x_tilt = angular_mesh[0][k_x_idx][k_y_idx]  # In rads.
                y_tilt = angular_mesh[1][k_x_idx][k_y_idx]  # In rads.
                tilt_series.append((x_tilt*1000, y_tilt*1000))  # In mrads.

    tilt_series = np.array(tilt_series)
    reordered_indices = np.lexsort((tilt_series[:, 1], tilt_series[:, 0]))
    tilt_series = tilt_series[reordered_indices]

    return tilt_series



###########################
## Define error messages ##
###########################

_check_and_convert_window_err_msg_1 = \
    ("The object ``window`` must be an array-like object of length 2 or 4, "
     "with each element being of type `float`, and the array as a whole "
     "satisfying: ``0<=window[0]<=window[1]`` and "
     "``0<=window[-2]<=window[-1]``.")
