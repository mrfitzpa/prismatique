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
r"""For specifying simulation parameters related to rectangular grid-like probe
scan patterns.

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



# For validating, pre-serializing, and de-pre-serializing integers used as seeds
# to random number generators.
import prismatique.thermal

# For determining the dimensions of sample supercells.
import prismatique.sample



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params"]



def _check_and_convert_step_size(params):
    obj_name = "step_size"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    step_size = czekitout.convert.to_pair_of_positive_floats(**kwargs)
    
    return step_size



def _pre_serialize_step_size(step_size):
    obj_to_pre_serialize = step_size
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_step_size(serializable_rep):
    step_size = serializable_rep

    return step_size



def _check_and_convert_window(params):
    obj_name = "window"

    current_func_name = "_check_and_convert_window"

    try:
        kwargs = {"obj": params[obj_name], "obj_name": obj_name}
        window = czekitout.convert.to_tuple_of_floats(**kwargs)

        if len(window) != 4:
            raise        
        if ((not (window[0] < window[1])) or (not (window[2] < window[3]))):
            raise
    except:
        err_msg = globals()[current_func_name+"_err_msg_1"]
        raise TypeError(err_msg)
    
    return window



def _pre_serialize_window(window):
    obj_to_pre_serialize = window
    module_alias = prismatique.tilt
    func_alias = module_alias._pre_serialize_window
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_window(serializable_rep):
    module_alias = prismatique.tilt
    func_alias = module_alias._de_pre_serialize_window
    window = func_alias(serializable_rep)

    return window



def _check_and_convert_jitter(params):
    obj_name = "jitter"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    jitter = czekitout.convert.to_nonnegative_float(**kwargs)
    
    return jitter



def _pre_serialize_jitter(jitter):
    obj_to_pre_serialize = jitter
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_jitter(serializable_rep):
    jitter = serializable_rep

    return jitter



def _check_and_convert_rng_seed(params):
    module_alias = prismatique.thermal
    func_alias = module_alias._check_and_convert_rng_seed
    rng_seed = func_alias(params)
    
    return rng_seed



def _pre_serialize_rng_seed(rng_seed):
    obj_to_pre_serialize = rng_seed
    module_alias = prismatique.thermal
    func_alias = module_alias._pre_serialize_rng_seed
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_rng_seed(serializable_rep):
    module_alias = prismatique.thermal
    func_alias = module_alias._de_pre_serialize_rng_seed
    rng_seed = func_alias(serializable_rep)

    return rng_seed



_module_alias = \
    prismatique.thermal
_default_step_size = \
    (0.25, 0.25)
_default_window = \
    (0, 1, 0, 1)
_default_jitter = \
    0
_default_rng_seed = \
    _module_alias._default_rng_seed
_default_skip_validation_and_conversion = \
    _module_alias._default_skip_validation_and_conversion



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related to rectangular grid-like probe scan
    patterns.

    A common scanning pattern in STEM is a rectangular grid. Sometimes it is
    useful to introduce small random deviations to such a scanning pattern. This
    is done in e.g. ptychography.

    The class :class:`prismatique.scan.rectangular.Params` specifies a scanning
    pattern that is obtained by generating an underlying rectangular grid to
    which to apply a random positional deviation to each point of the original
    grid, thus yielding an irregular scanning pattern that is rectangular
    grid-like, assuming the deviations are small. The radial components of the
    random positional deviations are sampled from a uniform distribution on an
    interval determined by the parameters defined below, and the polar angular
    components of the random position deviations are sampled from a uniform
    distribution on the interval :math:`\left[0, 2\pi\right)`.

    Both the
    :math:`x`- and :math:`y`-components of the random positional deviations are
    sampled from a normal distribution with zero-mean and a standard deviation
    which we refer to below as the "jitter" [the :math:`xy`-plane is assumed to
    coincide with the object plane].

    Parameters
    ----------
    step_size : `array_like` (`float`, shape=(``2``,)), optional
        ``step_size[0]`` specifies the probe scanning step size in the
        :math:`x`-direction in units of angstroms of the underlying rectangular
        grid; ``step_size[1]`` specifies the probe scanning step size in the
        :math:`y`-direction in units of angstroms of the underlying rectangular
        grid.
    window : `array_like` (`float`, shape=(``4``,)), optional
        ``window`` specifies a set of fractional coordinates that define the
        area of the underlying rectangular grid. The fractional coordinates are
        defined with respect to the spatial dimensions of the sample's supercell
        [see the documentation for the class
        :class:`prismatique.discretization.Params` for a discussion on sample
        supercells]. E.g. a fractional :math:`x`-coordinate of 0.5 corresponds
        to the midpoint of some supercell, along the :math:`x`-axis; a
        fractional :math:`x`-coordinate of 1.5 corresponds to the midpoint of
        the supercell to the right of the previous supercell, along the
        :math:`x`-axis; and a fractional :math:`x`-coordinate of -0.5
        corresponds to the midpoint of the supercell to the left of the first
        supercell mentioned, along the :math:`x`-axis.

        ``window[0]`` specifies the minimum fractional :math:`x`-coordinate of
        the scanning window of the underlying rectangular grid; ``window[1]``
        specifies the maximum fractional :math:`x`-coordinate of the scanning
        window of the underlying rectangular grid; ``window[2]`` specifies the
        minimum fractional :math:`y`-coordinate of the scanning window of the
        underlying rectangular grid; ``window[3]`` specifies the maximum
        fractional :math:`y`-coordinate of the scanning window of the underlying
        rectangular grid. Note that the tiling starts at the fractional
        coordinates ``(window[0], window[2])``.
    jitter : `int`, optional
        ``jitter``, along with the parameter ``step_size`` above, determine the
        interval from which the radial components of the random positional 
        deviations are sampled uniformly. Let :math:`\mathcal{J}`, 
        :math:`\Delta x`, and :math:`\Delta y` denote ``jitter``, 
        ``step_size[0]``, and ``step_size[1]`` respectively, then the interval 
        from which the radial components of the random positional deviations are
        sampled uniformly is 
        :math:`\left[0, \mathcal{J}\max\left(\Delta x, \Delta y\right)\right)`.
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
    ctor_param_names = ("step_size",
                        "window",
                        "jitter",
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
                 step_size=\
                 _default_step_size,
                 window=\
                 _default_window,
                 jitter=\
                 _default_jitter,
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



def _check_and_convert_rectangular_scan_params(params):
    obj_name = "rectangular_scan_params"
    obj = params[obj_name]

    accepted_types = (Params, type(None))
    
    if isinstance(obj, accepted_types[-1]):
        rectangular_scan_params = accepted_types[0]()
    else:
        kwargs = {"obj": obj,
                  "obj_name": obj_name,
                  "accepted_types": accepted_types}
        czekitout.check.if_instance_of_any_accepted_types(**kwargs)

        kwargs = obj.get_core_attrs(deep_copy=False)
        rectangular_scan_params = accepted_types[0](**kwargs)

    return rectangular_scan_params



def _pre_serialize_rectangular_scan_params(rectangular_scan_params):
    obj_to_pre_serialize = rectangular_scan_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_rectangular_scan_params(serializable_rep):
    rectangular_scan_params = Params.de_pre_serialize(serializable_rep)
    
    return rectangular_scan_params



def _generate_probe_positions(rectangular_scan_params, sample_specification):
    rectangular_scan_params_core_attrs = \
        rectangular_scan_params.get_core_attrs(deep_copy=False)

    step_size = rectangular_scan_params_core_attrs["step_size"]
    window = rectangular_scan_params_core_attrs["window"]
    jitter = rectangular_scan_params_core_attrs["jitter"]
    rng_seed = rectangular_scan_params_core_attrs["rng_seed"]

    Delta_X, Delta_Y, _ = \
        prismatique.sample._supercell_dims(sample_specification)

    min_x_probe_coord = Delta_X * window[0]
    max_x_probe_coord = Delta_X * window[1]
    min_y_probe_coord = Delta_Y * window[2]
    max_y_probe_coord = Delta_Y * window[3]

    x_probe_coord_step = step_size[0]
    y_probe_coord_step = step_size[1]
    tol = 1e-10

    x_coords = np.arange(min_x_probe_coord,
                         max_x_probe_coord+tol,
                         x_probe_coord_step)
    y_coords = np.arange(min_y_probe_coord,
                         max_y_probe_coord+tol,
                         y_probe_coord_step)

    np.random.seed(rng_seed)

    probe_positions = tuple()
    for x_coord in x_coords:
        for y_coord in y_coords:
            dr = np.random.uniform(0, jitter*max(step_size))
            theta = np.random.uniform(0, 2*np.pi)
            dx = dr * np.cos(theta)
            dy = dr * np.sin(theta)
            probe_positions += ((x_coord+dx, y_coord+dy),)

    np.random.seed()

    return probe_positions
    


###########################
## Define error messages ##
###########################

_check_and_convert_window_err_msg_1 = \
    ("The object ``window`` must be an array-like object of length 4, with "
     "each element being of type `float`, and the array as a whole satisfying: "
     "``window[0]<window[1]`` and ``window[2]<window[3]``.")
