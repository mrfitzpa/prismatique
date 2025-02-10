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
r"""For specifying simulation parameters related to the objective aperture.

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
import prismatique.tilt



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params"]



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

        if len(window) != 2:
            raise
        
        if not (window[0] <= window[1]):
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



_module_alias = \
    prismatique.tilt
_default_offset = \
    _module_alias._default_offset
_default_window = \
    (0, float("inf"))
_default_skip_validation_and_conversion = \
    _module_alias._default_skip_validation_and_conversion



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related objective aperture.

    The objective aperture is assumed to be annular, and is used in HRTEM
    simulations.

    The parameters below define an angular window within which all scattering
    angles that are not blocked by the objective aperture.

    Parameters
    ----------
    offset : array_like` (`float`, shape=(``2``,)), optional
        ``offset`` specifies the offset of the angular window of the objective
        aperture discussed above: ``offset[0]`` specifies the
        :math:`x`-coordinate of the offset in mrads; ``offset[1]`` specifies the
        :math:`y`-coordinate of the offset in mrads. 
    window : array_like` (`float`, shape=(``2``,)), optional
        If ``window`` is an array of length 2, then ``window`` specifies a
        radial window: ``window[0]`` and ``window[1]`` specify the minimum and
        maximum radial angles with respect to the offset angle in mrads.
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
                        "window")
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



def _check_and_convert_objective_aperture_params(params):
    obj_name = "objective_aperture_params"
    obj = params[obj_name]

    accepted_types = (Params, type(None))
    
    if isinstance(obj, accepted_types[-1]):
        objective_aperture_params = accepted_types[0]()
    else:
        kwargs = {"obj": obj,
                  "obj_name": obj_name,
                  "accepted_types": accepted_types}
        czekitout.check.if_instance_of_any_accepted_types(**kwargs)

        kwargs = obj.get_core_attrs(deep_copy=False)
        objective_aperture_params = accepted_types[0](**kwargs)

    return objective_aperture_params



def _pre_serialize_objective_aperture_params(objective_aperture_params):
    obj_to_pre_serialize = objective_aperture_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_objective_aperture_params(serializable_rep):
    objective_aperture_params = Params.de_pre_serialize(serializable_rep)
    
    return objective_aperture_params



_default_objective_aperture_params = None
    


###########################
## Define error messages ##
###########################

_check_and_convert_window_err_msg_1 = \
    ("The object ``window`` must be a pair of real numbers satisfying: "
     "``0<=window[0]<=window[1]``.")
