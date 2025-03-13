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
r"""For specifying output parameters that are applicable only to the PRISM
implementation of STEM simulations.

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
import prismatique.sample



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params"]



def _check_and_convert_enable_S_matrix_refocus(params):
    obj_name = "enable_S_matrix_refocus"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    enable_S_matrix_refocus = czekitout.convert.to_bool(**kwargs)
    
    return enable_S_matrix_refocus



def _pre_serialize_enable_S_matrix_refocus(enable_S_matrix_refocus):
    obj_to_pre_serialize = enable_S_matrix_refocus
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_enable_S_matrix_refocus(serializable_rep):
    enable_S_matrix_refocus = serializable_rep

    return enable_S_matrix_refocus



def _check_and_convert_save_S_matrices(params):
    obj_name = "save_S_matrices"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    save_S_matrices = czekitout.convert.to_bool(**kwargs)
    
    return save_S_matrices



def _pre_serialize_save_S_matrices(save_S_matrices):
    obj_to_pre_serialize = save_S_matrices
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_save_S_matrices(serializable_rep):
    save_S_matrices = serializable_rep

    return save_S_matrices



_module_alias = \
    prismatique.sample
_default_enable_S_matrix_refocus = \
    False
_default_save_S_matrices = \
    False
_default_skip_validation_and_conversion = \
    _module_alias._default_skip_validation_and_conversion



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The output parameters that are applicable only to the PRISM
    implementation of STEM simulations.

    For a general discussion on the possible output that can be generated from
    STEM simulations, see the documentation for the class
    :class:`prismatique.stem.output.Params`. For a general discussion on the
    PRISM algorithm and :math:`S`-matrices, see the documentation for the module
    :mod:`prismatique.stem`. These discussions provide important context to the
    description of the parameter set below.

    Parameters
    ----------
    enable_S_matrix_refocus : `bool`, optional
        If ``enable_S_matrix_refocus`` is set to ``True``, then :math:`S`-matrix
        refocusing is enabled. If ``enable_S_matrix_refocus`` is set to
        ``False``, then no :math:`S`-matrix refocusing is performed. 
    save_S_matrices : `bool`, optional
        If ``save_S_matrices`` is set to ``True``, then for each frozen phonon
        configuration subset, the corresponding :math:`S`-matrix data is written
        to the file with the basename ``"S_matrices_of_subset_"+str(i)+".h5"``,
        with ``i`` being the subset index. If ``save_S_matrices`` is set to
        ``False``, then no :math:`S`-matrix data is saved.
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
    ctor_param_names = ("enable_S_matrix_refocus",
                        "save_S_matrices")
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
                 enable_S_matrix_refocus=\
                 _default_enable_S_matrix_refocus,
                 save_S_matrices=\
                 _default_save_S_matrices,
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



def _check_and_convert_prism_output_params(params):
    obj_name = "prism_output_params"
    obj = params[obj_name]

    accepted_types = (Params,)
    
    kwargs = {"obj": obj,
              "obj_name": obj_name,
              "accepted_types": accepted_types}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    kwargs = obj.get_core_attrs(deep_copy=False)
    prism_output_params = accepted_types[0](**kwargs)

    return prism_output_params



def _pre_serialize_prism_output_params(prism_output_params):
    obj_to_pre_serialize = prism_output_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_prism_output_params(serializable_rep):
    prism_output_params = Params.de_pre_serialize(serializable_rep)
    
    return prism_output_params



_default_prism_output_params = None
    


###########################
## Define error messages ##
###########################
