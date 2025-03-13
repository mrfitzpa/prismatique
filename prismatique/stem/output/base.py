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
r"""For specifying the base output parameters for STEM simulations.

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
import prismatique.hrtem.output

# For validating, pre-serializing, and de-pre-serializing instances of the
# class :class:`prismatique.cbed.Params`.
import prismatique.cbed



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params"]



def _check_and_convert_output_dirname(params):
    module_alias = prismatique.hrtem.output
    func_alias = module_alias._check_and_convert_output_dirname
    output_dirname = func_alias(params)
    
    return output_dirname



def _pre_serialize_output_dirname(output_dirname):
    obj_to_pre_serialize = output_dirname
    module_alias = prismatique.hrtem.output
    func_alias = module_alias._pre_serialize_output_dirname
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_output_dirname(serializable_rep):
    module_alias = prismatique.hrtem.output
    func_alias = module_alias._de_pre_serialize_output_dirname
    output_dirname = func_alias(serializable_rep)

    return output_dirname



def _check_and_convert_max_data_size(params):
    module_alias = prismatique.hrtem.output
    func_alias = module_alias._check_and_convert_max_data_size
    max_data_size = func_alias(params)
    
    return max_data_size



def _pre_serialize_max_data_size(max_data_size):
    obj_to_pre_serialize = max_data_size
    module_alias = prismatique.hrtem.output
    func_alias = module_alias._pre_serialize_max_data_size
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_max_data_size(serializable_rep):
    module_alias = prismatique.hrtem.output
    func_alias = module_alias._de_pre_serialize_max_data_size
    max_data_size = func_alias(serializable_rep)

    return max_data_size



def _check_and_convert_cbed_params(params):
    module_alias = prismatique.cbed
    func_alias = module_alias._check_and_convert_cbed_params
    cbed_params = func_alias(params)
    
    return cbed_params



def _pre_serialize_cbed_params(cbed_params):
    obj_to_pre_serialize = cbed_params
    module_alias = prismatique.cbed
    func_alias = module_alias._pre_serialize_cbed_params
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_cbed_params(serializable_rep):
    module_alias = prismatique.cbed
    func_alias = module_alias._de_pre_serialize_cbed_params
    cbed_params = func_alias(serializable_rep)

    return cbed_params



def _check_and_convert_radial_step_size_for_3d_stem(params):
    obj_name = \
        "radial_step_size_for_3d_stem"
    kwargs = \
        {"obj": params[obj_name], "obj_name": obj_name}
    radial_step_size_for_3d_stem = \
        czekitout.convert.to_nonnegative_float(**kwargs)
    
    return radial_step_size_for_3d_stem



def _pre_serialize_radial_step_size_for_3d_stem(radial_step_size_for_3d_stem):
    obj_to_pre_serialize = radial_step_size_for_3d_stem
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_radial_step_size_for_3d_stem(serializable_rep):
    radial_step_size_for_3d_stem = serializable_rep

    return radial_step_size_for_3d_stem



def _check_and_convert_radial_range_for_2d_stem(params):
    obj_name = "radial_range_for_2d_stem"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    radial_range_for_2d_stem = czekitout.convert.to_pair_of_floats(**kwargs)

    current_func_name = "_check_and_convert_radial_range_for_2d_stem"

    if not (0 <= radial_range_for_2d_stem[0] <= radial_range_for_2d_stem[1]):
        err_msg = globals()[current_func_name+"_err_msg_1"]
        raise TypeError(err_msg)
    
    return radial_range_for_2d_stem



def _pre_serialize_radial_range_for_2d_stem(radial_range_for_2d_stem):
    obj_to_pre_serialize = radial_range_for_2d_stem
    serializable_rep = obj_to_pre_serialize
    
    return serializable_rep



def _de_pre_serialize_radial_range_for_2d_stem(serializable_rep):
    radial_range_for_2d_stem = serializable_rep
    
    return radial_range_for_2d_stem



def _check_and_convert_save_com(params):
    obj_name = "save_com"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    save_com = czekitout.convert.to_bool(**kwargs)
    
    return save_com



def _pre_serialize_save_com(save_com):
    obj_to_pre_serialize = save_com
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_save_com(serializable_rep):
    save_com = serializable_rep

    return save_com



def _check_and_convert_save_potential_slices(params):
    module_alias = prismatique.hrtem.output
    func_alias = module_alias._check_and_convert_save_potential_slices
    save_potential_slices = func_alias(params)
    
    return save_potential_slices



def _pre_serialize_save_potential_slices(save_potential_slices):
    obj_to_pre_serialize = save_potential_slices
    module_alias = prismatique.hrtem.output
    func_alias = module_alias._pre_serialize_save_potential_slices
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_save_potential_slices(serializable_rep):
    module_alias = prismatique.hrtem.output
    func_alias = module_alias._de_pre_serialize_save_potential_slices
    save_potential_slices = func_alias(serializable_rep)

    return save_potential_slices



_module_alias_1 = \
    prismatique.hrtem.output
_module_alias_2 = \
    prismatique.cbed
_default_output_dirname = \
    _module_alias_1._default_output_dirname
_default_max_data_size = \
    _module_alias_1._default_max_data_size
_default_cbed_params = \
    _module_alias_2._default_cbed_params
_default_radial_step_size_for_3d_stem = \
    1
_default_radial_range_for_2d_stem = \
    (0, 0)
_default_save_com = \
    False
_default_save_potential_slices = \
    _module_alias_1._default_save_potential_slices
_default_skip_validation_and_conversion = \
    _module_alias_1._default_skip_validation_and_conversion



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The base output parameters for STEM simulations.

    For a general discussion on the possible output that can be generated from
    STEM simulations, see the documentation for the module
    :mod:`prismatique.stem.output`. That discussion provides important context
    to the description of the parameter set below.

    Parameters
    ----------
    output_dirname : `str`, optional
        The relative or absolute path to the directory in which all output files
        are to be saved. If the directory doesn't exist upon saving the output 
        files, it will be created if possible.
    max_data_size : `int`, optional
        The data size limit, in bytes, of the STEM simulation output to be
        generated. If the output to be generated would require a data size
        larger than the aforementioned limit, then an exception is raised and
        the STEM simulation is not performed. Note that data size due to HDF5
        file overhead and metadata are not taken into account.
    cbed_params : :class:`prismatique.cbed.Params` | `None`, optional
        The simulation parameters related to convengent beam electron
        diffraction patterns, which includes parameters specifying what kind of
        4D-STEM output should be saved, if any at all. If ``cbed_params`` is set
        to `None` [i.e. the default value], then the aforementioned parameters
        are set to default values.
    radial_step_size_for_3d_stem : `float`, optional
        The bin width in mrads of the annular detectors used in 3D-STEM data
        collection. If set to zero, then no 3D-STEM data is saved. Note that 
        ``radial_step_size_for_3d_stem`` must be nonnegative.
    radial_range_for_2d_stem : `array_like` (`float`, shape=(``2``,)), optional
        ``radial_range_for_2d_stem[0]`` and ``radial_range_for_2d_stem[1]`` are
        the lower and upper radial integration limits respectively in
        angular/Fourier space in units of mrads used to obtain the 2D-STEM
        intensity data. If
        ``0<=radial_range_for_2d_stem[0]<radial_range_for_2d_stem[1]``, then the
        intensity 2D-STEM data is written to a file with the basename
        ``"stem_sim_intensity_output.h5"``. If
        ``0<=radial_range_for_2d_stem[0]==radial_range_for_2d_stem[1]``, then no
        2D-STEM data is saved. In all other scenarios an error is raised.
    save_com : `bool`, optional
        If ``save_com`` is set to ``True``, then the center-of-mass (COM)
        momentum averaged over atomic configurations is written to a file with
        the basename ``"stem_sim_intensity_output.h5"``. If ``save_com`` is set 
        to ``False``, then no COM data is saved.
    save_potential_slices : `bool`, optional
        If ``save_potential_slices`` is set to ``True``, then for each frozen
        phonon configuration subset, the corresponding potential slice data is
        written to a file with the basename
        ``"potential_slices_of_subset_"+str(i)+".h5"``, with ``i`` being the 
        subset index. If ``save_potential_slices`` is set to ``False``, then no
        potential slice data is saved.
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
    ctor_param_names = ("output_dirname",
                        "max_data_size",
                        "cbed_params",
                        "radial_step_size_for_3d_stem",
                        "radial_range_for_2d_stem",
                        "save_com",
                        "save_potential_slices")
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
                 output_dirname=\
                 _default_output_dirname,
                 max_data_size=\
                 _default_max_data_size,
                 cbed_params=\
                 _default_cbed_params,
                 radial_step_size_for_3d_stem=\
                 _default_radial_step_size_for_3d_stem,
                 radial_range_for_2d_stem=\
                 _default_radial_range_for_2d_stem,
                 save_com=\
                 _default_save_com,
                 save_potential_slices=\
                 _default_save_potential_slices,
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



def _check_and_convert_base_params(params):
    obj_name = "base_params"
    obj = params[obj_name]

    accepted_types = (Params, type(None))
    
    if isinstance(obj, accepted_types[-1]):
        base_params = accepted_types[0]()
    else:
        kwargs = {"obj": obj,
                  "obj_name": obj_name,
                  "accepted_types": accepted_types}
        czekitout.check.if_instance_of_any_accepted_types(**kwargs)

        kwargs = obj.get_core_attrs(deep_copy=False)
        base_params = accepted_types[0](**kwargs)

    return base_params



def _pre_serialize_base_params(base_params):
    obj_to_pre_serialize = base_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_base_params(serializable_rep):
    base_params = Params.de_pre_serialize(serializable_rep)
    
    return base_params



_default_base_params = None



###########################
## Define error messages ##
###########################

_check_and_convert_radial_range_for_2d_stem_err_msg_1 = \
    ("The object ``radial_range_for_2d_stem`` must either be a pair of "
     "real numbers satisfying: "
     "``0<=radial_range_for_2d_stem[0]<=radial_range_for_2d_stem[1]``.")
