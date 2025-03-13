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
r"""For specifying output parameters that are applicable only to the multislice
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



def _check_and_convert_num_slices_per_output(params):
    obj_name = "num_slices_per_output"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    num_slices_per_output = czekitout.convert.to_positive_int(**kwargs)
    
    return num_slices_per_output



def _pre_serialize_num_slices_per_output(num_slices_per_output):
    obj_to_pre_serialize = num_slices_per_output
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_num_slices_per_output(serializable_rep):
    num_slices_per_output = serializable_rep

    return num_slices_per_output



def _check_and_convert_z_start_output(params):
    obj_name = "z_start_output"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    z_start_output = czekitout.convert.to_nonnegative_float(**kwargs)
    
    return z_start_output



def _pre_serialize_z_start_output(z_start_output):
    obj_to_pre_serialize = z_start_output
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_z_start_output(serializable_rep):
    z_start_output = serializable_rep

    return z_start_output



_module_alias = \
    prismatique.sample
_default_num_slices_per_output = \
    1
_default_z_start_output = \
    float("inf")
_default_skip_validation_and_conversion = \
    _module_alias._default_skip_validation_and_conversion



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The output parameters that are applicable only to the multislice
    implementation of STEM simulations.

    For a general discussion on the possible output that can be generated from
    STEM simulations, see the documentation for the class
    :class:`prismatique.stem.output.Params`. That discussion provides important
    context to the description of the parameter set below.

    Parameters
    ----------
    num_slices_per_output : `int`, optional
        In the multislice algorithm, the sample is partitioned into slices along
        the :math:`z`-axis [i.e. the optic axis]. One can optionally save output
        for exit waves that emerge immediately after select intermediate slices,
        in addition to the output for the final exit wave. The output generated
        after a given slice is saved into an "output layer", i.e. a given output
        layer will contain the output for an exit wave that emerged from an
        intermediate slice or the final slice. Note that the exit wave emerging
        from the :math:`n^{\text{th}}` slice, where :math:`n=0` is the first
        slice, emerges at a depth of :math:`\left(n+1\right)\delta z`, where
        :math:`\delta z` is the slice thickness, given by
        Eq. :eq:`slice_thickness_in_potential_params`. ``num_slices_per_output``
        specifies the number of slices between intermediate outputs. For
        example, if ``num_slices_per_output`` is set to ``1``, then the output
        from each slice is saved; if ``num_slices_per_output`` is set to ``2``,
        then the output from every second slice is saved; and so on. Note that
        the output for the exit wave that emerges from the final slice is always
        saved to the last output layer, irrespective of the value of
        ``num_slices_per_output``. Moreover, note that ``num_slices_per_output``
        must be a positive `int`.
    z_start_output : `float`, optional
        Continuing from above, ``z_start_output`` specifies the depth in
        angstroms along the :math:`z`-axis at which intermediate output
        collection begins. Let ``m = int(math.ceil(z_start_output /
        (num_slices_per_output*slice_thickness))`` and ``n = max(1,
        m)*slices_per_output - 1``, where ``slice_thickness`` is the slice
        thickness :math:`\delta z` in angstroms. If ``n<=N_slice-1``, where
        ``N_slice`` is the total number of slices used to partition the sample,
        then intermediate output is stored with the output from the exit wave
        that emerges from the ``n`` th slice being saved into the first output
        layer. Note that the output for the exit wave that emerges from the
        final slice is always saved to the last output layer, irrespective of
        the value of ``z_start_output``. Note that ``z_start_output`` must be a
        nonnegative `float`.
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
    ctor_param_names = ("num_slices_per_output",
                        "z_start_output")
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
                 num_slices_per_output=\
                 _default_num_slices_per_output,
                 z_start_output=\
                 _default_z_start_output,
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



def _check_and_convert_multislice_output_params(params):
    obj_name = "multislice_output_params"
    obj = params[obj_name]

    accepted_types = (Params, type(None))
    
    if isinstance(obj, accepted_types[-1]):
        multislice_output_params = accepted_types[0]()
    else:
        kwargs = {"obj": obj,
                  "obj_name": obj_name,
                  "accepted_types": accepted_types}
        czekitout.check.if_instance_of_any_accepted_types(**kwargs)

        kwargs = obj.get_core_attrs(deep_copy=False)
        multislice_output_params = accepted_types[0](**kwargs)

    return multislice_output_params



def _pre_serialize_multislice_output_params(multislice_output_params):
    obj_to_pre_serialize = multislice_output_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_multislice_output_params(serializable_rep):
    multislice_output_params = Params.de_pre_serialize(serializable_rep)
    
    return multislice_output_params



_default_multislice_output_params = None
    


###########################
## Define error messages ##
###########################
