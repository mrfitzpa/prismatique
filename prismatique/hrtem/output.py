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
r"""For specifying the output parameters for HRTEM simulations.

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



# For validating instances of the class
# :class:`prismatique.hrtem.system.ModelParams`.
import prismatique.hrtem.system

# For validating instances of the classes
# :class:`prismatique.sample.ModelParams`, and
# :class:`prismatique.sample.PotentialSliceSubsetIDs`; and for calculating
# quantities related to the modelling of the sample.
import prismatique.sample

# For validating, pre-serializing, and de-pre-serializing instances of the class
# :class:`prismatique.hrtem.image.Params`.
import prismatique.hrtem.image

# For postprocessing HRTEM intensity images.
import prismatique._signal

# For generating tilts used in HRTEM simulations.
import prismatique.tilt



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params",
           "data_size"]



def _check_and_convert_output_dirname(params):
    module_alias = prismatique.sample
    func_alias = module_alias._check_and_convert_output_dirname
    output_dirname = func_alias(params)
    
    return output_dirname



def _pre_serialize_output_dirname(output_dirname):
    obj_to_pre_serialize = output_dirname
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_output_dirname(serializable_rep):
    output_dirname = serializable_rep

    return output_dirname



def _check_and_convert_max_data_size(params):
    module_alias = prismatique.sample
    func_alias = module_alias._check_and_convert_max_data_size
    max_data_size = func_alias(params)
    
    return max_data_size



def _pre_serialize_max_data_size(max_data_size):
    obj_to_pre_serialize = max_data_size
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_max_data_size(serializable_rep):
    max_data_size = serializable_rep

    return max_data_size



def _check_and_convert_image_params(params):
    module_alias = prismatique.hrtem.image
    func_alias = module_alias._check_and_convert_image_params
    image_params = func_alias(params)

    return image_params



def _pre_serialize_image_params(image_params):
    obj_to_pre_serialize = image_params
    module_alias = prismatique.hrtem.image
    func_alias = module_alias._pre_serialize_image_params
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_image_params(serializable_rep):
    module_alias = prismatique.hrtem.image
    func_alias = module_alias._de_pre_serialize_image_params
    image_params = func_alias(serializable_rep)

    return image_params



def _check_and_convert_save_potential_slices(params):
    obj_name = "save_potential_slices"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    save_potential_slices = czekitout.convert.to_bool(**kwargs)
    
    return save_potential_slices



def _pre_serialize_save_potential_slices(save_potential_slices):
    obj_to_pre_serialize = save_potential_slices
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_save_potential_slices(serializable_rep):
    save_potential_slices = serializable_rep

    return save_potential_slices



_module_alias_1 = \
    prismatique.sample
_default_output_dirname = \
    "sim_output_files"
_default_max_data_size = \
    _module_alias_1._default_max_data_size
_default_image_params = \
    None
_default_save_potential_slices = \
    False
_default_skip_validation_and_conversion = \
    _module_alias_1._default_skip_validation_and_conversion



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The output parameters for HRTEM simulations.

    See the documentation for the class :class:`prismatique.hrtem.image.Params`
    for a discussion on HRTEM image wavefunctions and intensities that is
    relevant to the discussion on this page.  

    Upon the completion of a HRTEM simulation, ``prismatique`` can optionally
    save a variety of different HRTEM data to a set of HDF5 files based on the
    specifications of the user. This is done by taking the original output file
    generated by ``prismatic``, and then restructuring it into one or more
    output files. If the output parameters specify that intensity data be saved,
    then said data is extracted from the original ``prismatic`` output files,
    postprocessed, and written to a new file with basename
    ``"hrtem_sim_intensity_output.h5"`` in a more readable layout. This new
    output file has the following structure:

    - metadata: <HDF5 group>
    
      * r_x: <HDF5 1D dataset>
        
        + dim 1: "r_x idx"
        + units: "Å"

      * r_y: <HDF5 1D dataset>
        
        + dim 1: "r_y idx"
        + units: "Å"

    - data: <HDF5 group>

      * intensity_image: <HDF5 2D dataset>
        
        + dim 1: "r_y idx"
        + dim 2: "r_x idx"
        + units: "dimensionless"

    Note that the sub-bullet points listed immediately below a given HDF5
    dataset display the HDF5 attributes associated with said HDF5 dataset. Each
    HDF5 scalar and dataset has a ``"units"`` attribute which, as the name
    suggests, indicates the units in which said data [i.e. the scalar or
    dataset] is expressed. Each HDF5 dataset will also have a set of attributes
    with names of the form ``"dim {}".format(i)`` with ``i`` being an integer
    ranging from 1 to the rank of said HDF5 dataset. Attribute ``"dim
    {}".format(i)`` of a given HDF5 dataset labels the :math:`i^{\text{th}}`
    dimension of the underlying array of the dataset. Most of these dimension
    labels should be self-explanatory but for clarification: "idx" is short for
    "index"; "avg" is short for "average"; and "r_x" and "r_y" refer to the
    :math:`x`- and :math:`y`-coordinates in the discretized real-space.

    If the output parameters specify that complex-valued wavefunction data be
    saved, then said data is extracted from the original ``prismatic`` output
    files, and written to a new set of files: one file per frozen phonon/atomic
    configuration subset. Note that, unlike the intensity data, the
    complex-valued wavefunction data is not postprocessed. See the documentation
    for the class :class:`prismatique.thermal.Params` for a discussion on frozen
    phonon configurations and their grouping into subsets. For the ``i`` th
    frozen phonon configuration subset, the corresponding wavefunction data is
    saved to the file with basename
    ``"hrtem_sim_wavefunction_output_of_subset_"+str(i)+".h5"``. Each one of
    these new output files has the following structure:

    - metadata: <HDF5 group>
    
      * tilts: <HDF5 2D dataset>
        
        + dim 1: "tilt idx"
        + dim 2: "vector component idx [0->x, 1->y]"
        + units: "mrad"

      * defocii: <HDF5 1D dataset>
        
        + dim 1: "defocus idx"
        + units: "Å"

      * r_x: <HDF5 1D dataset>
        
        + dim 1: "r_x idx"
        + units: "Å"

      * r_y: <HDF5 1D dataset>
        
        + dim 1: "r_y idx"
        + units: "Å"

    - data: <HDF5 group>

      * image_wavefunctions: <HDF5 5D dataset>
        
        + dim 1: "atomic config idx"
        + dim 2: "defocus idx"
        + dim 3: "tilt idx"
        + dim 4: "r_y idx"
        + dim 5: "r_x idx"
        + units: "dimensionless"

    ``prismatique`` can also optionally save the "potential slice" [i.e.
    Eq. :eq:`coarse_grained_potential_1`] data for each subset of frozen phonon
    configurations into separate HDF5 output files by extracting said data from
    the original ``prismatic`` output files. See the documentation for the class
    :class:`prismatique.thermal.Params` for a discussion on frozen phonon
    configurations their grouping into subsets. For the ``i`` th subset, the
    corresponding potential slice data is saved to an HDF5 output file with
    basename ``"potential_slices_"+str(i)+".h5"``. Unlike the output data in the
    file ``"hrtem_simulation_output.h5"``, the layout of the potential slice
    data is kept the same as that found in the original file [i.e. the same HDF5
    paths are used].  The same layout needs to be used in order for
    ``prismatic`` to be able to successfully import/load pre-calculated
    potential slice data for a future simulation.

    It is beyond the scope of the documentation to describe the structure of the
    potential slice output files. Users can analyze the data in these output
    files with the help of the tools found in the module
    :mod:`prismatique.load`.

    The last file that is always generated after running a HRTEM simulation is a
    JSON file that contains, in a serialized format, the simulation parameters
    used. This file has the basename ``"hrtem_simulation_parameters.json"``.

    Parameters
    ----------
    output_dirname : `str`, optional
        The relative or absolute path to the directory in which all output files
        are to be saved. If the directory doesn't exist upon saving the output 
        files, it will be created if possible.
    max_data_size : `int`, optional
        The data size limit, in bytes, of the HRTEM simulation output to be
        generated. If the output to be generated would require a data size
        larger than the aforementioned limit, then an exception is raised and
        the HRTEM simulation is not performed. Note that data size due to HDF5
        file overhead and metadata are not taken into account.
    image_params : :class:`prismatique.hrtem.image.Params` | `None`, optional
        The simulation parameters related to image wavefunctions and
        intensities, which includes parameters specifying what kind of image
        output should be saved, if any at all. If ``image_params`` is set to
        `None` [i.e. the default value], then the aforementioned parameters are
        set to default values.
    save_potential_slices : `bool`, optional
        If ``save_potential_slices`` is set to ``True``, then for each frozen
        phonon configuration subset, the corresponding potential slice data is
        written to the file with the basename
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
                        "image_params",
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
                 image_params=\
                 _default_image_params,
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



def _check_and_convert_output_params(params):
    obj_name = "output_params"
    obj = params[obj_name]

    accepted_types = (Params, type(None))
    
    if isinstance(obj, accepted_types[-1]):
        output_params = accepted_types[0]()
    else:
        kwargs = {"obj": obj,
                  "obj_name": obj_name,
                  "accepted_types": accepted_types}
        czekitout.check.if_instance_of_any_accepted_types(**kwargs)

        kwargs = obj.get_core_attrs(deep_copy=False)
        output_params = accepted_types[0](**kwargs)

    return output_params



def _pre_serialize_output_params(output_params):
    obj_to_pre_serialize = output_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_output_params(serializable_rep):
    output_params = Params.de_pre_serialize(serializable_rep)
    
    return output_params



_default_output_params = None



def _check_and_convert_hrtem_system_model_params(params):
    module_alias = prismatique.hrtem.system
    func_alias = module_alias._check_and_convert_hrtem_system_model_params
    hrtem_system_model_params = func_alias(params)

    return hrtem_system_model_params



def _check_and_convert_skip_validation_and_conversion(params):
    module_alias = prismatique.sample
    func_alias = module_alias._check_and_convert_skip_validation_and_conversion
    skip_validation_and_conversion = func_alias(params)

    return skip_validation_and_conversion



def data_size(hrtem_system_model_params,
              output_params=\
              _default_output_params,
              skip_validation_and_conversion=\
              _default_skip_validation_and_conversion):
    r"""Calculate the data size of the HRTEM simulation output that one could
    generate according to a given HRTEM system model, and output parameter set.

    Note that data size due to HDF5 file overhead and metadata are not taken
    into account.

    Parameters
    ----------
    hrtem_system_model_params : :class:`prismatique.hrtem.system.ModelParams`
        The simulation parameters related to the modelling of the HRTEM
        system. See the documentation for the class
        :class:`prismatique.hrtem.system.ModelParams` for a discussion on said
        parameters.
    output_params : :class:`prismatique.hrtem.output.Params` | `None`, optional
        The output parameters for the HRTEM simulation. See the documentation
        for the class :class:`prismatique.hrtem.output.Params` for a discussion
        on said parameters. If ``output_params`` is set to `None` [i.e. the
        default value], then the aforementioned simulation parameters are set to
        default values.
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
    output_data_size : `int`
        The data size in units of bytes.

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
    output_data_size  = _data_size(**kwargs)
    
    return None



def _data_size(hrtem_system_model_params, output_params):
    output_params_core_attrs = output_params.get_core_attrs(deep_copy=False)

    image_params = output_params_core_attrs["image_params"]
    image_params_core_attrs = image_params.get_core_attrs(deep_copy=False)

    output_data_size = 0

    kwargs = {"hrtem_system_model_params": hrtem_system_model_params,
              "output_params": output_params}
    if image_params_core_attrs["save_final_intensity"]:
        output_data_size += _data_size_of_intensity_output(**kwargs)
    if image_params_core_attrs["save_wavefunctions"]:
        del kwargs["output_params"]
        output_data_size += _data_size_of_wavefunction_output(**kwargs)

    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)

    sample_specification = \
        hrtem_system_model_params_core_attrs["sample_specification"]
    
    if output_params_core_attrs["save_potential_slices"]:
        kwargs = \
            {"sample_specification": sample_specification}
        output_data_size += \
            prismatique.sample._potential_slice_set_data_size(**kwargs)
            
    return output_data_size



def _data_size_of_intensity_output(hrtem_system_model_params, output_params):
    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)
    output_params_core_attrs = \
        output_params.get_core_attrs(deep_copy=False)

    sample_specification = \
        hrtem_system_model_params_core_attrs["sample_specification"]
    
    image_params = output_params_core_attrs["image_params"]
    image_params_core_attrs = image_params.get_core_attrs(deep_copy=False)
    
    postprocessing_seq = \
        image_params_core_attrs["postprocessing_seq"]
    num_pixels_in_postprocessed_2d_signal_space = \
        prismatique._signal._num_pixels_in_postprocessed_2d_signal_space

    kwargs = \
        {"sample_specification": sample_specification,
         "signal_is_cbed_pattern_set": False,
         "postprocessing_seq": postprocessing_seq}
    num_pixels_in_postprocessed_image = \
        num_pixels_in_postprocessed_2d_signal_space(**kwargs)

    size_of_single = 4  # In bytes.

    data_size_of_intensity_output = (size_of_single
                                     * num_pixels_in_postprocessed_image)

    return data_size_of_intensity_output



def _data_size_of_wavefunction_output(hrtem_system_model_params):
    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)

    sample_specification = \
        hrtem_system_model_params_core_attrs["sample_specification"]
    tilt_params = \
        hrtem_system_model_params_core_attrs["tilt_params"]
    gun_model_params = \
        hrtem_system_model_params_core_attrs["gun_model_params"]
    defocal_offset_supersampling = \
        hrtem_system_model_params_core_attrs["defocal_offset_supersampling"]

    gun_model_params_core_attrs = \
        gun_model_params.get_core_attrs(deep_copy=False)
    
    mean_beam_energy = \
        gun_model_params_core_attrs["mean_beam_energy"]

    kwargs = \
        {"sample_specification": sample_specification}
    total_num_frozen_phonon_configs = \
        prismatique.sample._total_num_frozen_phonon_configs(**kwargs)

    tilt_series = prismatique.tilt._series(sample_specification,
                                           mean_beam_energy,
                                           tilt_params)

    kwargs["signal_is_cbed_pattern_set"] = \
        False
    num_pixels_in_unprocessed_image = \
        prismatique._signal._num_pixels_in_unprocessed_2d_signal_space(**kwargs)

    size_of_single = 4  # In bytes.
    size_of_complex_single = 2 * size_of_single  # In bytes.

    data_size_of_wavefunction_output = (total_num_frozen_phonon_configs
                                        * defocal_offset_supersampling
                                        * len(tilt_series)
                                        * num_pixels_in_unprocessed_image
                                        * size_of_complex_single)

    return data_size_of_wavefunction_output



###########################
## Define error messages ##
###########################
