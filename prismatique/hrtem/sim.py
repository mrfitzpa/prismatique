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
r"""For running HRTEM simulations.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For performing operations on file and directory paths.
import pathlib

# For explicit garbage collection.
import gc



# For general array handling.
import numpy as np

# For validating and converting objects.
import czekitout.check
import czekitout.convert

# For defining classes that support enforced validation, updatability,
# pre-serialization, and de-serialization.
import fancytypes

# For loading objects from and saving objects to HDF5 files.
import h5pywrappers

# For conveniently extracting certain properties of the electron beam.
import embeam

# For postprocessing ``hyperspy`` signals.
import empix

# For creating a :obj:`pyprismatic.Metadata` object that is responsible for
# running the ``prismatic`` simulation.
import pyprismatic



# For validating, pre-serializing, and de-pre-serializing instances of the
# classes :class:`prismatique.worker.Params`,
# :class:`prismatique.hrtem.system.ModelParams`, and
# :class:`prismatique.hrtem.output.Params`.
import prismatique.worker
import prismatique.hrtem.system
import prismatique.hrtem.output

# For validating instances of the classes
# :class:`prismatique.sample.ModelParams`, and
# :class:`prismatique.sample.PotentialSliceSubsetIDs`; for calculating
# quantities related to the modelling of the sample; for validating certain
# filenames; and for importing various other helper functions.
import prismatique.sample

# For postprocessing HRTEM intensity images.
import prismatique._signal

# For generating tilts used in HRTEM simulations.
import prismatique.tilt



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params",
           "run"]



def _check_and_convert_hrtem_system_model_params(params):
    module_alias = prismatique.hrtem.system
    func_alias = module_alias._check_and_convert_hrtem_system_model_params
    hrtem_system_model_params = func_alias(params)

    return hrtem_system_model_params



def _pre_serialize_hrtem_system_model_params(hrtem_system_model_params):
    obj_to_pre_serialize = hrtem_system_model_params
    module_alias = prismatique.hrtem.system
    func_alias = module_alias._pre_serialize_hrtem_system_model_params
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_hrtem_system_model_params(serializable_rep):
    module_alias = prismatique.hrtem.system
    func_alias = module_alias._de_pre_serialize_hrtem_system_model_params
    hrtem_system_model_params = func_alias(serializable_rep)

    return hrtem_system_model_params



def _check_and_convert_output_params(params):
    module_alias = prismatique.hrtem.output
    func_alias = module_alias._check_and_convert_output_params
    output_params = func_alias(params)

    return output_params



def _pre_serialize_output_params(output_params):
    obj_to_pre_serialize = output_params
    module_alias = prismatique.hrtem.output
    func_alias = module_alias._pre_serialize_output_params
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_output_params(serializable_rep):
    module_alias = prismatique.hrtem.output
    func_alias = module_alias._de_pre_serialize_output_params
    output_params = func_alias(serializable_rep)

    return output_params



def _check_and_convert_worker_params(params):
    module_alias = prismatique.worker
    func_alias = module_alias._check_and_convert_worker_params
    worker_params = func_alias(params)

    return worker_params



def _pre_serialize_worker_params(worker_params):
    obj_to_pre_serialize = worker_params
    module_alias = prismatique.worker
    func_alias = module_alias._pre_serialize_worker_params
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_worker_params(serializable_rep):
    module_alias = prismatique.worker
    func_alias = module_alias._de_pre_serialize_worker_params
    worker_params = func_alias(serializable_rep)

    return worker_params



_module_alias = \
    prismatique.worker
_default_output_params = \
    None
_default_worker_params = \
    _module_alias._default_worker_params
_default_skip_validation_and_conversion = \
    _module_alias._default_skip_validation_and_conversion



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The HRTEM simulation parameters.

    Parameters
    ----------
    hrtem_system_model_params : :class:`prismatique.hrtem.system.ModelParams`
        The simulation parameters related to the modelling of HRTEM systems. See
        the documentation for the class
        :class:`prismatique.hrtem.system.ModelParams` for a discussion on said
        parameters.  
    output_params : :class:`prismatique.hrtem.output.Params` | `None`, optional
        The output parameters for the HRTEM simulation. See the documentation
        for the class :class:`prismatique.hrtem.output.Params` for a discussion
        on said parameters. If ``output_params`` is set to `None` [i.e. the
        default value], then the aforementioned simulation parameters are set to
        default values.
    worker_params : :class:`prismatique.worker.Params` | `None`, optional
        The simulation parameters related to GPU and CPU workers. See the
        documentation for the class :class:`prismatique.worker.Params` for a
        discussion on said parameters. If ``worker_params`` is set to `None`
        [i.e. the default value], then the aforementioned simulation parameters 
        are set to default values.
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
    ctor_param_names = ("hrtem_system_model_params",
                        "output_params",
                        "worker_params")
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
                 hrtem_system_model_params,
                 output_params=\
                 _default_output_params,
                 worker_params=\
                 _default_worker_params,
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



def _check_and_convert_sim_params(params):
    obj_name = "sim_params"
    obj = params[obj_name]

    accepted_types = (Params,)

    kwargs = {"obj": obj,
              "obj_name": obj_name,
              "accepted_types": accepted_types}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    kwargs = obj.get_core_attrs(deep_copy=False)
    sim_params = accepted_types[0](**kwargs)

    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]

    kwargs = {"params": dict()}
    kwargs["params"]["hrtem_system_model_params"] = hrtem_system_model_params
    _check_and_convert_hrtem_system_model_params(**kwargs)

    _pre_save(sim_params)
    _check_data_size(sim_params)

    return sim_params



def _pre_save(sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]
    output_params = \
        sim_params_core_attrs["output_params"]

    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)
    
    sample_specification = \
        hrtem_system_model_params_core_attrs["sample_specification"]

    output_params_core_attrs = \
        output_params.get_core_attrs(deep_copy=False)
    
    if output_params_core_attrs["save_potential_slices"]:
        kwargs = {"sample_specification": sample_specification,
                  "output_dirname": output_params_core_attrs["output_dirname"],
                  "unformatted_basename": "potential_slices_of_subset_{}.h5"}
        prismatique.sample._pre_save(**kwargs)

    filenames = tuple()
    if _intensity_output_is_to_be_saved(sim_params):
        filenames += (_intensity_output_filename(sim_params),)
    if _wavefunction_output_is_to_be_saved(sim_params):
        filenames += _wavefunction_output_filenames(sim_params)

    for filename in filenames:
        if pathlib.Path(filename).is_file():
            pathlib.Path(filename).unlink(missing_ok=True)

    prismatique.sample._check_hdf5_filenames(filenames)

    return None



def _intensity_output_is_to_be_saved(sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    output_params = \
        sim_params_core_attrs["output_params"]

    output_params_core_attrs = \
        output_params.get_core_attrs(deep_copy=False)

    image_params = \
        output_params_core_attrs["image_params"]
    
    image_params_core_attrs = \
        image_params.get_core_attrs(deep_copy=False)

    intensity_output_is_to_be_saved = \
        image_params_core_attrs["save_final_intensity"]

    return intensity_output_is_to_be_saved



def _intensity_output_filename(sim_params):
    output_param_subset = _output_param_subset(sim_params)
    output_dirname = output_param_subset["output_dirname"]
    filename = output_dirname + "/hrtem_sim_intensity_output.h5"

    return filename



def _output_param_subset(sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    output_params = \
        sim_params_core_attrs["output_params"]

    output_params_core_attrs = \
        output_params.get_core_attrs(deep_copy=False)

    output_dirname = output_params_core_attrs["output_dirname"]
    save_potential_slices = output_params_core_attrs["save_potential_slices"]

    output_param_subset = {"output_dirname": output_dirname,
                           "save_potential_slices": save_potential_slices,
                           "save_S_matrices": False}

    return output_param_subset



def _wavefunction_output_is_to_be_saved(sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    output_params = \
        sim_params_core_attrs["output_params"]

    output_params_core_attrs = \
        output_params.get_core_attrs(deep_copy=False)

    image_params = \
        output_params_core_attrs["image_params"]
    
    image_params_core_attrs = \
        image_params.get_core_attrs(deep_copy=False)
    
    wavefunction_output_is_to_be_saved = \
        image_params_core_attrs["save_wavefunctions"]

    return wavefunction_output_is_to_be_saved



def _wavefunction_output_filenames(sim_params):
    sim_params_core_attrs = sim_params.get_core_attrs(deep_copy=False)
    
    output_params = sim_params_core_attrs["output_params"]

    output_params_core_attrs = output_params.get_core_attrs(deep_copy=False)

    output_dirname = output_params_core_attrs["output_dirname"]

    num_atomic_config_subsets = _num_atomic_config_subsets(sim_params)

    filenames = tuple()
    for atomic_config_subset_idx in range(num_atomic_config_subsets):
        unformatted_basename = "hrtem_sim_wavefunction_output_of_subset_{}.h5"
        basename = unformatted_basename.format(atomic_config_subset_idx)
        filenames += (output_dirname + "/" + basename,)

    return filenames



def _num_atomic_config_subsets(sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]

    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)

    sample_specification = \
        hrtem_system_model_params_core_attrs["sample_specification"]

    kwargs = \
        {"sample_specification": sample_specification}
    num_atomic_config_subsets = \
        prismatique.sample._num_frozen_phonon_config_subsets(**kwargs)

    return num_atomic_config_subsets



def _check_data_size(sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]

    output_params = \
        sim_params_core_attrs["output_params"]

    output_params_core_attrs = \
        output_params.get_core_attrs(deep_copy=False)

    output_params = sim_params_core_attrs["output_params"]
    max_data_size = output_params_core_attrs["max_data_size"]

    kwargs = {"hrtem_system_model_params": hrtem_system_model_params,
              "output_params": output_params}
    output_data_size = prismatique.hrtem.output._data_size(**kwargs)

    current_func_name = "_check_data_size"

    if max_data_size < output_data_size:
        unformatted_err_msg = globals()[current_func_name+"_err_msg_1"]
        err_msg = unformatted_err_msg.format(output_data_size, max_data_size)
        raise MemoryError(err_msg)

    return None



def run(sim_params):
    r"""Run HRTEM simulation.

    Parameters
    ----------
    sim_params : :class:`prismatique.hrtem.sim.Params`
        The HRTEM simulation parameters. See the documentation for the class
        :class:`prismatique.hrtem.sim.Params` for a discussion on said 
        parameters.

    Returns
    -------

    """
    params = locals()

    global_symbol_table = globals()
    for param_name in params:
        func_name = "_check_and_convert_" + param_name
        func_alias = global_symbol_table[func_name]
        params[param_name] = func_alias(params)

    kwargs = params
    _run(**kwargs)

    return None



def _run(sim_params):
    rng_seeds = _generate_thermal_rng_seeds_from_sim_params(sim_params)    
    num_atomic_config_subsets = len(rng_seeds)
    _remove_temp_files(sim_params, subset_idx=0, first_or_last_call=True)
    _initialize_output_files(sim_params)

    for atomic_config_subset_idx in range(num_atomic_config_subsets):
        kwargs = {"sim_params": sim_params,
                  "atomic_config_subset_idx": atomic_config_subset_idx,
                  "rng_seeds": rng_seeds}
        _run_prismatic_sims_and_postprocess_output_for_subset(**kwargs)

    _remove_temp_files(sim_params, subset_idx=0, first_or_last_call=True)
            
    _serialize_sim_params(sim_params)

    print("\n\n\n")

    return None



def _generate_thermal_rng_seeds_from_sim_params(sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]

    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)
    
    sample_specification = \
        hrtem_system_model_params_core_attrs["sample_specification"]
    
    kwargs = \
        {"sample_specification": sample_specification}
    module_alias = \
        prismatique.sample
    rng_seeds = \
        module_alias._generate_rng_seeds_from_sample_specification(**kwargs)

    return rng_seeds



def _remove_temp_files(sim_params, subset_idx, first_or_last_call):
    output_param_subset = _output_param_subset(sim_params)
    prismatique.sample._remove_temp_files(output_param_subset,
                                          subset_idx,
                                          first_or_last_call)

    intensity_output_is_not_to_be_saved = \
        not _intensity_output_is_to_be_saved(sim_params)

    if first_or_last_call:
        if intensity_output_is_not_to_be_saved:
            filename = _intensity_output_filename(sim_params)
            if pathlib.Path(filename).is_file():
                pathlib.Path(filename).unlink(missing_ok=True)

    return None



def _initialize_output_files(sim_params):
    _write_metadata_to_output_files(sim_params)
    _initialize_data_in_output_files(sim_params)

    return None



def _write_metadata_to_output_files(sim_params):
    filenames = tuple()
    idx_offset = 0
    if _intensity_output_is_to_be_saved(sim_params):
        filenames += (_intensity_output_filename(sim_params),)
        idx_offset += 1
    if _wavefunction_output_is_to_be_saved(sim_params):
        filenames += _wavefunction_output_filenames(sim_params)

    for filename_idx, filename in enumerate(filenames):
        kwargs = {"sim_params": sim_params, "filename": filename}
        _write_r_x_and_r_y_metadata_to_output_file(**kwargs)
        if filename_idx >= idx_offset:
            _write_tilt_metadata_to_output_file(**kwargs)
            _write_defocus_metadata_to_output_file(**kwargs)

    return None



def _write_r_x_and_r_y_metadata_to_output_file(sim_params, filename):
    path_in_file = "/metadata"
    group_id = h5pywrappers.obj.ID(filename, path_in_file)
    h5pywrappers.group.save(None, group_id, write_mode="w")
    group = h5pywrappers.group.load(group_id, read_only=False)

    for_postprocessed_image = \
        (True
         if (filename == _intensity_output_filename(sim_params))
         else False)

    r_x = _r_x(sim_params, for_postprocessed_image)
    dataset = group.create_dataset(name="r_x", data=r_x, dtype="float32")
    dataset.attrs["dim 1"] = "r_x idx"
    dataset.attrs["units"] = "Å"

    r_y = _r_y(sim_params, for_postprocessed_image)
    dataset = group.create_dataset(name="r_y", data=r_y, dtype="float32")
    dataset.attrs["dim 1"] = "r_y idx"
    dataset.attrs["units"] = "Å"

    group.file.close()

    return None



def _r_x(sim_params, for_postprocessed_image):
    kwargs = {"sim_params": sim_params,
              "navigation_dims": tuple(),
              "signal_dtype": "float"}

    image_set_signal = (_blank_postprocessed_image_set_signal(**kwargs)
                        if for_postprocessed_image
                        else _blank_unprocessed_image_set_signal(**kwargs))

    offset = image_set_signal.axes_manager[0].offset
    size = image_set_signal.axes_manager[0].size
    scale = image_set_signal.axes_manager[0].scale

    r_x = offset + scale*np.arange(size)

    return r_x



def _r_y(sim_params, for_postprocessed_image):
    kwargs = {"sim_params": sim_params,
              "navigation_dims": tuple(),
              "signal_dtype": "float"}

    image_set_signal = (_blank_postprocessed_image_set_signal(**kwargs)
                        if for_postprocessed_image
                        else _blank_unprocessed_image_set_signal(**kwargs))

    offset = image_set_signal.axes_manager[1].offset
    size = image_set_signal.axes_manager[1].size
    scale = image_set_signal.axes_manager[1].scale

    r_y = offset + scale*np.arange(size)

    return r_y



def _blank_postprocessed_image_set_signal(sim_params,
                                          navigation_dims,
                                          signal_dtype):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]
    output_params = \
        sim_params_core_attrs["output_params"]

    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)
    
    sample_specification = \
        hrtem_system_model_params_core_attrs["sample_specification"]

    output_params_core_attrs = \
        output_params.get_core_attrs(deep_copy=False)

    image_params = \
        output_params_core_attrs["image_params"]

    image_params_core_attrs = \
        image_params.get_core_attrs(deep_copy=False)

    postprocessing_seq = \
        image_params_core_attrs["postprocessing_seq"]

    kwargs = \
        {"sample_specification": sample_specification,
         "postprocessing_seq": postprocessing_seq,
         "navigation_dims": navigation_dims,
         "signal_is_cbed_pattern_set": False,
         "signal_dtype": signal_dtype}
    blank_postprocessed_image_set_signal = \
        prismatique._signal._blank_postprocessed_2d_signal(**kwargs)

    return blank_postprocessed_image_set_signal



def _blank_unprocessed_image_set_signal(sim_params,
                                        navigation_dims,
                                        signal_dtype):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]
    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)
    
    sample_specification = \
        hrtem_system_model_params_core_attrs["sample_specification"]

    kwargs = \
        {"sample_specification": sample_specification,
         "navigation_dims": navigation_dims,
         "signal_is_cbed_pattern_set": False,
         "signal_dtype": signal_dtype}
    blank_unprocessed_image_set_signal = \
        prismatique._signal._blank_unprocessed_2d_signal(**kwargs)

    return blank_unprocessed_image_set_signal



def _write_tilt_metadata_to_output_file(sim_params, filename):
    path_in_file = "/metadata"
    group_id = h5pywrappers.obj.ID(filename, path_in_file)
    group = h5pywrappers.group.load(group_id, read_only=False)

    tilt_series = _tilt_series(sim_params)
    dataset = group.create_dataset(name="tilts",
                                   data=tilt_series,
                                   dtype="float32")
    dataset.attrs["dim 1"] = "tilt idx"
    dataset.attrs["units"] = "mrad"

    group.file.close()

    return None



def _tilt_series(sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]
    
    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)
    
    sample_specification = \
        hrtem_system_model_params_core_attrs["sample_specification"]
    tilt_params = \
        hrtem_system_model_params_core_attrs["tilt_params"]
    gun_model_params = \
        hrtem_system_model_params_core_attrs["gun_model_params"]

    gun_model_params_core_attrs = \
        gun_model_params.get_core_attrs(deep_copy=False)
    
    mean_beam_energy = \
        gun_model_params_core_attrs["mean_beam_energy"]
    
    tilt_series = prismatique.tilt._series(sample_specification,
                                           mean_beam_energy,
                                           tilt_params)

    return tilt_series



def _write_defocus_metadata_to_output_file(sim_params, filename):
    path_in_file = "/metadata"
    group_id = h5pywrappers.obj.ID(filename, path_in_file)
    group = h5pywrappers.group.load(group_id, read_only=False)

    defocii = _defocii(sim_params)
    dataset = group.create_dataset(name="defocii",
                                   data=defocii,
                                   dtype="float32")
    dataset.attrs["dim 1"] = "defocus idx"
    dataset.attrs["units"] = "Å"

    group.file.close()

    return None



def _defocii(sim_params):
    # Build probe model simply to extract HRTEM beam defocii more easily.
    probe_model_params = _probe_model_params(sim_params)
    
    wavelength = _wavelength(sim_params)

    C_2_0_mag = _get_C_2_0_mag_from_probe_model_params(probe_model_params)
    Delta_f = wavelength * C_2_0_mag / np.pi

    gauss_hermite_points, _ = \
        _gauss_hermite_points_and_weights(probe_model_params)
    sigma_f = \
        _sigma_f(probe_model_params)
    
    defocal_offsets = np.sqrt(2) * gauss_hermite_points * sigma_f

    defocii = Delta_f + defocal_offsets  # In Å.

    return defocii



def _probe_model_params(sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]

    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)
    
    gun_model_params = \
        hrtem_system_model_params_core_attrs["gun_model_params"]
    lens_model_params = \
        hrtem_system_model_params_core_attrs["lens_model_params"]
    defocal_offset_supersampling = \
        hrtem_system_model_params_core_attrs["defocal_offset_supersampling"]

    kwargs = {"lens_model_params": lens_model_params,
              "gun_model_params": gun_model_params,
              "defocal_offset_supersampling": defocal_offset_supersampling}
    probe_model_params = embeam.stem.probe.ModelParams(**kwargs)

    return probe_model_params



def _wavelength(sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]

    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)
    
    gun_model_params = \
        hrtem_system_model_params_core_attrs["gun_model_params"]

    gun_model_params_core_attrs = \
        gun_model_params.get_core_attrs(deep_copy=False)
    
    mean_beam_energy = \
        gun_model_params_core_attrs["mean_beam_energy"]
    
    wavelength = \
        embeam.wavelength(mean_beam_energy)

    return wavelength



def _get_C_2_0_mag_from_probe_model_params(probe_model_params):
    probe_model_params_core_attrs = \
        probe_model_params.get_core_attrs(deep_copy=False)
    lens_model_params = \
        probe_model_params_core_attrs["lens_model_params"]

    lens_model_params_core_attrs = \
        lens_model_params.get_core_attrs(deep_copy=False)
    coherent_aberrations = \
        lens_model_params_core_attrs["coherent_aberrations"]

    C_2_0_mag = 0
    
    for coherent_aberration in coherent_aberrations:
        coherent_aberration_core_attrs = \
            coherent_aberration.get_core_attrs(deep_copy=False)
        
        mn_pair = (coherent_aberration_core_attrs["m"],
                   coherent_aberration_core_attrs["n"])
        if mn_pair == (2, 0):
            C_2_0_mag = coherent_aberration_core_attrs["C_mag"]
            break

    return C_2_0_mag



def _gauss_hermite_points_and_weights(probe_model_params):
    if probe_model_params.is_coherent:
        gauss_hermite_points = \
            np.array([0])
        gauss_hermite_weights = \
            np.sqrt([np.pi])
    else:
        probe_model_params_core_attrs = \
            probe_model_params.get_core_attrs(deep_copy=False)
        defocal_offset_supersampling = \
            probe_model_params_core_attrs["defocal_offset_supersampling"]
        gauss_hermite_points, gauss_hermite_weights = \
            np.polynomial.hermite.hermgauss(defocal_offset_supersampling)

    return gauss_hermite_points, gauss_hermite_weights



def _sigma_f(probe_model_params):
    probe_model_params_core_attrs = \
        probe_model_params.get_core_attrs(deep_copy=False)

    lens_model_params = probe_model_params_core_attrs["lens_model_params"]
    gun_model_params = probe_model_params_core_attrs["gun_model_params"]

    lens_model_params_core_attrs = \
        lens_model_params.get_core_attrs(deep_copy=False)
    gun_model_params_core_attrs = \
        gun_model_params.get_core_attrs(deep_copy=False)

    C_c = (lens_model_params_core_attrs["chromatic_aberration_coef"]
           * (1e-3/1e-10))
    I = lens_model_params_core_attrs["mean_current"]
    sigma_I = lens_model_params_core_attrs["std_dev_current"]

    E = gun_model_params_core_attrs["mean_beam_energy"]
    sigma_E = gun_model_params_core_attrs["intrinsic_energy_spread"]
    sigma_V = gun_model_params_core_attrs["accel_voltage_spread"]

    sigma_f = C_c * np.sqrt((sigma_E/E)**2
                                  + (2*sigma_I/I)**2
                                  + (sigma_V/E)**2)
    
    return sigma_f



def _initialize_data_in_output_files(sim_params):
    filenames = tuple()
    idx_offset = -1
    if _intensity_output_is_to_be_saved(sim_params):
        filenames += (_intensity_output_filename(sim_params),)
        idx_offset += 1
    if _wavefunction_output_is_to_be_saved(sim_params):
        filenames += _wavefunction_output_filenames(sim_params)

    for filename_idx, filename in enumerate(filenames):
        if filename_idx == idx_offset:
            _initialize_intensity_data_in_output_file(sim_params, filename)
        else:
            atomic_config_subset_idx = filename_idx - 1 - idx_offset
            kwargs = {"sim_params": sim_params,
                      "filename": filename,
                      "atomic_config_subset_idx": atomic_config_subset_idx}
            _initialize_wavefunction_data_in_output_file(**kwargs)

    return None



def _initialize_intensity_data_in_output_file(sim_params, filename):
    r_x = _r_x(sim_params, for_postprocessed_image=True)
    r_y = _r_y(sim_params, for_postprocessed_image=True)

    dataset_shape = (len(r_y), len(r_x))

    path_in_file = "/data"
    group_id = h5pywrappers.obj.ID(filename, path_in_file)
    h5pywrappers.group.save(None, group_id, write_mode="a")
    group = h5pywrappers.group.load(group_id, read_only=False)
    
    dataset = group.create_dataset(name="intensity_image",
                                   shape=dataset_shape,
                                   dtype="float32",
                                   fillvalue=0)

    dataset.attrs["dim 1"] = "r_y idx"
    dataset.attrs["dim 2"] = "r_x idx"
    dataset.attrs["units"] = "dimensionless"

    group.file.close()

    return None



def _initialize_wavefunction_data_in_output_file(sim_params,
                                                 filename,
                                                 atomic_config_subset_idx):
    kwargs = {"sim_params": sim_params,
              "atomic_config_subset_idx": atomic_config_subset_idx}
    num_atomic_configs_in_subset = _num_atomic_configs_in_subset(**kwargs)

    defocii = _defocii(sim_params)

    tilt_series = _tilt_series(sim_params)

    r_x = _r_x(sim_params, for_postprocessed_image=False)
    r_y = _r_y(sim_params, for_postprocessed_image=False)

    dataset_shape = (num_atomic_configs_in_subset,
                     len(defocii),
                     len(tilt_series),
                     len(r_y),
                     len(r_x))

    path_in_file = "/data"
    group_id = h5pywrappers.obj.ID(filename, path_in_file)
    h5pywrappers.group.save(None, group_id, write_mode="a")
    group = h5pywrappers.group.load(group_id, read_only=False)
    
    dataset = group.create_dataset(name="image_wavefunctions",
                                   shape=dataset_shape,
                                   dtype="complex64",
                                   fillvalue=0j)

    dataset.attrs["dim 1"] = "atomic config idx"
    dataset.attrs["dim 2"] = "defocus idx"
    dataset.attrs["dim 3"] = "tilt idx"
    dataset.attrs["dim 4"] = "r_y idx"
    dataset.attrs["dim 5"] = "r_x idx"
    dataset.attrs["units"] = "dimensionless"

    group.file.close()

    return None



def _num_atomic_configs_in_subset(sim_params, atomic_config_subset_idx):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]
    
    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)
    
    sample_specification = \
        hrtem_system_model_params_core_attrs["sample_specification"]
    
    kwargs = \
        {"sample_specification": sample_specification,
         "subset_idx": atomic_config_subset_idx}
    num_atomic_configs_in_subset = \
        prismatique.sample._num_frozen_phonon_configs_in_subset(**kwargs)

    return num_atomic_configs_in_subset



def _initialize_pyprismatic_sim_obj(sim_params):
    pyprismatic_sim_obj = pyprismatic.Metadata()

    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]
    output_params = \
        sim_params_core_attrs["output_params"]
    worker_params = \
        sim_params_core_attrs["worker_params"]

    output_params_core_attrs = \
        output_params.get_core_attrs(deep_copy=False)
    output_dirname = \
        output_params_core_attrs["output_dirname"]

    module_alias = prismatique.sample
    func_alias = module_alias._set_pyprismatic_sim_obj_attrs_to_default_values
    func_alias(pyprismatic_sim_obj, output_dirname)
    
    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "hrtem_system_model_params": hrtem_system_model_params,
              "output_dirname": output_dirname}
    _unpack_hrtem_system_model_params_into_pyprismatic_sim_obj(**kwargs)

    unpack_output_params_into_pyprismatic_sim_obj = \
        prismatique.sample._unpack_output_params_into_pyprismatic_sim_obj
    unpack_worker_params_into_pyprismatic_sim_obj = \
        prismatique.sample._unpack_worker_params_into_pyprismatic_sim_obj

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "output_params": output_params,
              "sample_specification": None}
    unpack_output_params_into_pyprismatic_sim_obj(**kwargs)

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "worker_params": worker_params}
    unpack_worker_params_into_pyprismatic_sim_obj(**kwargs)

    return pyprismatic_sim_obj



def _unpack_hrtem_system_model_params_into_pyprismatic_sim_obj(
        pyprismatic_sim_obj, hrtem_system_model_params, output_dirname):
    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)
    
    sample_specification = \
        hrtem_system_model_params_core_attrs["sample_specification"]
    gun_model_params = \
        hrtem_system_model_params_core_attrs["gun_model_params"]
    lens_model_params = \
        hrtem_system_model_params_core_attrs["lens_model_params"]
    tilt_params = \
        hrtem_system_model_params_core_attrs["tilt_params"]

    unpack_sample_specification_into_pyprismatic_sim_obj = \
        prismatique.sample._unpack_sample_specification_into_pyprismatic_sim_obj
    unpack_gun_model_params_into_pyprismatic_sim_obj = \
        prismatique.sample._unpack_gun_model_params_into_pyprismatic_sim_obj
    unpack_lens_model_params_into_pyprismatic_sim_obj = \
        prismatique.sample._unpack_lens_model_params_into_pyprismatic_sim_obj
    unpack_tilt_params_into_pyprismatic_sim_obj = \
        prismatique.sample._unpack_tilt_params_into_pyprismatic_sim_obj

    unpack_sample_specification_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                         output_dirname,
                                                         sample_specification)
    unpack_gun_model_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                     gun_model_params)
    unpack_lens_model_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                      lens_model_params,
                                                      output_dirname)
    unpack_tilt_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                tilt_params)

    return None



def _run_prismatic_sims_and_postprocess_output_for_subset(
        sim_params,
        atomic_config_subset_idx,
        rng_seeds):
    defocii = _defocii(sim_params)  # In Å.
    num_defocii = len(defocii)

    for defocus_idx in range(num_defocii):
        func_alias = \
            _run_prismatic_sims_and_postprocess_output_for_subset_and_defocus

        kwargs = {"sim_params": sim_params,
                  "atomic_config_subset_idx": atomic_config_subset_idx,
                  "defocus_idx": defocus_idx,
                  "rng_seeds": rng_seeds}
        func_alias(**kwargs)

    _remove_temp_files(sim_params,
                       subset_idx=atomic_config_subset_idx,
                       first_or_last_call=False)

    return None



def _run_prismatic_sims_and_postprocess_output_for_subset_and_defocus(
        sim_params,
        atomic_config_subset_idx,
        defocus_idx,
        rng_seeds):
    pyprismatic_sim_obj = _initialize_pyprismatic_sim_obj(sim_params)
    _update_pyprismatic_sim_obj_for_next_prismatic_sim(pyprismatic_sim_obj,
                                                       sim_params,
                                                       atomic_config_subset_idx,
                                                       defocus_idx,
                                                       rng_seeds)

    
    # Run prismatic simulation.
    _call_pyprismatic_sim_obj_go(pyprismatic_sim_obj)
    gc.collect()

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "sim_params": sim_params,
              "atomic_config_subset_idx": atomic_config_subset_idx,
              "defocus_idx": defocus_idx}
    _postprocess_and_reorganize_current_prismatic_sim_output(**kwargs)
    
    return None



def _update_pyprismatic_sim_obj_for_next_prismatic_sim(pyprismatic_sim_obj,
                                                       sim_params,
                                                       atomic_config_subset_idx,
                                                       defocus_idx,
                                                       rng_seeds):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]
    
    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)
    
    sample_specification = \
        hrtem_system_model_params_core_attrs["sample_specification"]

    output_param_subset = _output_param_subset(sim_params)
    defocii = _defocii(sim_params)  # In Å.

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "sample_specification": sample_specification,
              "output_param_subset": output_param_subset,
              "subset_idx": atomic_config_subset_idx,
              "defocus_idx": defocus_idx,
              "defocii": defocii,
              "rng_seeds": rng_seeds}
    module_alias = prismatique.sample
    module_alias._update_pyprismatic_sim_obj_for_next_prismatic_sim(**kwargs)

    return None



def _call_pyprismatic_sim_obj_go(pyprismatic_sim_obj):
    pyprismatic_sim_obj.go()

    return None



def _postprocess_and_reorganize_current_prismatic_sim_output(
        pyprismatic_sim_obj, sim_params, atomic_config_subset_idx, defocus_idx):
    if defocus_idx == 0:
        output_param_subset = _output_param_subset(sim_params)

        module_alias = \
            prismatique.sample
        func_alias = \
            module_alias._move_sample_specification_output_to_separate_file
        
        kwargs = {"output_dirname": output_param_subset["output_dirname"],
                  "new_output_filename": ""}
        if pyprismatic_sim_obj.savePotentialSlices:
            kwargs["sample_specification_type"] = "potential_slice_subset"
            func_alias(**kwargs)
            gc.collect()
    
    kwargs = {"sim_params": sim_params,
              "atomic_config_subset_idx": atomic_config_subset_idx}
    num_atomic_configs_in_subset = _num_atomic_configs_in_subset(**kwargs)

    for atomic_config_idx in range(num_atomic_configs_in_subset):
        triplet_of_indices = \
            {"atomic_config_subset_idx": atomic_config_subset_idx,
             "atomic_config_idx": atomic_config_idx,
             "defocus_idx": defocus_idx}
            
        kwargs = {"sim_params": sim_params,
                  "triplet_of_indices": triplet_of_indices}
        func_alias = _postprocess_and_reorganize_image_subset
        func_alias(**kwargs)
        gc.collect()

    return None



def _postprocess_and_reorganize_image_subset(sim_params, triplet_of_indices):
    kwargs = \
        {"sim_params": sim_params,
         "triplet_of_indices": triplet_of_indices}
    unprocessed_complex_image_subset_signal = \
        _load_unprocessed_complex_image_subset_signal(**kwargs)

    _apply_objective_aperture(unprocessed_complex_image_subset_signal,
                              sim_params)

    if _wavefunction_output_is_to_be_saved(sim_params):
        output_data = unprocessed_complex_image_subset_signal.data
        filenames = _wavefunction_output_filenames(sim_params)
        filename = filenames[triplet_of_indices["atomic_config_subset_idx"]]
        kwargs = {"sim_params": sim_params,
                  "unprocessed_complex_image_subset": output_data,
                  "filename": filename,
                  "triplet_of_indices": triplet_of_indices}
        _write_unprocessed_complex_image_subset_to_output_file(**kwargs)

    if _intensity_output_is_to_be_saved(sim_params):
        kwargs = \
            {"unprocessed_intensity_image_subset_signal": \
             empix.abs_sq(unprocessed_complex_image_subset_signal),
             "sim_params": \
             sim_params}
        postprocessed_intensity_image_subset_signal = \
            _postprocess_intensity_image_subset_signal(**kwargs)

        filename = _intensity_output_filename(sim_params)
        path_in_file = "/data/intensity_image"
        dataset_id = h5pywrappers.obj.ID(filename, path_in_file)
        dataset = h5pywrappers.dataset.load(dataset_id, read_only=False)

        kwargs = {"triplet_of_indices": \
                  triplet_of_indices,
                  "postprocessed_intensity_image_subset_signal": \
                  postprocessed_intensity_image_subset_signal,
                  "sim_params": \
                  sim_params,
                  "dataset": \
                  dataset}
        _update_intensity_data_in_output_dataset(**kwargs)

        dataset.file.close()

    return None



def _load_unprocessed_complex_image_subset_signal(sim_params,
                                                  triplet_of_indices):
    output_dirname = _output_param_subset(sim_params)["output_dirname"]
    filename = output_dirname + "/prismatic_output.h5"

    atomic_config_idx = triplet_of_indices["atomic_config_idx"]
    atomic_config_subset_idx = triplet_of_indices["atomic_config_subset_idx"]

    atomic_config_idx_str = str(atomic_config_idx).rjust(4, "0")

    unformatted_path_in_file = ("4DSTEM_simulation/data/realslices"
                                "/HRTEM_fp{}/data")
    path_in_file = unformatted_path_in_file.format(atomic_config_idx_str)
    dataset_id = h5pywrappers.obj.ID(filename, path_in_file)
    dataset = h5pywrappers.dataset.load(dataset_id, read_only=True)

    kwargs = {"sim_params": sim_params,
              "atomic_config_subset_idx": atomic_config_subset_idx}
    num_atomic_configs_in_subset = _num_atomic_configs_in_subset(**kwargs)
    
    unprocessed_complex_image_subset_data = \
        np.transpose(dataset[()], axes=(2, 1, 0))[:, ::-1, :]

    dataset.file.close()

    kwargs = \
        {"sim_params": sim_params,
         "navigation_dims": (unprocessed_complex_image_subset_data.shape[0],),
         "signal_dtype": "complex"}
    unprocessed_complex_image_subset_signal = \
        _blank_unprocessed_image_set_signal(**kwargs)
    unprocessed_complex_image_subset_signal.data = \
        unprocessed_complex_image_subset_data

    return unprocessed_complex_image_subset_signal



def _apply_objective_aperture(unprocessed_complex_image_subset_signal,
                              sim_params):
    angular_mask = _angular_mask(sim_params)

    signal_data = unprocessed_complex_image_subset_signal.data
    signal_rank = len(signal_data.shape)
    multi_dim_slice = ((slice(None),)*(signal_rank-2)
                       + (slice(None, None, -1), slice(None)))

    temp_data = \
        np.fft.fft2(signal_data[multi_dim_slice]) * angular_mask
    unprocessed_complex_image_subset_signal.data = \
        np.fft.ifft2(temp_data)[multi_dim_slice]

    return None



def _angular_mask(sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]

    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)
    
    objective_aperture_params = \
        hrtem_system_model_params_core_attrs["objective_aperture_params"]

    objective_aperture_params_core_attrs = \
        objective_aperture_params.get_core_attrs(deep_copy=False)
    
    offset = \
        objective_aperture_params_core_attrs["offset"]
    window = \
        objective_aperture_params_core_attrs["window"]
    
    angular_mesh = _angular_mesh(sim_params)

    scattering_angles_relative_to_aperture_offset = \
        prismatique.sample._rel_tilts
    rel_angles = \
        scattering_angles_relative_to_aperture_offset(angular_mesh, offset)

    min_r_angle, max_r_angle = np.array(window) / 1000  # In rads.

    angular_mask = (rel_angles <= max_r_angle) * (rel_angles >= min_r_angle)

    return angular_mask



def _angular_mesh(sim_params):
    # Note that the angular mesh generated here is different from that in the
    # module :mod:`prismatique.sample`.
    
    r_x = _r_x(sim_params, for_postprocessed_image=False)
    r_y = _r_y(sim_params, for_postprocessed_image=False)
    wavelength = _wavelength(sim_params)

    n_x = len(r_x)
    n_y = len(r_y)
    
    Delta_tilde_x = r_x[1] - r_x[0]
    Delta_tilde_y = -(r_y[1] - r_y[0])

    k_x = prismatique.sample._FFT_1D_freqs(n_x, Delta_tilde_x)
    k_y = prismatique.sample._FFT_1D_freqs(n_y, Delta_tilde_y)
    k_mesh = np.meshgrid(k_x, k_y, indexing="xy")

    angular_mesh = (k_mesh[0] * wavelength, k_mesh[1] * wavelength)

    return angular_mesh



def _write_unprocessed_complex_image_subset_to_output_file(
        sim_params,
        unprocessed_complex_image_subset,
        filename,
        triplet_of_indices):
    path_in_file = "/data/image_wavefunctions"
    dataset_id = h5pywrappers.obj.ID(filename, path_in_file)
    
    multi_dim_slice = (triplet_of_indices["atomic_config_idx"],
                       triplet_of_indices["defocus_idx"],
                       slice(None),
                       slice(None),
                       slice(None))
    datasubset_id = h5pywrappers.datasubset.ID(dataset_id, multi_dim_slice)

    datasubset = unprocessed_complex_image_subset
    h5pywrappers.datasubset.save(datasubset, datasubset_id)

    return None



def _postprocess_intensity_image_subset_signal(
        unprocessed_intensity_image_subset_signal, sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    output_params = \
        sim_params_core_attrs["output_params"]

    output_params_core_attrs = \
        output_params.get_core_attrs(deep_copy=False)

    image_params = \
        output_params_core_attrs["image_params"]

    image_params_core_attrs = \
        image_params.get_core_attrs(deep_copy=False)
    
    postprocessing_seq = \
        image_params_core_attrs["postprocessing_seq"]

    kwargs = \
        {"input_signal": unprocessed_intensity_image_subset_signal,
         "postprocessing_seq": postprocessing_seq}
    postprocessed_intensity_image_subset_signal = \
        prismatique._signal._postprocess_2d_signal(**kwargs)

    return postprocessed_intensity_image_subset_signal



def _update_intensity_data_in_output_dataset(
        triplet_of_indices,
        postprocessed_intensity_image_subset_signal,
        sim_params,
        dataset):
    atomic_config_idx = triplet_of_indices["atomic_config_idx"]
    atomic_config_subset_idx = triplet_of_indices["atomic_config_subset_idx"]
    defocus_idx = triplet_of_indices["defocus_idx"]

    signal_data = postprocessed_intensity_image_subset_signal.data
    tilt_weights = _tilt_weights(sim_params)

    dataset[()] += (np.einsum('ijk,i->jk', signal_data, tilt_weights)
                    * _w_f_l(sim_params, l=defocus_idx)
                    / _total_num_frozen_phonon_configs(sim_params)
                    / np.sqrt(np.pi))

    kwargs = {"sim_params": sim_params,
              "atomic_config_subset_idx": atomic_config_subset_idx}
    num_atomic_configs_in_subset = _num_atomic_configs_in_subset(**kwargs)
    
    num_atomic_config_subsets = _num_atomic_config_subsets(sim_params)
    num_defocii = len(_defocii(sim_params))

    sim_params_core_attrs = sim_params.get_core_attrs(deep_copy=False)
    output_params = sim_params_core_attrs["output_params"]
    
    output_params_core_attrs = output_params.get_core_attrs(deep_copy=False)
    image_params = output_params_core_attrs["image_params"]

    if ((atomic_config_subset_idx == num_atomic_config_subsets-1)
        and (atomic_config_idx == num_atomic_configs_in_subset-1)
        and (defocus_idx == num_defocii-1)):
        image_params_core_attrs = \
            image_params.get_core_attrs(deep_copy=False)
        avg_num_electrons_per_postprocessed_image = \
            image_params_core_attrs["avg_num_electrons_per_postprocessed_image"]
        
        datasubset = dataset[()].clip(min=0)
        datasubset *= (avg_num_electrons_per_postprocessed_image
                       / np.sum(datasubset, axis=(-2, -1)))
        if image_params_core_attrs["apply_shot_noise"]:
            datasubset = np.random.poisson(datasubset)
        dataset[()] = datasubset

    return None



def _tilt_weights(sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]

    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)
    
    tilt_params = \
        hrtem_system_model_params_core_attrs["tilt_params"]

    tilt_params_core_attrs = \
        tilt_params.get_core_attrs(deep_copy=False)

    tilt_offset = tilt_params_core_attrs["offset"]  # In mrads.
    tilt_spread = tilt_params_core_attrs["spread"]  # In mrads.
    tilt_series = np.array(_tilt_series(sim_params))  # In mrads.
    rel_tilts = np.linalg.norm(tilt_series - tilt_offset, axis=1)  # In mrads.

    if tilt_spread > 0:
        exp_arg = -0.5*(rel_tilts/tilt_spread)*(rel_tilts/tilt_spread)
        tilt_weights = np.exp(exp_arg)
    else:
        tilt_weights = (rel_tilts == 0.0).astype(float)
    tilt_weights /= np.linalg.norm(tilt_weights)

    return tilt_weights



def _w_f_l(sim_params, l):
    # Build probe model simply to extract HRTEM beam defocus weights more
    # easily.
    probe_model_params = _probe_model_params(sim_params)

    _, gauss_hermite_weights = \
        _gauss_hermite_points_and_weights(probe_model_params)

    w_f_l = gauss_hermite_weights[l]

    return w_f_l



def _total_num_frozen_phonon_configs(sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    hrtem_system_model_params = \
        sim_params_core_attrs["hrtem_system_model_params"]

    hrtem_system_model_params_core_attrs = \
        hrtem_system_model_params.get_core_attrs(deep_copy=False)
    
    sample_specification = \
        hrtem_system_model_params_core_attrs["sample_specification"]

    kwargs = \
        {"sample_specification": sample_specification}
    total_num_frozen_phonon_configs = \
        prismatique.sample._total_num_frozen_phonon_configs(**kwargs)
    
    return total_num_frozen_phonon_configs



def _serialize_sim_params(sim_params):
    sim_params_core_attrs = \
        sim_params.get_core_attrs(deep_copy=False)
    
    output_params = \
        sim_params_core_attrs["output_params"]

    output_params_core_attrs = \
        output_params.get_core_attrs(deep_copy=False)
    
    output_dirname = \
        output_params_core_attrs["output_dirname"]

    filename = \
        output_dirname + "/hrtem_sim_params.json"

    kwargs = {"filename": filename, "overwrite": True}
    sim_params.dump(**kwargs)

    return None
    


###########################
## Define error messages ##
###########################

_check_data_size_err_msg_1 = \
    ("The data size of the output of the HRTEM simulation to be performed is "
     "{} bytes, exceeding the maximum allowed size of {} bytes specified by "
     "the simulation parameters.")
