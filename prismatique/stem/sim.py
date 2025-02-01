r"""For running STEM simulations.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For performing deep copies of objects.
import copy

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

# For conveniently extracting certain properties of probe models.
import embeam

# For postprocessing ``hyperspy`` signals.
import empix

# For creating a :obj:`pyprismatic.Metadata` object that is responsible for
# running the ``prismatic`` simulation.
import pyprismatic



# For validating, pre-serializing, and de-pre-serializing instances of the
# classes :class:`prismatique.worker.Params`,
# :class:`prismatique.stem.system.ModelParams`, and
# :class:`prismatique.stem.output.Params`.
import prismatique.worker
import prismatique.stem.system
import prismatique.stem.output

# For validating instances of the classes
# :class:`prismatique.sample.ModelParams`,
# :class:`prismatique.sample.PotentialSliceSubsetIDs`,
# :class:`prismatique.sample.SMatrixSubsetIDs`, and
# :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`; for calculating
# quantities related to the modelling of the sample; for validating certain
# filenames; and for importing various other helper functions.
import prismatique.sample

# For generating probe positions.
import prismatique.scan

# For generating and postprocessing CBED patterns.
import prismatique._signal



############################
## Authorship information ##
############################

__author__     = "Matthew Fitzpatrick"
__copyright__  = "Copyright 2023"
__credits__    = ["Matthew Fitzpatrick"]
__maintainer__ = "Matthew Fitzpatrick"
__email__      = "mrfitzpa@uvic.ca"
__status__     = "Development"



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params",
           "run"]



def _check_and_convert_stem_system_model_params(ctor_params):
    func_alias = _check_and_convert_stem_system_model_params_and_output_params
    stem_system_model_params, _ = func_alias(ctor_params)
    
    return stem_system_model_params



def _check_and_convert_stem_system_model_params_and_output_params(ctor_params):
    check_and_convert_stem_system_model_params = \
        prismatique.stem.system._check_and_convert_stem_system_model_params
    check_and_convert_output_params = \
        prismatique.stem.output._check_and_convert_output_params
    check_sample_specification_wrt_output_params = \
        prismatique.stem.output._check_sample_specification_wrt_output_params

    stem_system_model_params = \
        check_and_convert_stem_system_model_params(ctor_params)
    
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]
    output_params = \
        check_and_convert_output_params(ctor_params)
    format_arg = \
        "stem_system_model_params.core_attrs['sample_specification']"

    check_sample_specification_wrt_output_params(sample_specification,
                                                 output_params,
                                                 format_arg)
    
    return stem_system_model_params, output_params



def _pre_serialize_stem_system_model_params(stem_system_model_params):
    pre_serialize_stem_system_model_params = \
        prismatique.stem.system._pre_serialize_stem_system_model_params
    serializable_rep = \
        pre_serialize_stem_system_model_params(stem_system_model_params)

    return serializable_rep



def _de_pre_serialize_stem_system_model_params(serializable_rep):
    de_pre_serialize_stem_system_model_params = \
        prismatique.stem.system._de_pre_serialize_stem_system_model_params
    stem_system_model_params = \
        de_pre_serialize_stem_system_model_params(serializable_rep)

    return stem_system_model_params



def _check_and_convert_output_params(ctor_params):
    func_alias = _check_and_convert_stem_system_model_params_and_output_params
    _, output_params = func_alias(ctor_params)
    
    return output_params



def _pre_serialize_output_params(output_params):
    pre_serialize_output_params = \
        prismatique.stem.output._pre_serialize_output_params
    serializable_rep = \
        pre_serialize_output_params(output_params)

    return serializable_rep



def _de_pre_serialize_output_params(serializable_rep):
    de_pre_serialize_output_params = \
        prismatique.stem.output._de_pre_serialize_output_params
    output_params = \
        de_pre_serialize_output_params(serializable_rep)

    return output_params



def _check_and_convert_worker_params(ctor_params):
    check_and_convert_worker_params = \
        prismatique.worker._check_and_convert_worker_params
    worker_params = \
        check_and_convert_worker_params(ctor_params)
    
    return worker_params



def _pre_serialize_worker_params(worker_params):
    pre_serialize_worker_params = \
        prismatique.worker._pre_serialize_worker_params
    serializable_rep = \
        pre_serialize_worker_params(worker_params)

    return serializable_rep



def _de_pre_serialize_worker_params(serializable_rep):
    de_pre_serialize_worker_params = \
        prismatique.worker._de_pre_serialize_worker_params
    worker_params = \
        de_pre_serialize_worker_params(serializable_rep)

    return worker_params



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The STEM simulation parameters.

    Parameters
    ----------
    stem_system_model_params : :class:`prismatique.stem.system.ModelParams`
        The simulation parameters related to the modelling of the STEM
        system. See the documentation for the class
        :class:`prismatique.stem.system.ModelParams` for a discussion on said
        parameters.
    output_params : :class:`prismatique.stem.output.Params` | `None`, optional
        The output parameters for the STEM simulation. See the documentation for
        the class :class:`prismatique.stem.output.Params` for a discussion on
        said parameters. If ``output_params`` is set to `None` [i.e. the default
        value], then the aforementioned simulation parameters are set to default
        values.
    worker_params : :class:`prismatique.worker.Params` | `None`, optional
        The simulation parameters related to GPU and CPU workers. See the
        documentation for the class :class:`prismatique.worker.Params` for a
        discussion on said parameters. If ``worker_params`` is set to `None`
        [i.e. the default value], then the aforementioned simulation parameters 
        are set to default values.

    Attributes
    ----------
    core_attrs : `dict`, read-only
        A `dict` representation of the core attributes: each `dict` key is a
        `str` representing the name of a core attribute, and the corresponding
        `dict` value is the object to which said core attribute is set. The core
        attributes are the same as the construction parameters, except that 
        their values might have been updated since construction.

    """
    _validation_and_conversion_funcs = \
        {"stem_system_model_params": \
         _check_and_convert_stem_system_model_params,
         "output_params": _check_and_convert_output_params,
         "worker_params": _check_and_convert_worker_params}

    _pre_serialization_funcs = \
        {"stem_system_model_params": _pre_serialize_stem_system_model_params,
         "output_params": _pre_serialize_output_params,
         "worker_params": _pre_serialize_worker_params}

    _de_pre_serialization_funcs = \
        {"stem_system_model_params": _de_pre_serialize_stem_system_model_params,
         "output_params": _de_pre_serialize_output_params,
         "worker_params": _de_pre_serialize_worker_params}
    
    def __init__(self,
                 stem_system_model_params,
                 output_params=None,
                 worker_params=None):
        ctor_params = {"stem_system_model_params": stem_system_model_params,
                       "output_params": output_params,
                       "worker_params": worker_params}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_sim_params(ctor_params):
    sim_params = copy.deepcopy(ctor_params["sim_params"])    
    kwargs = {"obj": sim_params,
              "obj_name": "sim_params",
              "accepted_types": (Params,)}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return sim_params



def _pre_serialize_sim_params(sim_params):
    serializable_rep = sim_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_sim_params(serializable_rep):
    sim_params = Params.de_pre_serialize(serializable_rep)

    return sim_params



def run(sim_params):
    r"""Run STEM simulation.

    Parameters
    ----------
    sim_params : :class:`prismatique.stem.sim.Params`
        The STEM simulation parameters. See the documentation for the class
        :class:`prismatique.stem.sim.Params` for a discussion on said 
        parameters.

    Returns
    -------

    """
    _check_sim_params(sim_params)
    _run(sim_params)
    _serialize_sim_params(sim_params)
    print("\n\n\n")

    return None



def _check_sim_params(sim_params):
    temp_ctor_params = {"sim_params": sim_params}
    sim_params = _check_and_convert_sim_params(temp_ctor_params)
    
    mod_alias = prismatique.stem.output
    stem_system_model_params = sim_params.core_attrs["stem_system_model_params"]
    mod_alias._check_stem_system_model_params(stem_system_model_params)
    
    _pre_save(sim_params)
    _check_data_size(sim_params)

    return None



def _pre_save(sim_params):
    stem_system_model_params = sim_params.core_attrs["stem_system_model_params"]
    
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]
    
    output_params = sim_params.core_attrs["output_params"]
    base_output_params = output_params.core_attrs["base_params"]
    alg_specific_output_params = output_params.core_attrs["alg_specific_params"]
    
    output_dirname = base_output_params.core_attrs["output_dirname"]

    if base_output_params.core_attrs["save_potential_slices"]:
        unformatted_basename = "potential_slices_of_subset_{}.h5"
        prismatique.sample._pre_save(sample_specification,
                                     output_dirname,
                                     unformatted_basename)
    if "save_S_matrices" in alg_specific_output_params.core_attrs:
        if alg_specific_output_params.core_attrs["save_S_matrices"]:
            unformatted_basename = "S_matrices_of_subset_{}.h5"
            prismatique.sample._pre_save(sample_specification,
                                         output_dirname,
                                         unformatted_basename)

    filenames = []
    filenames.append(_intensity_output_filename(sim_params))
    if _wavefunction_output_is_to_be_saved(sim_params):
        filenames += _wavefunction_output_filenames(sim_params)

    for filename in filenames:
        if pathlib.Path(filename).is_file():
            pathlib.Path(filename).unlink(missing_ok=True)

    prismatique.sample._check_hdf5_filenames(filenames)

    return None



def _intensity_output_is_to_be_saved(sim_params):
    stem_system_model_params = sim_params.core_attrs["stem_system_model_params"]
    output_params = sim_params.core_attrs["output_params"]
    base_output_params = output_params.core_attrs["base_params"]

    kwargs = \
        {"base_output_params": base_output_params}
    intensity_output_is_to_be_saved = \
        prismatique.stem.output._intensity_output_is_to_be_saved(**kwargs)

    return intensity_output_is_to_be_saved



def _intensity_output_filename(sim_params):
    output_dirname = _output_param_subset(sim_params)["output_dirname"]
    filename = output_dirname + "/stem_sim_intensity_output.h5"

    return filename



def _wavefunction_output_is_to_be_saved(sim_params):
    stem_system_model_params = sim_params.core_attrs["stem_system_model_params"]
    output_params = sim_params.core_attrs["output_params"]
    base_output_params = output_params.core_attrs["base_params"]

    kwargs = \
        {"base_output_params": base_output_params}
    wavefunction_output_is_to_be_saved = \
        prismatique.stem.output._wavefunction_output_is_to_be_saved(**kwargs)

    return wavefunction_output_is_to_be_saved



def _wavefunction_output_filenames(sim_params):
    output_params = sim_params.core_attrs["output_params"]
    base_output_params = output_params.core_attrs["base_params"]    
    output_dirname = base_output_params.core_attrs["output_dirname"]

    num_atomic_config_subsets = _num_atomic_config_subsets(sim_params)

    filenames = []
    for atomic_config_subset_idx in range(num_atomic_config_subsets):
        unformatted_basename = "stem_sim_wavefunction_output_of_subset_{}.h5"
        basename = unformatted_basename.format(atomic_config_subset_idx)
        filenames.append(output_dirname + "/" + basename)

    return filenames



def _num_atomic_config_subsets(sim_params):
    stem_system_model_params = \
        sim_params.core_attrs["stem_system_model_params"]
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]

    kwargs = \
        {"sample_specification": sample_specification}
    num_atomic_config_subsets = \
        prismatique.sample._num_frozen_phonon_config_subsets(**kwargs)

    return num_atomic_config_subsets



def _check_data_size(sim_params):
    stem_system_model_params = sim_params.core_attrs["stem_system_model_params"]

    output_params = sim_params.core_attrs["output_params"]
    base_output_params = output_params.core_attrs["base_params"]
    max_data_size = base_output_params.core_attrs["max_data_size"]

    kwargs = {"stem_system_model_params": stem_system_model_params,
              "output_params": output_params}
    output_data_size = prismatique.stem.output._data_size(**kwargs)

    if max_data_size < output_data_size:
        unformatted_err_msg = _check_data_size_err_msg_1
        err_msg = unformatted_err_msg.format(output_data_size, max_data_size)
        raise MemoryError(err_msg)

    return None



def _run(sim_params):
    rng_seeds = _generate_thermal_rng_seeds_from_sim_params(sim_params)    
    num_atomic_config_subsets = len(rng_seeds)
    _remove_temp_files(sim_params, subset_idx=0, first_or_last_call=True)
    _initialize_output_files(sim_params)

    try:
        for atomic_config_subset_idx in range(num_atomic_config_subsets):
            kwargs = {"sim_params": sim_params,
                      "atomic_config_subset_idx": atomic_config_subset_idx,
                      "rng_seeds": rng_seeds}
            _run_prismatic_sims_and_postprocess_output_for_subset(**kwargs)

        _remove_temp_files(sim_params, subset_idx=0, first_or_last_call=True)
            
    except BaseException as err:
        _remove_temp_files(sim_params, subset_idx=0, first_or_last_call=True)
        raise err

    return None



def _generate_thermal_rng_seeds_from_sim_params(sim_params):
    stem_system_model_params = \
        sim_params.core_attrs["stem_system_model_params"]
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]
    kwargs = \
        {"sample_specification": sample_specification}
    mod_alias = \
        prismatique.sample
    rng_seeds = \
        mod_alias._generate_rng_seeds_from_sample_specification(**kwargs)

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



def _output_param_subset(sim_params):
    output_params = sim_params.core_attrs["output_params"]
    base_output_params = output_params.core_attrs["base_params"]
    alg_specific_output_params = output_params.core_attrs["alg_specific_params"]
    
    output_dirname = \
        base_output_params.core_attrs["output_dirname"]
    save_potential_slices = \
        base_output_params.core_attrs["save_potential_slices"]

    if "save_S_matrices" in alg_specific_output_params.core_attrs:
        save_S_matrices = \
            alg_specific_output_params.core_attrs["save_S_matrices"]
    else:
        save_S_matrices = False

    output_param_subset = {"output_dirname": output_dirname,
                           "save_potential_slices": save_potential_slices,
                           "save_S_matrices": save_S_matrices}

    return output_param_subset



def _initialize_output_files(sim_params):
    _write_metadata_to_output_files(sim_params)
    _initialize_data_in_output_files(sim_params)

    return None



def _write_metadata_to_output_files(sim_params):
    filenames = []
    filenames.append(_intensity_output_filename(sim_params))
    if _wavefunction_output_is_to_be_saved(sim_params):
        filenames += _wavefunction_output_filenames(sim_params)

    for filename_idx, filename in enumerate(filenames):
        kwargs = {"sim_params": sim_params, "filename": filename}
        _write_probe_position_metadata_to_output_file(**kwargs)
        _write_output_layer_depth_metadata_to_output_file(**kwargs)
        _write_k_x_and_k_y_metadata_to_output_file(**kwargs)
        if filename_idx == 0:
            if _intensity_3d_stem_output_is_to_be_saved(sim_params):
                _write_k_xy_metadata_to_output_file(**kwargs)
        else:
            _write_defocus_metadata_to_output_file(**kwargs)

    return None



def _write_probe_position_metadata_to_output_file(sim_params, filename):
    path_in_file = "/metadata"
    group_id = h5pywrappers.obj.ID(filename, path_in_file)
    h5pywrappers.group.save(None, group_id, write_mode="a")
    group = h5pywrappers.group.load(group_id, read_only=False)

    probe_positions = _generate_probe_positions(sim_params, save=False)
    dataset = group.create_dataset(name="probe_positions",
                                   data=probe_positions,
                                   dtype="float32")
    dataset.attrs["dim 1"] = "probe idx"
    dataset.attrs["dim 2"] = "vector component idx [0->x, 1->y]"
    dataset.attrs["units"] = "Å"
    dataset.attrs["pattern type"] = _scan_pattern_type(sim_params)
    if len(_grid_dims_in_units_of_probe_shifts(sim_params)) != 3:
        dataset.attrs["grid dims in units of probe shifts"] = \
            _grid_dims_in_units_of_probe_shifts(sim_params)

    group.file.close()

    return None



def _generate_probe_positions(sim_params, save):
    output_params = \
        sim_params.core_attrs["output_params"]
    stem_system_model_params = \
        sim_params.core_attrs["stem_system_model_params"]
    
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]
    scan_config = \
        stem_system_model_params.core_attrs["scan_config"]

    base_output_params = output_params.core_attrs["base_params"]
    output_dirname = base_output_params.core_attrs["output_dirname"]

    temp_scan_config_filename = \
        prismatique.sample._generate_temp_scan_config_filename(output_dirname)
    filename = \
        temp_scan_config_filename if save else None

    kwargs = {"sample_specification": sample_specification,
              "scan_config": scan_config,
              "filename": filename}
    probe_positions = prismatique.scan._generate_probe_positions(**kwargs)

    return probe_positions



def _scan_pattern_type(sim_params):
    stem_system_model_params = sim_params.core_attrs["stem_system_model_params"]
    scan_config = stem_system_model_params.core_attrs["scan_config"]
    scan_pattern_type = prismatique.scan._pattern_type(scan_config)

    return scan_pattern_type



def _grid_dims_in_units_of_probe_shifts(sim_params):
    stem_system_model_params = \
        sim_params.core_attrs["stem_system_model_params"]
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]
    scan_config = \
        stem_system_model_params.core_attrs["scan_config"]

    kwargs = {"sample_specification": sample_specification,
              "scan_config": scan_config}
    grid_dims = prismatique.scan._grid_dims_in_units_of_probe_shifts(**kwargs)

    return grid_dims



def _write_output_layer_depth_metadata_to_output_file(sim_params, filename):
    path_in_file = "/metadata"
    group_id = h5pywrappers.obj.ID(filename, path_in_file)
    group = h5pywrappers.group.load(group_id, read_only=False)

    output_layer_depths = _output_layer_depths(sim_params)
    dataset = group.create_dataset(name="output_layer_depths",
                                   data=output_layer_depths,
                                   dtype="float32")
    dataset.attrs["dim 1"] = "output layer idx"
    dataset.attrs["units"] = "Å"

    group.file.close()

    return None



def _output_layer_depths(sim_params):
    stem_system_model_params = sim_params.core_attrs["stem_system_model_params"]
    output_params = sim_params.core_attrs["output_params"]

    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]

    alg_specific_params = output_params.core_attrs["alg_specific_params"]

    kwargs = {"sample_specification": sample_specification,
              "alg_specific_params": alg_specific_params}
    output_layer_depths = prismatique.stem.output._layer_depths(**kwargs)

    return output_layer_depths



def _write_k_x_and_k_y_metadata_to_output_file(sim_params, filename):
    path_in_file = "/metadata"
    group_id = h5pywrappers.obj.ID(filename, path_in_file)
    group = h5pywrappers.group.load(group_id, read_only=False)

    if filename == _intensity_output_filename(sim_params):
        for_postprocessed_dp = True
    else:
        for_postprocessed_dp = False

    k_x = _k_x(sim_params, for_postprocessed_dp)
    dataset = group.create_dataset(name="k_x", data=k_x, dtype="float32")
    dataset.attrs["dim 1"] = "k_x idx"
    dataset.attrs["units"] = "1/Å"

    k_y = _k_y(sim_params, for_postprocessed_dp)
    dataset = group.create_dataset(name="k_y", data=k_y, dtype="float32")
    dataset.attrs["dim 1"] = "k_y idx"
    dataset.attrs["units"] = "1/Å"

    group.file.close()

    return None



def _k_x(sim_params, for_postprocessed_dp):
    kwargs = {"sim_params": sim_params,
              "navigation_dims": tuple(),
              "signal_dtype": "float"}
    if for_postprocessed_dp:
        dp_set_signal = _blank_postprocessed_dp_set_signal(**kwargs)
    else:
        dp_set_signal = _blank_unprocessed_dp_set_signal(**kwargs)

    offset = dp_set_signal.axes_manager[0].offset
    size = dp_set_signal.axes_manager[0].size
    scale = dp_set_signal.axes_manager[0].scale

    k_x = offset + scale*np.arange(size)

    return k_x



def _k_y(sim_params, for_postprocessed_dp):
    kwargs = {"sim_params": sim_params,
              "navigation_dims": tuple(),
              "signal_dtype": "float"}
    if for_postprocessed_dp:
        dp_set_signal = _blank_postprocessed_dp_set_signal(**kwargs)
    else:
        dp_set_signal = _blank_unprocessed_dp_set_signal(**kwargs)

    offset = dp_set_signal.axes_manager[1].offset
    size = dp_set_signal.axes_manager[1].size
    scale = dp_set_signal.axes_manager[1].scale

    k_y = offset + scale*np.arange(size)

    return k_y



def _blank_postprocessed_dp_set_signal(sim_params,
                                       navigation_dims,
                                       signal_dtype):
    stem_system_model_params = sim_params.core_attrs["stem_system_model_params"]
    output_params = sim_params.core_attrs["output_params"]

    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]

    base_output_params = output_params.core_attrs["base_params"]
    cbed_params = base_output_params.core_attrs["cbed_params"]
    postprocessing_seq = cbed_params.core_attrs["postprocessing_seq"]

    kwargs = \
        {"sample_specification": sample_specification,
         "postprocessing_seq": postprocessing_seq,
         "navigation_dims": navigation_dims,
         "signal_is_cbed_pattern_set": True,
         "signal_dtype": signal_dtype}
    blank_postprocessed_dp_set_signal = \
        prismatique._signal._blank_postprocessed_2d_signal(**kwargs)

    return blank_postprocessed_dp_set_signal



def _blank_unprocessed_dp_set_signal(sim_params, navigation_dims, signal_dtype):
    stem_system_model_params = \
        sim_params.core_attrs["stem_system_model_params"]
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]

    kwargs = \
        {"sample_specification": sample_specification,
         "navigation_dims": navigation_dims,
         "signal_is_cbed_pattern_set": True,
         "signal_dtype": signal_dtype}
    blank_unprocessed_dp_set_signal = \
        prismatique._signal._blank_unprocessed_2d_signal(**kwargs)

    return blank_unprocessed_dp_set_signal



def _intensity_2d_stem_output_is_to_be_saved(sim_params):
    output_params = \
        sim_params.core_attrs["output_params"]
    base_output_params = \
        output_params.core_attrs["base_params"]
    radial_step_size_for_3d_stem = \
        base_output_params.core_attrs["radial_step_size_for_3d_stem"]

    if radial_step_size_for_3d_stem > 0:
        intensity_2d_stem_output_is_to_be_saved = True
    else:
        intensity_2d_stem_output_is_to_be_saved = False

    return intensity_2d_stem_output_is_to_be_saved



def _write_k_xy_metadata_to_output_file(sim_params, filename):
    path_in_file = "/metadata"
    group_id = h5pywrappers.obj.ID(filename, path_in_file)
    group = h5pywrappers.group.load(group_id, read_only=False)

    k_xy = _k_xy(sim_params)
    dataset = group.create_dataset(name="k_xy", data=k_xy, dtype="float32")
    dataset.attrs["dim 1"] = "k_xy idx"
    dataset.attrs["units"] = "1/Å"

    group.file.close()

    return None



def _k_xy(sim_params):
    kwargs = {"sim_params": sim_params,
              "navigation_dims": tuple(),
              "signal_dtype": "float"}
    blank_postprocessed_dp_signal = _blank_postprocessed_dp_set_signal(**kwargs)
    
    kwargs = {"input_signal": blank_postprocessed_dp_signal,
              "sim_params": sim_params}
    integrated_signal = _integrate_dp_to_3d_stem_signal(**kwargs)

    offset = integrated_signal.axes_manager[0].offset
    size = integrated_signal.axes_manager[0].size
    scale = integrated_signal.axes_manager[0].scale

    k_xy = offset + scale*np.arange(size)

    return k_xy



def _integrate_dp_to_3d_stem_signal(input_signal, sim_params):
    pixel_area = np.abs(input_signal.axes_manager[-2].scale
                        * input_signal.axes_manager[-1].scale)
    input_signal.data /= pixel_area

    num_bins = _num_bins_in_3d_stem_integration(sim_params, input_signal)
    
    radial_k_range = _radial_k_range_in_3d_stem_integration(sim_params,
                                                            input_signal)

    kwargs = {"center": (0, 0),
              "radial_range": radial_k_range,
              "num_bins": num_bins}
    optional_params = empix.OptionalAzimuthalIntegrationParams(**kwargs)
    integrated_dp_signal = empix.azimuthally_integrate(input_signal,
                                                       optional_params)

    return integrated_dp_signal



def _num_bins_in_3d_stem_integration(sim_params, input_signal):
    stem_system_model_params = sim_params.core_attrs["stem_system_model_params"]
    output_params = sim_params.core_attrs["output_params"]

    kwargs = \
        {"stem_system_model_params": stem_system_model_params,
         "output_params": output_params,
         "input_signal": input_signal}
    num_bins = \
        prismatique.stem.output._num_bins_in_3d_stem_integration(**kwargs)

    return num_bins



def _radial_k_range_in_3d_stem_integration(sim_params, input_signal):
    wavelength = _wavelength(sim_params)
    
    output_params = \
        sim_params.core_attrs["output_params"]
    base_output_params = \
        output_params.core_attrs["base_params"]
    radial_angular_step_size = \
        base_output_params.core_attrs["radial_step_size_for_3d_stem"]
    radial_k_step_size = \
        (radial_angular_step_size / 1000) / wavelength

    num_bins = _num_bins_in_3d_stem_integration(sim_params, input_signal)

    radial_k_range = (0, num_bins*radial_k_step_size)

    return radial_k_range



def _wavelength(sim_params):
    stem_system_model_params = \
        sim_params.core_attrs["stem_system_model_params"]
    probe_model_params = \
        stem_system_model_params.core_attrs["probe_model_params"]
    gun_model_params = \
        probe_model_params.core_attrs["gun_model_params"]
    mean_beam_energy = \
        gun_model_params.core_attrs["mean_beam_energy"]
    wavelength = \
        embeam.wavelength(mean_beam_energy)

    return wavelength



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
    stem_system_model_params = \
        sim_params.core_attrs["stem_system_model_params"]
    probe_model_params = \
        stem_system_model_params.core_attrs["probe_model_params"]
    
    wavelength = _wavelength(sim_params)

    get_C_2_0_mag_from_probe_model_params = \
        embeam.stem.probe._get_C_2_0_mag_from_probe_model_params
    
    C_2_0_mag = get_C_2_0_mag_from_probe_model_params(probe_model_params)
    Delta_f = wavelength * C_2_0_mag / np.pi
    
    sigma_f = probe_model_params._sigma_f
    defocal_offsets = (np.sqrt(2)
                       * probe_model_params._gauss_hermite_points
                       * sigma_f)

    defocii = Delta_f + defocal_offsets  # In Å.

    return defocii



def _initialize_data_in_output_files(sim_params):
    filenames = []
    filenames.append(_intensity_output_filename(sim_params))
    if _wavefunction_output_is_to_be_saved(sim_params):
        filenames += _wavefunction_output_filenames(sim_params)

    for filename_idx, filename in enumerate(filenames):
        if filename_idx == 0:
            _initialize_intensity_data_in_output_file(sim_params, filename)
        else:
            atomic_config_subset_idx = filename_idx - 1
            kwargs = {"sim_params": sim_params,
                      "filename": filename,
                      "atomic_config_subset_idx": atomic_config_subset_idx}
            _initialize_wavefunction_data_in_output_file(**kwargs)

    return None



def _initialize_intensity_data_in_output_file(sim_params, filename):
    _initialize_intensity_4d_stem_data_in_output_file(sim_params, filename)

    if _com_output_is_to_be_saved(sim_params):
        _initialize_com_data_in_output_file(sim_params, filename)
    if _intensity_3d_stem_output_is_to_be_saved(sim_params):
        _initialize_intensity_3d_stem_data_in_output_file(sim_params, filename)
    if _intensity_2d_stem_output_is_to_be_saved(sim_params):
        _initialize_intensity_2d_stem_data_in_output_file(sim_params, filename)

    return None



def _initialize_intensity_4d_stem_data_in_output_file(sim_params, filename):
    output_layer_depths = _output_layer_depths(sim_params)
    probe_positions = _generate_probe_positions(sim_params, save=False)
    k_x = _k_x(sim_params, for_postprocessed_dp=True)
    k_y = _k_y(sim_params, for_postprocessed_dp=True)

    dataset_shape = (len(output_layer_depths),
                     len(probe_positions),
                     len(k_y),
                     len(k_x))

    path_in_file = "/data/4D_STEM"
    group_id = h5pywrappers.obj.ID(filename, path_in_file)
    h5pywrappers.group.save(None, group_id, write_mode="a")
    group = h5pywrappers.group.load(group_id, read_only=False)
    
    dataset = group.create_dataset(name="intensity_DPs",
                                   shape=dataset_shape,
                                   dtype="float32",
                                   fillvalue=0)

    dataset.attrs["dim 1"] = "output layer idx"
    dataset.attrs["dim 2"] = "probe idx"
    dataset.attrs["dim 3"] = "k_y idx"
    dataset.attrs["dim 4"] = "k_x idx"
    dataset.attrs["units"] = "dimensionless"

    group.file.close()

    return None



def _com_output_is_to_be_saved(sim_params):
    output_params = \
        sim_params.core_attrs["output_params"]
    base_output_params = \
        output_params.core_attrs["base_params"]
    com_output_is_to_be_saved = \
        base_output_params.core_attrs["save_com"]

    return com_output_is_to_be_saved



def _initialize_com_data_in_output_file(sim_params, filename):
    output_layer_depths = _output_layer_depths(sim_params)
    probe_positions = _generate_probe_positions(sim_params, save=False)

    dataset_shape = (len(output_layer_depths), 2, len(probe_positions))

    path_in_file = "/data"
    group_id = h5pywrappers.obj.ID(filename, path_in_file)
    group = h5pywrappers.group.load(group_id, read_only=False)
    
    dataset = group.create_dataset(name="center_of_mass_momentum",
                                   shape=dataset_shape,
                                   dtype="float32",
                                   fillvalue=0)

    dataset.attrs["dim 1"] = "output layer idx"
    dataset.attrs["dim 2"] = "vector component idx [0->x, 1->y]"
    dataset.attrs["dim 3"] = "probe idx"
    dataset.attrs["units"] = "1/Å"

    group.file.close()
    
    return None



def _initialize_intensity_3d_stem_data_in_output_file(sim_params, filename):
    output_layer_depths = _output_layer_depths(sim_params)
    probe_positions = _generate_probe_positions(sim_params, save=False)

    k_xy = _k_xy(sim_params)

    dataset_shape = (len(output_layer_depths), len(probe_positions), len(k_xy))

    path_in_file = "/data/3D_STEM"
    group_id = h5pywrappers.obj.ID(filename, path_in_file)
    h5pywrappers.group.save(None, group_id, write_mode="a")
    group = h5pywrappers.group.load(group_id, read_only=False)
    
    dataset = group.create_dataset(name="integrated_intensities",
                                   shape=dataset_shape,
                                   dtype="float32",
                                   fillvalue=0)

    dataset.attrs["dim 1"] = "output layer idx"
    dataset.attrs["dim 2"] = "probe idx"
    dataset.attrs["dim 3"] = "k_xy idx"
    dataset.attrs["units"] = "dimensionless"

    group.file.close()
    
    return None



def _intensity_2d_stem_output_is_to_be_saved(sim_params):
    output_params = \
        sim_params.core_attrs["output_params"]
    base_output_params = \
        output_params.core_attrs["base_params"]
    radial_range_for_2d_stem = \
        base_output_params.core_attrs["radial_range_for_2d_stem"]

    if radial_range_for_2d_stem[0] != radial_range_for_2d_stem[1]:
        intensity_2d_stem_output_is_to_be_saved = True
    else:
        intensity_2d_stem_output_is_to_be_saved = False

    return intensity_2d_stem_output_is_to_be_saved



def _initialize_intensity_2d_stem_data_in_output_file(sim_params, filename):
    output_layer_depths = _output_layer_depths(sim_params)
    probe_positions = _generate_probe_positions(sim_params, save=False)

    dataset_shape = (len(output_layer_depths), len(probe_positions))

    path_in_file = "/data/2D_STEM"
    group_id = h5pywrappers.obj.ID(filename, path_in_file)
    h5pywrappers.group.save(None, group_id, write_mode="a")
    group = h5pywrappers.group.load(group_id, read_only=False)
    
    dataset = group.create_dataset(name="integrated_intensities",
                                   shape=dataset_shape,
                                   dtype="float32",
                                   fillvalue=0)

    dataset.attrs["dim 1"] = "output layer idx"
    dataset.attrs["dim 2"] = "probe idx"
    dataset.attrs["units"] = "dimensionless"

    output_params = \
        sim_params.core_attrs["output_params"]
    base_output_params = \
        output_params.core_attrs["base_params"]
    radial_range_for_2d_stem = \
        base_output_params.core_attrs["radial_range_for_2d_stem"]
    dataset.attrs["lower integration limit in mrads"] = \
        radial_range_for_2d_stem[0]
    dataset.attrs["upper integration limit in mrads"] = \
        radial_range_for_2d_stem[1]

    group.file.close()
    
    return None



def _initialize_wavefunction_data_in_output_file(sim_params,
                                                 filename,
                                                 atomic_config_subset_idx):
    output_layer_depths = _output_layer_depths(sim_params)

    kwargs = {"sim_params": sim_params,
              "atomic_config_subset_idx": atomic_config_subset_idx}
    num_atomic_configs_in_subset = _num_atomic_configs_in_subset(**kwargs)

    defocii = _defocii(sim_params)

    probe_positions = _generate_probe_positions(sim_params, save=False)

    k_x = _k_x(sim_params, for_postprocessed_dp=False)
    k_y = _k_y(sim_params, for_postprocessed_dp=False)

    dataset_shape = (len(output_layer_depths),
                     num_atomic_configs_in_subset,
                     len(defocii),
                     len(probe_positions),
                     len(k_y),
                     len(k_x))

    path_in_file = "/data/4D_STEM"
    group_id = h5pywrappers.obj.ID(filename, path_in_file)
    h5pywrappers.group.save(None, group_id, write_mode="a")
    group = h5pywrappers.group.load(group_id, read_only=False)
    
    dataset = group.create_dataset(name="complex_valued_DPs",
                                   shape=dataset_shape,
                                   dtype="complex64",
                                   fillvalue=0j)

    dataset.attrs["dim 1"] = "output layer idx"
    dataset.attrs["dim 2"] = "atomic config idx"
    dataset.attrs["dim 3"] = "defocus idx"
    dataset.attrs["dim 4"] = "probe idx"
    dataset.attrs["dim 5"] = "k_y idx"
    dataset.attrs["dim 6"] = "k_x idx"
    dataset.attrs["units"] = "dimensionless"

    group.file.close()

    return None



def _num_atomic_configs_in_subset(sim_params, atomic_config_subset_idx):
    stem_system_model_params = \
        sim_params.core_attrs["stem_system_model_params"]
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]
    kwargs = \
        {"sample_specification": sample_specification,
         "subset_idx": atomic_config_subset_idx}
    num_atomic_configs_in_subset = \
        prismatique.sample._num_frozen_phonon_configs_in_subset(**kwargs)

    return num_atomic_configs_in_subset



def _initialize_pyprismatic_sim_obj(sim_params):
    pyprismatic_sim_obj = pyprismatic.Metadata()

    stem_system_model_params = sim_params.core_attrs["stem_system_model_params"]
    output_params = sim_params.core_attrs["output_params"]
    worker_params = sim_params.core_attrs["worker_params"]

    base_output_params = output_params.core_attrs["base_params"]
    output_dirname = base_output_params.core_attrs["output_dirname"]

    mod_alias = prismatique.sample
    func_alias = mod_alias._set_pyprismatic_sim_obj_attrs_to_default_values
    func_alias(pyprismatic_sim_obj, output_dirname)

    _ = _generate_probe_positions(sim_params, save=True)

    stem_system_model_params = sim_params.core_attrs["stem_system_model_params"]
    output_params = sim_params.core_attrs["output_params"]
    worker_params = sim_params.core_attrs["worker_params"]

    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]
    
    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "stem_system_model_params": stem_system_model_params,
              "output_dirname": output_dirname}
    _unpack_stem_system_model_params_into_pyprismatic_sim_obj(**kwargs)

    unpack_stem_output_params_into_pyprismatic_sim_obj = \
        prismatique.sample._unpack_stem_output_params_into_pyprismatic_sim_obj
    unpack_worker_params_into_pyprismatic_sim_obj = \
        prismatique.sample._unpack_worker_params_into_pyprismatic_sim_obj

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "stem_output_params": output_params,
              "sample_specification": sample_specification}
    unpack_stem_output_params_into_pyprismatic_sim_obj(**kwargs)

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "worker_params": worker_params}
    unpack_worker_params_into_pyprismatic_sim_obj(**kwargs)

    return pyprismatic_sim_obj



def _unpack_stem_system_model_params_into_pyprismatic_sim_obj(
        pyprismatic_sim_obj, stem_system_model_params, output_dirname):
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]
    probe_model_params = \
        stem_system_model_params.core_attrs["probe_model_params"]
    specimen_tilt = \
        stem_system_model_params.core_attrs["specimen_tilt"]

    unpack_sample_specification_into_pyprismatic_sim_obj = \
        prismatique.sample._unpack_sample_specification_into_pyprismatic_sim_obj
    unpack_probe_model_params_into_pyprismatic_sim_obj = \
        prismatique.sample._unpack_probe_model_params_into_pyprismatic_sim_obj

    unpack_sample_specification_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                         output_dirname,
                                                         sample_specification)
    unpack_probe_model_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                       probe_model_params,
                                                       output_dirname)

    pyprismatic_sim_obj.probeXtilt = specimen_tilt[0]
    pyprismatic_sim_obj.probeYtilt = specimen_tilt[1]

    pyprismatic_sim_obj.probes_file = \
        prismatique.sample._generate_temp_scan_config_filename(output_dirname)

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

    num_atomic_config_subsets = len(rng_seeds)
    if atomic_config_subset_idx == num_atomic_config_subsets-1:
        kwargs = {"sim_params": sim_params,
                  "filename": _intensity_output_filename(sim_params)}
        _calc_and_write_remaining_intensity_data_to_output_file(**kwargs)
        _remove_any_temp_intensity_4d_stem_data_in_output_file(**kwargs)

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
    stem_system_model_params = \
        sim_params.core_attrs["stem_system_model_params"]
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]

    output_param_subset = _output_param_subset(sim_params)
    defocii = _defocii(sim_params)  # In Å.

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "sample_specification": sample_specification,
              "output_param_subset": output_param_subset,
              "subset_idx": atomic_config_subset_idx,
              "defocus_idx": defocus_idx,
              "defocii": defocii,
              "rng_seeds": rng_seeds}
    mod_alias = prismatique.sample
    mod_alias._update_pyprismatic_sim_obj_for_next_prismatic_sim(**kwargs)

    return None



def _call_pyprismatic_sim_obj_go(pyprismatic_sim_obj):
    pyprismatic_sim_obj.go()

    return None



def _postprocess_and_reorganize_current_prismatic_sim_output(
        pyprismatic_sim_obj, sim_params, atomic_config_subset_idx, defocus_idx):
    if defocus_idx == 0:
        output_param_subset = _output_param_subset(sim_params)
        func_alias = _move_sample_specification_output_to_separate_file
        kwargs = {"output_dirname": output_param_subset["output_dirname"]}
        if pyprismatic_sim_obj.savePotentialSlices:
            kwargs["sample_specification_type"] = "potential_slice_subset"
            func_alias(**kwargs)
            gc.collect()
        if pyprismatic_sim_obj.saveSMatrix:
            kwargs["sample_specification_type"] = "S_matrix_subset"
            func_alias(**kwargs)
            gc.collect()
    
    kwargs = {"sim_params": sim_params,
              "atomic_config_subset_idx": atomic_config_subset_idx}
    num_atomic_configs_in_subset = _num_atomic_configs_in_subset(**kwargs)
    output_layer_depths = _output_layer_depths(sim_params)

    for output_layer_idx, _ in enumerate(output_layer_depths):
        for atomic_config_idx in range(num_atomic_configs_in_subset):
            quartet_of_indices = \
                {"atomic_config_subset_idx": atomic_config_subset_idx,
                 "output_layer_idx": output_layer_idx,
                 "atomic_config_idx": atomic_config_idx,
                 "defocus_idx": defocus_idx}
            
            kwargs = {"sim_params": sim_params,
                      "quartet_of_indices": quartet_of_indices}
            func_alias = _postprocess_and_reorganize_dp_subset
            func_alias(**kwargs)
            gc.collect()

    return None



def _move_sample_specification_output_to_separate_file(
        output_dirname, sample_specification_type):
    mod_alias = prismatique.sample
    kwargs = {"new_output_filename": "",
              "output_dirname": output_dirname,
              "sample_specification_type": sample_specification_type}
    mod_alias._move_sample_specification_output_to_separate_file(**kwargs)
    
    return None



def _postprocess_and_reorganize_dp_subset(sim_params, quartet_of_indices):
    kwargs = \
        {"sim_params": sim_params,
         "quartet_of_indices": quartet_of_indices}
    unprocessed_complex_dp_subset_signal = \
        _load_unprocessed_complex_dp_subset_signal(**kwargs)

    if _wavefunction_output_is_to_be_saved(sim_params):
        output_data = unprocessed_complex_dp_subset_signal.data
        filenames = _wavefunction_output_filenames(sim_params)
        filename = filenames[quartet_of_indices["atomic_config_subset_idx"]]
        kwargs = {"sim_params": sim_params,
                  "unprocessed_complex_dp_subset": output_data,
                  "filename": filename,
                  "quartet_of_indices": quartet_of_indices}
        _write_unprocessed_complex_dp_subset_to_output_file(**kwargs)

    kwargs = \
        {"unprocessed_intensity_dp_subset_signal": \
         empix.abs_sq(unprocessed_complex_dp_subset_signal),
         "sim_params": sim_params}
    postprocessed_intensity_dp_subset_signal = \
         _postprocess_intensity_dp_subset_signal(**kwargs)

    kwargs = {"sim_params": sim_params,
              "filename": _intensity_output_filename(sim_params),
              "quartet_of_indices": quartet_of_indices,
              "postprocessed_intensity_dp_subset_signal": \
              postprocessed_intensity_dp_subset_signal}
    _update_intensity_4d_stem_data_in_output_file(**kwargs)

    return None



def _load_unprocessed_complex_dp_subset_signal(sim_params, quartet_of_indices):
    output_dirname = _output_param_subset(sim_params)["output_dirname"]
    filename = output_dirname + "/prismatic_output.h5"

    output_layer_idx = quartet_of_indices["output_layer_idx"]
    atomic_config_idx = quartet_of_indices["atomic_config_idx"]
    atomic_config_subset_idx = quartet_of_indices["atomic_config_subset_idx"]

    output_layer_idx_str = str(output_layer_idx).rjust(4, "0")
    atomic_config_idx_str = str(atomic_config_idx).rjust(4, "0")

    unformatted_path_in_file = ("4DSTEM_simulation/data/datacubes"
                                "/CBED_array_depth{}_fp{}/data")
    path_in_file = unformatted_path_in_file.format(output_layer_idx_str,
                                                   atomic_config_idx_str)
    dataset_id = h5pywrappers.obj.ID(filename, path_in_file)
    dataset = h5pywrappers.dataset.load(dataset_id, read_only=True)

    kwargs = {"sim_params": sim_params,
              "atomic_config_subset_idx": atomic_config_subset_idx}
    num_atomic_configs_in_subset = _num_atomic_configs_in_subset(**kwargs)
    
    unprocessed_complex_dp_subset_data = \
        np.transpose(dataset[:, 0, :, :], axes=(0, 2, 1))[:, ::-1, :]
    unprocessed_complex_dp_subset_data *= \
        num_atomic_configs_in_subset  # Fixes prismatic bug.

    dataset.file.close()

    kwargs = \
        {"sim_params": sim_params,
         "navigation_dims": (unprocessed_complex_dp_subset_data.shape[0],),
         "signal_dtype": "complex"}
    unprocessed_complex_dp_subset_signal = \
        _blank_unprocessed_dp_set_signal(**kwargs)
    unprocessed_complex_dp_subset_signal.data = \
        unprocessed_complex_dp_subset_data

    return unprocessed_complex_dp_subset_signal



def _write_unprocessed_complex_dp_subset_to_output_file(
        sim_params,
        unprocessed_complex_dp_subset,
        filename,
        quartet_of_indices):
    path_in_file = "/data/4D_STEM/complex_valued_DPs"
    dataset_id = h5pywrappers.obj.ID(filename, path_in_file)
    
    multi_dim_slice = (quartet_of_indices["output_layer_idx"],
                       quartet_of_indices["atomic_config_idx"],
                       quartet_of_indices["defocus_idx"],
                       slice(None),
                       slice(None),
                       slice(None))
    datasubset_id = h5pywrappers.datasubset.ID(dataset_id, multi_dim_slice)

    datasubset = unprocessed_complex_dp_subset
    h5pywrappers.datasubset.save(datasubset, datasubset_id)

    return None



def _postprocess_intensity_dp_subset_signal(
        unprocessed_intensity_dp_subset_signal, sim_params):
    output_params = sim_params.core_attrs["output_params"]
    base_output_params = output_params.core_attrs["base_params"]
    cbed_params = base_output_params.core_attrs["cbed_params"]
    postprocessing_seq = cbed_params.core_attrs["postprocessing_seq"]

    kwargs = \
        {"input_signal": unprocessed_intensity_dp_subset_signal,
         "postprocessing_seq": postprocessing_seq}
    postprocessed_intensity_dp_subset_signal = \
        prismatique._signal._postprocess_2d_signal(**kwargs)

    return postprocessed_intensity_dp_subset_signal



def _update_intensity_4d_stem_data_in_output_file(
        sim_params,
        filename,
        quartet_of_indices,
        postprocessed_intensity_dp_subset_signal):
    output_layer_idx = quartet_of_indices["output_layer_idx"]
    atomic_config_idx = quartet_of_indices["atomic_config_idx"]
    atomic_config_subset_idx = quartet_of_indices["atomic_config_subset_idx"]
    defocus_idx = quartet_of_indices["defocus_idx"]

    path_in_file = "/data/4D_STEM/intensity_DPs"
    dataset_id = h5pywrappers.obj.ID(filename, path_in_file)
    dataset = h5pywrappers.dataset.load(dataset_id, read_only=False)

    dataset[output_layer_idx] += (postprocessed_intensity_dp_subset_signal.data
                                  * _w_f_l(sim_params, l=defocus_idx)
                                  / _total_num_frozen_phonon_configs(sim_params)
                                  / np.sqrt(np.pi))

    kwargs = {"sim_params": sim_params,
              "atomic_config_subset_idx": atomic_config_subset_idx}
    num_atomic_configs_in_subset = _num_atomic_configs_in_subset(**kwargs)
    num_atomic_config_subsets = _num_atomic_config_subsets(sim_params)
    num_defocii = len(_defocii(sim_params))

    output_params = sim_params.core_attrs["output_params"]
    base_output_params = output_params.core_attrs["base_params"]
    cbed_params = base_output_params.core_attrs["cbed_params"]

    if ((atomic_config_subset_idx == num_atomic_config_subsets-1)
        and (atomic_config_idx == num_atomic_configs_in_subset-1)
        and (defocus_idx == num_defocii-1)):
        avg_num_electrons_per_postprocessed_dp = \
            cbed_params.core_attrs["avg_num_electrons_per_postprocessed_dp"]

        datasubset = dataset[output_layer_idx].clip(min=0)
        datasubset *= (avg_num_electrons_per_postprocessed_dp
                       / np.mean(np.sum(datasubset, axis=(-2, -1))))
        if cbed_params.core_attrs["apply_shot_noise"]:
            datasubset = np.random.poisson(datasubset)
        dataset[output_layer_idx] = datasubset

    dataset.file.close()

    return None



def _w_f_l(sim_params, l):
    stem_system_model_params = \
        sim_params.core_attrs["stem_system_model_params"]
    probe_model_params = \
        stem_system_model_params.core_attrs["probe_model_params"]

    w_f_l = probe_model_params._gauss_hermite_weights[l]

    return w_f_l



def _total_num_frozen_phonon_configs(sim_params):
    stem_system_model_params = \
        sim_params.core_attrs["stem_system_model_params"]
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]

    kwargs = \
        {"sample_specification": sample_specification}
    total_num_frozen_phonon_configs = \
        prismatique.sample._total_num_frozen_phonon_configs(**kwargs)
    
    return total_num_frozen_phonon_configs



def _calc_and_write_remaining_intensity_data_to_output_file(sim_params,
                                                            filename):
    func_alias_1 = \
        _calc_and_write_com_data_of_output_layer_to_output_file
    func_alias_2 = \
        _calc_and_write_intensity_3d_stem_data_of_output_layer_to_output_file
    func_alias_3 = \
        _calc_and_write_intensity_2d_stem_data_of_output_layer_to_output_file

    output_layer_depths = _output_layer_depths(sim_params)
    for output_layer_idx, _ in enumerate(output_layer_depths):
        kwargs = {"sim_params": sim_params,
                  "output_layer_idx": output_layer_idx,
                  "filename": filename}
        
        func_aliases = []
        if _com_output_is_to_be_saved(sim_params):
            func_aliases.append(func_alias_1)
        if _intensity_3d_stem_output_is_to_be_saved(sim_params):
            func_aliases.append(func_alias_2)
        if _intensity_2d_stem_output_is_to_be_saved(sim_params):
            func_aliases.append(func_alias_3)

        for func_alias in func_aliases:
            func_alias(**kwargs)
            gc.collect()

    return None



def _intensity_3d_stem_output_is_to_be_saved(sim_params):
    output_params = \
        sim_params.core_attrs["output_params"]
    base_output_params = \
        output_params.core_attrs["base_params"]
    radial_step_size_for_3d_stem = \
        base_output_params.core_attrs["radial_step_size_for_3d_stem"]

    if radial_step_size_for_3d_stem > 0:
        intensity_3d_stem_output_is_to_be_saved = True
    else:
        intensity_3d_stem_output_is_to_be_saved = False

    return intensity_3d_stem_output_is_to_be_saved



def _calc_and_write_com_data_of_output_layer_to_output_file(sim_params,
                                                            output_layer_idx,
                                                            filename):
    k_mesh = _k_mesh(sim_params, for_postprocessed_dp=True)
    
    intensity_4d_stem_signal = _load_intensity_4d_stem_signal(sim_params,
                                                              output_layer_idx,
                                                              filename)

    kwargs = {"center": (0, 0), "radial_range": None}
    optional_params = empix.OptionalAnnularIntegrationParams(**kwargs)

    normalization_signal = empix.annularly_integrate(intensity_4d_stem_signal,
                                                     optional_params)
    
    integrand_signal = copy.deepcopy(intensity_4d_stem_signal)
    integrand_signal.data *= k_mesh[0]
    com_x_signal = empix.annularly_integrate(integrand_signal, optional_params)
    com_x_signal.data /= normalization_signal.data

    integrand_signal = copy.deepcopy(intensity_4d_stem_signal)
    integrand_signal.data *= k_mesh[1]
    com_y_signal = empix.annularly_integrate(integrand_signal, optional_params)
    com_y_signal.data /= normalization_signal.data

    com_data_of_output_layer = np.array([com_x_signal.data, com_y_signal.data])

    path_in_file = "/data/center_of_mass_momentum"
    dataset_id = h5pywrappers.obj.ID(filename, path_in_file)
    
    multi_dim_slice = (output_layer_idx, slice(None), slice(None))
    datasubset_id = h5pywrappers.datasubset.ID(dataset_id, multi_dim_slice)

    datasubset = com_data_of_output_layer
    h5pywrappers.datasubset.save(datasubset, datasubset_id)

    return None



def _k_mesh(sim_params, for_postprocessed_dp):
    k_x = _k_x(sim_params, for_postprocessed_dp)
    k_y = _k_y(sim_params, for_postprocessed_dp)
    k_mesh = np.meshgrid(k_x, k_y, indexing="xy")

    return k_mesh



def _load_intensity_4d_stem_signal(sim_params, output_layer_idx, filename):
    path_in_file = "/data/4D_STEM/intensity_DPs"
    dataset_id = h5pywrappers.obj.ID(filename, path_in_file)
    
    multi_dim_slice = (output_layer_idx, slice(None), slice(None), slice(None))
    datasubset_id = h5pywrappers.datasubset.ID(dataset_id, multi_dim_slice)
    intensity_4d_stem_data = h5pywrappers.datasubset.load(datasubset_id)

    kwargs = {"sim_params": sim_params,
              "navigation_dims": (intensity_4d_stem_data.shape[0],),
              "signal_dtype": "float"}
    intensity_4d_stem_signal = _blank_postprocessed_dp_set_signal(**kwargs)
    intensity_4d_stem_signal.data = intensity_4d_stem_data

    return intensity_4d_stem_signal



def _calc_and_write_intensity_3d_stem_data_of_output_layer_to_output_file(
        sim_params, output_layer_idx, filename):
    intensity_4d_stem_signal = \
        _load_intensity_4d_stem_signal(sim_params, output_layer_idx, filename)
    intensity_3d_stem_signal = \
        _integrate_dp_to_3d_stem_signal(intensity_4d_stem_signal, sim_params)
    intensity_3d_stem_data_of_output_layer = \
        intensity_3d_stem_signal.data

    path_in_file = "/data/3D_STEM/integrated_intensities"
    dataset_id = h5pywrappers.obj.ID(filename, path_in_file)
    
    multi_dim_slice = (output_layer_idx, slice(None), slice(None))
    datasubset_id = h5pywrappers.datasubset.ID(dataset_id, multi_dim_slice)

    datasubset = intensity_3d_stem_data_of_output_layer
    h5pywrappers.datasubset.save(datasubset, datasubset_id)

    return None



def _calc_and_write_intensity_2d_stem_data_of_output_layer_to_output_file(
        sim_params, output_layer_idx, filename):
    output_params = \
        sim_params.core_attrs["output_params"]
    base_output_params = \
        output_params.core_attrs["base_params"]
    radial_angular_range_for_2d_stem = \
        np.array(base_output_params.core_attrs["radial_range_for_2d_stem"])
    wavelength = \
        _wavelength(sim_params)
    radial_k_range_for_2d_stem = \
        (radial_angular_range_for_2d_stem / 1000) / wavelength
    
    intensity_4d_stem_signal = _load_intensity_4d_stem_signal(sim_params,
                                                              output_layer_idx,
                                                              filename)

    input_signal = intensity_4d_stem_signal
    pixel_area = np.abs(input_signal.axes_manager[-2].scale
                        * input_signal.axes_manager[-1].scale)
    input_signal.data /= pixel_area

    kwargs = {"center": (0, 0), "radial_range": radial_k_range_for_2d_stem}
    optional_params = empix.OptionalAnnularIntegrationParams(**kwargs)
    intensity_2d_stem_signal = empix.annularly_integrate(input_signal,
                                                         optional_params)
    
    intensity_2d_stem_data_of_output_layer = intensity_2d_stem_signal.data

    path_in_file = "/data/2D_STEM/integrated_intensities"
    dataset_id = h5pywrappers.obj.ID(filename, path_in_file)
    
    multi_dim_slice = (output_layer_idx, slice(None))
    datasubset_id = h5pywrappers.datasubset.ID(dataset_id, multi_dim_slice)

    datasubset = intensity_2d_stem_data_of_output_layer
    h5pywrappers.datasubset.save(datasubset, datasubset_id)

    return None



def _remove_any_temp_intensity_4d_stem_data_in_output_file(sim_params,
                                                           filename):
    output_params = sim_params.core_attrs["output_params"]
    base_output_params = output_params.core_attrs["base_params"]
    cbed_params = base_output_params.core_attrs["cbed_params"]

    intensity_4d_stem_data_is_temp = \
        not cbed_params.core_attrs["save_final_intensity"]

    if intensity_4d_stem_data_is_temp:
        path_in_file = "/data/4D_STEM"
        group_id = h5pywrappers.obj.ID(filename, path_in_file)
        group = h5pywrappers.group.load(group_id, read_only=False)
        del group["intensity_DPs"]
        group.file.close()

    return None



def _serialize_sim_params(sim_params):
    output_params = sim_params.core_attrs["output_params"]
    base_output_params = output_params.core_attrs["base_params"]    
    output_dirname = base_output_params.core_attrs["output_dirname"]
    
    sim_params.dump(filename=output_dirname+"/stem_sim_params.json",
                    overwrite=True)

    return None
    


###########################
## Define error messages ##
###########################

_check_data_size_err_msg_1 = \
    ("The data size of the output of the STEM simulation to be performed is {} "
     "bytes, exceeding the maximum allowed size of {} bytes specified by the "
     "simulation parameters.")
