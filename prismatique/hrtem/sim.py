r"""For running HRTEM simulations.

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



def _check_and_convert_hrtem_system_model_params(ctor_params):
    check_and_convert_hrtem_system_model_params = \
        prismatique.hrtem.system._check_and_convert_hrtem_system_model_params
    hrtem_system_model_params = \
        check_and_convert_hrtem_system_model_params(ctor_params)
    
    return hrtem_system_model_params



def _pre_serialize_hrtem_system_model_params(hrtem_system_model_params):
    pre_serialize_hrtem_system_model_params = \
        prismatique.hrtem.system._pre_serialize_hrtem_system_model_params
    serializable_rep = \
        pre_serialize_hrtem_system_model_params(hrtem_system_model_params)

    return serializable_rep



def _de_pre_serialize_hrtem_system_model_params(serializable_rep):
    de_pre_serialize_hrtem_system_model_params = \
        prismatique.hrtem.system._de_pre_serialize_hrtem_system_model_params
    hrtem_system_model_params = \
        de_pre_serialize_hrtem_system_model_params(serializable_rep)

    return hrtem_system_model_params



def _check_and_convert_output_params(ctor_params):
    check_and_convert_output_params = \
        prismatique.hrtem.output._check_and_convert_output_params
    output_params = \
        check_and_convert_output_params(ctor_params)
    
    return output_params



def _pre_serialize_output_params(output_params):
    pre_serialize_output_params = \
        prismatique.hrtem.output._pre_serialize_output_params
    serializable_rep = \
        pre_serialize_output_params(output_params)

    return serializable_rep



def _de_pre_serialize_output_params(serializable_rep):
    de_pre_serialize_output_params = \
        prismatique.hrtem.output._de_pre_serialize_output_params
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
        {"hrtem_system_model_params": \
         _check_and_convert_hrtem_system_model_params,
         "output_params": _check_and_convert_output_params,
         "worker_params": _check_and_convert_worker_params}

    _pre_serialization_funcs = \
        {"hrtem_system_model_params": _pre_serialize_hrtem_system_model_params,
         "output_params": _pre_serialize_output_params,
         "worker_params": _pre_serialize_worker_params}

    _de_pre_serialization_funcs = \
        {"hrtem_system_model_params": \
         _de_pre_serialize_hrtem_system_model_params,
         "output_params": _de_pre_serialize_output_params,
         "worker_params": _de_pre_serialize_worker_params}
    
    def __init__(self,
                 hrtem_system_model_params,
                 output_params=None,
                 worker_params=None):
        ctor_params = {"hrtem_system_model_params": hrtem_system_model_params,
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
    _check_sim_params(sim_params)
    _run(sim_params)
    _serialize_sim_params(sim_params)
    print("\n\n\n")

    return None



def _check_sim_params(sim_params):
    temp_ctor_params = {"sim_params": sim_params}
    sim_params = _check_and_convert_sim_params(temp_ctor_params)

    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    
    mod_alias = prismatique.hrtem.output
    mod_alias._check_hrtem_system_model_params(hrtem_system_model_params)
    
    _pre_save(sim_params)
    _check_data_size(sim_params)

    return None



def _pre_save(sim_params):
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    sample_specification = \
        hrtem_system_model_params.core_attrs["sample_specification"]
    
    output_params = sim_params.core_attrs["output_params"]
    output_dirname = output_params.core_attrs["output_dirname"]

    if output_params.core_attrs["save_potential_slices"]:
        unformatted_basename = "potential_slices_of_subset_{}.h5"
        prismatique.sample._pre_save(sample_specification,
                                     output_dirname,
                                     unformatted_basename)

    filenames = []
    if _intensity_output_is_to_be_saved(sim_params):
        filenames.append(_intensity_output_filename(sim_params))
    if _wavefunction_output_is_to_be_saved(sim_params):
        filenames += _wavefunction_output_filenames(sim_params)

    for filename in filenames:
        if pathlib.Path(filename).is_file():
            pathlib.Path(filename).unlink(missing_ok=True)

    prismatique.sample._check_hdf5_filenames(filenames)

    return None



def _intensity_output_is_to_be_saved(sim_params):
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    output_params = \
        sim_params.core_attrs["output_params"]
    image_params = \
        output_params.core_attrs["image_params"]
    intensity_output_is_to_be_saved = \
        image_params.core_attrs["save_final_intensity"]

    return intensity_output_is_to_be_saved



def _intensity_output_filename(sim_params):
    output_dirname = _output_param_subset(sim_params)["output_dirname"]
    filename = output_dirname + "/hrtem_sim_intensity_output.h5"

    return filename



def _wavefunction_output_is_to_be_saved(sim_params):
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    output_params = \
        sim_params.core_attrs["output_params"]
    image_params = \
        output_params.core_attrs["image_params"]
    wavefunction_output_is_to_be_saved = \
        image_params.core_attrs["save_wavefunctions"]

    return wavefunction_output_is_to_be_saved



def _wavefunction_output_filenames(sim_params):
    output_params = sim_params.core_attrs["output_params"]
    output_dirname = output_params.core_attrs["output_dirname"]

    num_atomic_config_subsets = _num_atomic_config_subsets(sim_params)

    filenames = []
    for atomic_config_subset_idx in range(num_atomic_config_subsets):
        unformatted_basename = "hrtem_sim_wavefunction_output_of_subset_{}.h5"
        basename = unformatted_basename.format(atomic_config_subset_idx)
        filenames.append(output_dirname + "/" + basename)

    return filenames



def _num_atomic_config_subsets(sim_params):
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    sample_specification = \
        hrtem_system_model_params.core_attrs["sample_specification"]

    kwargs = \
        {"sample_specification": sample_specification}
    num_atomic_config_subsets = \
        prismatique.sample._num_frozen_phonon_config_subsets(**kwargs)

    return num_atomic_config_subsets



def _check_data_size(sim_params):
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]

    output_params = sim_params.core_attrs["output_params"]
    max_data_size = output_params.core_attrs["max_data_size"]

    kwargs = {"hrtem_system_model_params": hrtem_system_model_params,
              "output_params": output_params}
    output_data_size = prismatique.hrtem.output._data_size(**kwargs)

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
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    sample_specification = \
        hrtem_system_model_params.core_attrs["sample_specification"]
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
    
    output_dirname = output_params.core_attrs["output_dirname"]
    save_potential_slices = output_params.core_attrs["save_potential_slices"]

    output_param_subset = {"output_dirname": output_dirname,
                           "save_potential_slices": save_potential_slices,
                           "save_S_matrices": False}

    return output_param_subset



def _initialize_output_files(sim_params):
    _write_metadata_to_output_files(sim_params)
    _initialize_data_in_output_files(sim_params)

    return None



def _write_metadata_to_output_files(sim_params):
    filenames = []
    idx_offset = 0
    if _intensity_output_is_to_be_saved(sim_params):
        filenames.append(_intensity_output_filename(sim_params))
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
    h5pywrappers.group.save(None, group_id, write_mode="a")
    group = h5pywrappers.group.load(group_id, read_only=False)

    if filename == _intensity_output_filename(sim_params):
        for_postprocessed_image = True
    else:
        for_postprocessed_image = False

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
    if for_postprocessed_image:
        image_set_signal = _blank_postprocessed_image_set_signal(**kwargs)
    else:
        image_set_signal = _blank_unprocessed_image_set_signal(**kwargs)

    offset = image_set_signal.axes_manager[0].offset
    size = image_set_signal.axes_manager[0].size
    scale = image_set_signal.axes_manager[0].scale

    r_x = offset + scale*np.arange(size)

    return r_x



def _r_y(sim_params, for_postprocessed_image):
    kwargs = {"sim_params": sim_params,
              "navigation_dims": tuple(),
              "signal_dtype": "float"}
    if for_postprocessed_image:
        image_set_signal = _blank_postprocessed_image_set_signal(**kwargs)
    else:
        image_set_signal = _blank_unprocessed_image_set_signal(**kwargs)

    offset = image_set_signal.axes_manager[1].offset
    size = image_set_signal.axes_manager[1].size
    scale = image_set_signal.axes_manager[1].scale

    r_y = offset + scale*np.arange(size)

    return r_y



def _blank_postprocessed_image_set_signal(sim_params,
                                          navigation_dims,
                                          signal_dtype):
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    output_params = \
        sim_params.core_attrs["output_params"]

    sample_specification = \
        hrtem_system_model_params.core_attrs["sample_specification"]

    image_params = output_params.core_attrs["image_params"]
    postprocessing_seq = image_params.core_attrs["postprocessing_seq"]

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
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    sample_specification = \
        hrtem_system_model_params.core_attrs["sample_specification"]

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
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    sample_specification = \
        hrtem_system_model_params.core_attrs["sample_specification"]
    tilt_params = \
        hrtem_system_model_params.core_attrs["tilt_params"]
    gun_model_params = \
        hrtem_system_model_params.core_attrs["gun_model_params"]
    mean_beam_energy = \
        gun_model_params.core_attrs["mean_beam_energy"]
    
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



def _probe_model_params(sim_params):
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    gun_model_params = \
        hrtem_system_model_params.core_attrs["gun_model_params"]
    lens_model_params = \
        hrtem_system_model_params.core_attrs["lens_model_params"]
    defocal_offset_supersampling = \
        hrtem_system_model_params.core_attrs["defocal_offset_supersampling"]

    kwargs = {"lens_model_params": lens_model_params,
              "gun_model_params": gun_model_params,
              "defocal_offset_supersampling": defocal_offset_supersampling}
    probe_model_params = embeam.stem.probe.ModelParams(**kwargs)

    return probe_model_params



def _wavelength(sim_params):
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    gun_model_params = \
        hrtem_system_model_params.core_attrs["gun_model_params"]
    mean_beam_energy = \
        gun_model_params.core_attrs["mean_beam_energy"]
    wavelength = \
        embeam.wavelength(mean_beam_energy)

    return wavelength



def _initialize_data_in_output_files(sim_params):
    filenames = []
    idx_offset = -1
    if _intensity_output_is_to_be_saved(sim_params):
        filenames.append(_intensity_output_filename(sim_params))
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
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    sample_specification = \
        hrtem_system_model_params.core_attrs["sample_specification"]
    kwargs = \
        {"sample_specification": sample_specification,
         "subset_idx": atomic_config_subset_idx}
    num_atomic_configs_in_subset = \
        prismatique.sample._num_frozen_phonon_configs_in_subset(**kwargs)

    return num_atomic_configs_in_subset



def _initialize_pyprismatic_sim_obj(sim_params):
    pyprismatic_sim_obj = pyprismatic.Metadata()

    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    output_params = \
        sim_params.core_attrs["output_params"]
    worker_params = \
        sim_params.core_attrs["worker_params"]

    output_dirname = output_params.core_attrs["output_dirname"]

    mod_alias = prismatique.sample
    func_alias = mod_alias._set_pyprismatic_sim_obj_attrs_to_default_values
    func_alias(pyprismatic_sim_obj, output_dirname)
    
    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "hrtem_system_model_params": hrtem_system_model_params,
              "output_dirname": output_dirname}
    _unpack_hrtem_system_model_params_into_pyprismatic_sim_obj(**kwargs)

    unpack_hrtem_output_params_into_pyprismatic_sim_obj = \
        prismatique.sample._unpack_hrtem_output_params_into_pyprismatic_sim_obj
    unpack_worker_params_into_pyprismatic_sim_obj = \
        prismatique.sample._unpack_worker_params_into_pyprismatic_sim_obj

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "hrtem_output_params": output_params}
    unpack_hrtem_output_params_into_pyprismatic_sim_obj(**kwargs)

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "worker_params": worker_params}
    unpack_worker_params_into_pyprismatic_sim_obj(**kwargs)

    return pyprismatic_sim_obj



def _unpack_hrtem_system_model_params_into_pyprismatic_sim_obj(
        pyprismatic_sim_obj, hrtem_system_model_params, output_dirname):
    sample_specification = \
        hrtem_system_model_params.core_attrs["sample_specification"]
    gun_model_params = \
        hrtem_system_model_params.core_attrs["gun_model_params"]
    lens_model_params = \
        hrtem_system_model_params.core_attrs["lens_model_params"]
    tilt_params = \
        hrtem_system_model_params.core_attrs["tilt_params"]

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
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    sample_specification = \
        hrtem_system_model_params.core_attrs["sample_specification"]

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
    move_sample_specification_output_to_separate_file = \
        prismatique.stem.sim._move_sample_specification_output_to_separate_file
    
    if defocus_idx == 0:
        output_param_subset = _output_param_subset(sim_params)
        func_alias = move_sample_specification_output_to_separate_file
        kwargs = {"output_dirname": output_param_subset["output_dirname"]}
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
             "sim_params": sim_params}
        postprocessed_intensity_image_subset_signal = \
            _postprocess_intensity_image_subset_signal(**kwargs)

        kwargs = {"sim_params": sim_params,
                  "filename": _intensity_output_filename(sim_params),
                  "triplet_of_indices": triplet_of_indices,
                  "postprocessed_intensity_image_subset_signal": \
                  postprocessed_intensity_image_subset_signal}
        _update_intensity_data_in_output_file(**kwargs)

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
    multi_dim_slice = tuple([slice(None)]*(signal_rank-2)
                            + [slice(None, None, -1), slice(None)])

    temp_data = \
        np.fft.fft2(signal_data[multi_dim_slice]) * angular_mask
    unprocessed_complex_image_subset_signal.data = \
        np.fft.ifft2(temp_data)[multi_dim_slice]

    return None



def _angular_mask(sim_params):
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    objective_aperture_params = \
        hrtem_system_model_params.core_attrs["objective_aperture_params"]
    offset = \
        objective_aperture_params.core_attrs["offset"]
    window = \
        objective_aperture_params.core_attrs["window"]
    
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
    output_params = sim_params.core_attrs["output_params"]
    image_params = output_params.core_attrs["image_params"]
    postprocessing_seq = image_params.core_attrs["postprocessing_seq"]

    kwargs = \
        {"input_signal": unprocessed_intensity_image_subset_signal,
         "postprocessing_seq": postprocessing_seq}
    postprocessed_intensity_image_subset_signal = \
        prismatique._signal._postprocess_2d_signal(**kwargs)

    return postprocessed_intensity_image_subset_signal



def _update_intensity_data_in_output_file(
        sim_params,
        filename,
        triplet_of_indices,
        postprocessed_intensity_image_subset_signal):
    atomic_config_idx = triplet_of_indices["atomic_config_idx"]
    atomic_config_subset_idx = triplet_of_indices["atomic_config_subset_idx"]
    defocus_idx = triplet_of_indices["defocus_idx"]

    signal_data = postprocessed_intensity_image_subset_signal.data
    tilt_weights = _tilt_weights(sim_params)

    path_in_file = "/data/intensity_image"
    dataset_id = h5pywrappers.obj.ID(filename, path_in_file)
    dataset = h5pywrappers.dataset.load(dataset_id, read_only=False)

    dataset[()] += (np.einsum('ijk,i->jk', signal_data, tilt_weights)
                    * _w_f_l(sim_params, l=defocus_idx)
                    / _total_num_frozen_phonon_configs(sim_params)
                    / np.sqrt(np.pi))

    kwargs = {"sim_params": sim_params,
              "atomic_config_subset_idx": atomic_config_subset_idx}
    num_atomic_configs_in_subset = _num_atomic_configs_in_subset(**kwargs)
    num_atomic_config_subsets = _num_atomic_config_subsets(sim_params)
    num_defocii = len(_defocii(sim_params))

    output_params = sim_params.core_attrs["output_params"]
    image_params = output_params.core_attrs["image_params"]

    if ((atomic_config_subset_idx == num_atomic_config_subsets-1)
        and (atomic_config_idx == num_atomic_configs_in_subset-1)
        and (defocus_idx == num_defocii-1)):
        avg_num_electrons_per_postprocessed_image = \
            image_params.core_attrs["avg_num_electrons_per_postprocessed_image"]
        
        datasubset = dataset[()].clip(min=0)
        datasubset *= (avg_num_electrons_per_postprocessed_image
                       / np.sum(datasubset, axis=(-2, -1)))
        if image_params.core_attrs["apply_shot_noise"]:
            datasubset = np.random.poisson(datasubset)
        dataset[()] = datasubset

    dataset.file.close()

    return None



def _tilt_weights(sim_params):
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    
    tilt_params = hrtem_system_model_params.core_attrs["tilt_params"]
    tilt_offset = tilt_params.core_attrs["offset"]  # In mrads.
    tilt_spread = tilt_params.core_attrs["spread"]  # In mrads.
    tilt_series = _tilt_series(sim_params)  # In mrads.
    rel_tilts = np.linalg.norm(tilt_series - tilt_offset, axis=1)  # In mrads.

    if tilt_spread > 0:
        exp_arg = -0.5*(rel_tilts/tilt_spread)*(rel_tilts/tilt_spread)
        tilt_weights = np.exp(exp_arg)
    else:
        tilt_weights = (rel_tilts == 0.0).astype(float)
    tilt_weights /= np.linalg.norm(tilt_weights)

    return tilt_weights



def _w_f_l(sim_params, l):
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]

    # Build probe model simply to extract HRTEM beam defocus weights more
    # easily.
    probe_model_params = _probe_model_params(sim_params)

    w_f_l = probe_model_params._gauss_hermite_weights[l]

    return w_f_l



def _total_num_frozen_phonon_configs(sim_params):
    hrtem_system_model_params = \
        sim_params.core_attrs["hrtem_system_model_params"]
    sample_specification = \
        hrtem_system_model_params.core_attrs["sample_specification"]

    kwargs = \
        {"sample_specification": sample_specification}
    total_num_frozen_phonon_configs = \
        prismatique.sample._total_num_frozen_phonon_configs(**kwargs)
    
    return total_num_frozen_phonon_configs



def _serialize_sim_params(sim_params):
    output_params = sim_params.core_attrs["output_params"]
    output_dirname = output_params.core_attrs["output_dirname"]
    
    sim_params.dump(filename=output_dirname+"/hrtem_sim_params.json",
                    overwrite=True)

    return None
    


###########################
## Define error messages ##
###########################

_check_data_size_err_msg_1 = \
    ("The data size of the output of the HRTEM simulation to be performed is "
     "{} bytes, exceeding the maximum allowed size of {} bytes specified by "
     "the simulation parameters.")
