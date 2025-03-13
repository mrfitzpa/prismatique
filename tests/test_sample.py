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
r"""Contains tests for the module :mod:`prismatique.sample`.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For moving files.
import shutil



# For operations related to unit tests.
import pytest

# For saving objects to HDF5 files.
import h5pywrappers

# For deleting objects stored in HDF5 files.
import h5py

# For general array handling.
import numpy as np



# For simulating electron microscopy experiments using Multislice algorithms.
import prismatique



##################################
## Define classes and functions ##
##################################



def test_1_of_ModelParams():
    with pytest.raises(TypeError) as err_info:
        kwargs = {"atomic_coords_filename": "foobar.txt",
                  "unit_cell_tiling": (1, 1)}
        prismatique.sample.ModelParams(**kwargs)

    return None



def test_1_of_check_atomic_coords_file_format(helpers):
    helpers.remove_output_files()
    
    helpers.generate_invalid_atomic_coords_file_1()
    helpers.generate_invalid_atomic_coords_file_2()
    helpers.generate_invalid_atomic_coords_file_3()

    filenames = (helpers.generate_atomic_coords_filename_3(),
                 helpers.generate_atomic_coords_filename_4(),
                 helpers.generate_atomic_coords_filename_5(),
                 "foobar")
    type_errors = 3*(IOError,) + (FileNotFoundError,)
    zip_obj = zip(filenames, type_errors)

    for filename, type_error in zip_obj:
        with pytest.raises(type_error) as err_info:
            kwargs = {"atomic_coords_filename": filename}
            prismatique.sample.check_atomic_coords_file_format(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_PotentialSliceSubsetIDs():
    kwargs = \
        {"filenames": ("foobar",)}
    potential_slice_subset_ids = \
        prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"serializable_rep": potential_slice_subset_ids.pre_serialize()}
    prismatique.sample.PotentialSliceSubsetIDs.de_pre_serialize(**kwargs)

    with pytest.raises(ValueError) as err_info:
        kwargs = {"filenames": tuple()}
        prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    with pytest.raises(ValueError) as err_info:
        kwargs = {"filenames": ("foobar",),
                  "max_num_frozen_phonon_configs_per_subset": 1}
        prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    return None



def test_1_of_SMatrixSubsetIDs():
    kwargs = \
        {"filenames": ("foobar",)}
    potential_slice_subset_ids = \
        prismatique.sample.SMatrixSubsetIDs(**kwargs)

    kwargs = {"serializable_rep": potential_slice_subset_ids.pre_serialize()}
    prismatique.sample.SMatrixSubsetIDs.de_pre_serialize(**kwargs)

    with pytest.raises(ValueError) as err_info:
        kwargs = {"filenames": tuple()}
        prismatique.sample.SMatrixSubsetIDs(**kwargs)

    return None



def test_1_of_PotentialSliceAndSMatrixSubsetIDs():
    kwargs = \
        {"potential_slice_subset_filenames": ("foobar",),
         "S_matrix_subset_filenames": ("foobar",)}
    sample_specification = \
        prismatique.sample.PotentialSliceAndSMatrixSubsetIDs(**kwargs)

    method_alias = \
        prismatique.sample.PotentialSliceAndSMatrixSubsetIDs.de_pre_serialize
    kwargs = \
        {"serializable_rep": sample_specification.pre_serialize()}
    _ = \
        method_alias(**kwargs)

    with pytest.raises(ValueError) as err_info:
        kwargs = {"potential_slice_subset_filenames": tuple(),
                  "S_matrix_subset_filenames": ("foobar",)}
        prismatique.sample.PotentialSliceAndSMatrixSubsetIDs(**kwargs)

    with pytest.raises(ValueError) as err_info:
        kwargs = {"potential_slice_subset_filenames": ("foobar",),
                  "S_matrix_subset_filenames": tuple()}
        prismatique.sample.PotentialSliceAndSMatrixSubsetIDs(**kwargs)

    return None



def test_1_of_unit_cell_dims(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_hrtem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "skip_validation_and_conversion": True}
    prismatique.sample.unit_cell_dims(**kwargs)

    path_in_file = ("/4DSTEM_simulation/metadata/metadata_0/original"
                    "/simulation_parameters")
    
    kwargs = {"filename": filename, "path_in_file": path_in_file}
    obj_id = h5pywrappers.obj.ID(**kwargs)

    hdf5_attr_names = ("c", "t")
    new_hdf5_attrs = ((5.0, 5.0), (1, 1))
    zip_obj = zip(hdf5_attr_names, new_hdf5_attrs)

    for hdf5_attr_name, new_hdf5_attr in zip_obj:
        kwargs = {"obj_id": obj_id, "attr_name": hdf5_attr_name}
        attr_id = h5pywrappers.attr.ID(**kwargs)
        
        old_hdf5_attr = h5pywrappers.attr.load(attr_id)

        kwargs = {"attr": new_hdf5_attr, "attr_id": attr_id, "write_mode": "a"}
        h5pywrappers.attr.save(**kwargs)

        with pytest.raises(TypeError) as err_info:
            kwargs = {"sample_specification": sample_specification,
                      "skip_validation_and_conversion": False}
            prismatique.sample.unit_cell_dims(**kwargs)

        kwargs = {"attr": old_hdf5_attr, "attr_id": attr_id, "write_mode": "a"}
        h5pywrappers.attr.save(**kwargs)
        
    helpers.remove_output_files()

    return None



def test_2_of_unit_cell_dims(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.SMatrixSubsetIDs(**kwargs)

    with pytest.raises(IOError) as err_info:
        kwargs = {"sample_specification": sample_specification,
                  "skip_validation_and_conversion": False}
        prismatique.sample.unit_cell_dims(**kwargs)

    path_in_file = "/4DSTEM_simulation/data/realslices/ppotential_fp0000/data"

    kwargs = {"filename": filename, "path_in_file": path_in_file}
    dataset_id = h5pywrappers.obj.ID(**kwargs)

    datasets = (np.zeros((3, 3, 3), dtype=np.complex64), np.zeros((3, 3)))
    for dataset in datasets:
        kwargs = {"dataset": dataset,
                  "dataset_id": dataset_id,
                  "write_mode": "a"}
        h5pywrappers.dataset.save(**kwargs)

        cls_alias = prismatique.sample.PotentialSliceSubsetIDs
        kwargs = {"filenames": (filename,)}
        sample_specification = cls_alias(**kwargs)

        with pytest.raises(IOError) as err_info:
            kwargs = {"sample_specification": sample_specification,
                      "skip_validation_and_conversion": False}
            prismatique.sample.unit_cell_dims(**kwargs)
        
    helpers.remove_output_files()

    return None



def test_3_of_unit_cell_dims(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    filename = output_dirname + "/S_matrices_of_subset_0.h5"
    path_in_file = "/4DSTEM_simulation/data/realslices/smatrix_fp0000/data"

    kwargs = {"filename": filename, "path_in_file": path_in_file}
    dataset_id = h5pywrappers.obj.ID(**kwargs)

    datasets = (np.zeros((3, 3, 3)), np.zeros((3, 3)))
    for dataset in datasets:
        kwargs = {"dataset": dataset,
                  "dataset_id": dataset_id,
                  "write_mode": "a"}
        h5pywrappers.dataset.save(**kwargs)

        kwargs = {"filenames": (filename,)}
        sample_specification = prismatique.sample.SMatrixSubsetIDs(**kwargs)

        with pytest.raises(IOError) as err_info:
            kwargs = {"sample_specification": sample_specification,
                      "skip_validation_and_conversion": False}
            prismatique.sample.unit_cell_dims(**kwargs)

    with h5py.File(filename, "a") as file_obj:
        path_in_file = "/4DSTEM_simulation/data/realslices"
        del file_obj[path_in_file]

    with pytest.raises(IOError) as err_info:
        kwargs = {"sample_specification": sample_specification,
                  "skip_validation_and_conversion": False}
        prismatique.sample.unit_cell_dims(**kwargs)
        
    helpers.remove_output_files()

    return None



def test_4_of_unit_cell_dims(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()
    helpers.generate_atomic_coords_file_2()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname_1 = helpers.generate_output_dirname_1()
    output_dirname_2 = helpers.generate_output_dirname_2()

    shutil.move(output_dirname_1, output_dirname_2)

    sim_params = helpers.generate_stem_sim_params_2()
    helpers.generate_dummy_sim_output_files(sim_params)

    filename_1 = output_dirname_1 + "/potential_slices_of_subset_0.h5"
    filename_2 = output_dirname_2 + "/potential_slices_of_subset_0.h5"
    filename_3 = output_dirname_2 + "/S_matrices_of_subset_0.h5"

    cls_alias_set = (prismatique.sample.PotentialSliceSubsetIDs,
                     prismatique.sample.PotentialSliceAndSMatrixSubsetIDs)
    
    for cls_alias in cls_alias_set:
        if cls_alias == prismatique.sample.PotentialSliceSubsetIDs:
            kwargs = {"filenames": (filename_1, filename_2)}
        else:
            kwargs = {"potential_slice_subset_filenames": (filename_1,),
                      "S_matrix_subset_filenames": (filename_3,)}
        sample_specification = cls_alias(**kwargs)

        with pytest.raises(ValueError) as err_info:
            kwargs = {"sample_specification": sample_specification,
                      "skip_validation_and_conversion": False}
            prismatique.sample.unit_cell_dims(**kwargs)

    shutil.move(output_dirname_2, output_dirname_1)

    for skip_validation_and_conversion in (False, True):
        kwargs = {"sample_specification": \
                  helpers.generate_sample_specification_1(),
                  "skip_validation_and_conversion": \
                  skip_validation_and_conversion}
        prismatique.sample.unit_cell_dims(**kwargs)
    
    helpers.remove_output_files()

    return None



def test_1_of_supercell_dims(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "skip_validation_and_conversion": False}
    prismatique.sample.supercell_dims(**kwargs)

    kwargs = {"sample_specification": helpers.generate_sample_specification_1(),
              "skip_validation_and_conversion": True}
    prismatique.sample.supercell_dims(**kwargs)
    
    helpers.remove_output_files()

    return None



def test_1_of_supercell_xy_dims_in_pixels(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "skip_validation_and_conversion": False}
    prismatique.sample.supercell_xy_dims_in_pixels(**kwargs)

    kwargs = {"sample_specification": helpers.generate_sample_specification_1(),
              "skip_validation_and_conversion": True}
    prismatique.sample.supercell_xy_dims_in_pixels(**kwargs)
    
    helpers.remove_output_files()

    return None



def test_1_of_supercell_lateral_pixel_size(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "skip_validation_and_conversion": False}
    prismatique.sample.supercell_lateral_pixel_size(**kwargs)

    kwargs = {"sample_specification": helpers.generate_sample_specification_1(),
              "skip_validation_and_conversion": True}
    prismatique.sample.supercell_lateral_pixel_size(**kwargs)
    
    helpers.remove_output_files()

    return None



def test_1_of_supercell_slice_thickness(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "skip_validation_and_conversion": False}
    prismatique.sample.supercell_slice_thickness(**kwargs)

    kwargs = {"sample_specification": helpers.generate_sample_specification_1(),
              "skip_validation_and_conversion": True}
    prismatique.sample.supercell_slice_thickness(**kwargs)
    
    helpers.remove_output_files()

    return None



def test_1_of_num_slices(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "skip_validation_and_conversion": False}
    prismatique.sample.num_slices(**kwargs)

    kwargs = {"sample_specification": helpers.generate_sample_specification_1(),
              "skip_validation_and_conversion": True}
    prismatique.sample.num_slices(**kwargs)
    
    helpers.remove_output_files()

    return None



def test_1_of_num_frozen_phonon_config_subsets(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "skip_validation_and_conversion": False}
    prismatique.sample.num_frozen_phonon_config_subsets(**kwargs)

    kwargs = {"sample_specification": helpers.generate_sample_specification_1(),
              "skip_validation_and_conversion": True}
    prismatique.sample.num_frozen_phonon_config_subsets(**kwargs)
    
    helpers.remove_output_files()

    return None



def test_1_of_num_frozen_phonon_configs_in_subset(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename_1 = output_dirname + "/potential_slices_of_subset_0.h5"
    filename_2 = output_dirname + "/S_matrices_of_subset_0.h5"

    kwargs = {"filenames": (filename_1,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "skip_validation_and_conversion": False}
    prismatique.sample.num_frozen_phonon_configs_in_subset(**kwargs)

    kwargs = \
        {"potential_slice_subset_filenames": (filename_1,),
         "S_matrix_subset_filenames": (filename_2,)}
    sample_specification = \
        prismatique.sample.PotentialSliceAndSMatrixSubsetIDs(**kwargs)

    for subset_idx in (0, 1):
        kwargs = {"sample_specification": sample_specification,
                  "subset_idx": subset_idx,
                  "skip_validation_and_conversion": False}
        prismatique.sample.num_frozen_phonon_configs_in_subset(**kwargs)

    with pytest.raises(ValueError) as err_info:
        kwargs = {"sample_specification": sample_specification,
                  "subset_idx": 10,
                  "skip_validation_and_conversion": False}
        prismatique.sample.num_frozen_phonon_configs_in_subset(**kwargs)

    kwargs = {"sample_specification": helpers.generate_sample_specification_1(),
              "skip_validation_and_conversion": True}
    prismatique.sample.num_frozen_phonon_configs_in_subset(**kwargs)
    
    helpers.remove_output_files()

    return None



def test_1_of_total_num_frozen_phonon_configs(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "skip_validation_and_conversion": False}
    prismatique.sample.total_num_frozen_phonon_configs(**kwargs)

    kwargs = {"sample_specification": helpers.generate_sample_specification_1(),
              "skip_validation_and_conversion": True}
    prismatique.sample.total_num_frozen_phonon_configs(**kwargs)
    
    helpers.remove_output_files()

    return None



def test_1_of_S_matrix_k_xy_vectors(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename_1 = output_dirname + "/potential_slices_of_subset_0.h5"
    filename_2 = output_dirname + "/S_matrices_of_subset_0.h5"

    kwargs = {"filenames": (filename_1,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    for skip_validation_and_conversion in (False, True):
        kwargs = {"sample_specification": \
                  sample_specification,
                  "probe_model_params": \
                  helpers.generate_probe_model_params_2(),
                  "skip_validation_and_conversion": \
                  skip_validation_and_conversion}
        prismatique.sample.S_matrix_k_xy_vectors(**kwargs)

    kwargs = {"filenames": (filename_2,)}
    sample_specification = prismatique.sample.SMatrixSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "probe_model_params": helpers.generate_probe_model_params_2(),
              "skip_validation_and_conversion": False}
    prismatique.sample.S_matrix_k_xy_vectors(**kwargs)

    with pytest.raises(IOError) as err_info:
        kwargs = {"sample_specification": sample_specification,
                  "probe_model_params": helpers.generate_probe_model_params_1(),
                  "skip_validation_and_conversion": False}
        prismatique.sample.S_matrix_k_xy_vectors(**kwargs)
    
    helpers.remove_output_files()

    return None



def test_1_of_potential_slice_subset_data_size(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "skip_validation_and_conversion": False}
    prismatique.sample.potential_slice_subset_data_size(**kwargs)

    kwargs = {"sample_specification": helpers.generate_sample_specification_1(),
              "skip_validation_and_conversion": True}
    prismatique.sample.potential_slice_subset_data_size(**kwargs)
    
    helpers.remove_output_files()

    return None



def test_1_of_potential_slice_set_data_size(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "skip_validation_and_conversion": False}
    prismatique.sample.potential_slice_set_data_size(**kwargs)

    kwargs = {"sample_specification": helpers.generate_sample_specification_1(),
              "skip_validation_and_conversion": True}
    prismatique.sample.potential_slice_set_data_size(**kwargs)
    
    helpers.remove_output_files()

    return None



def test_1_of_S_matrix_subset_data_size(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "probe_model_params": helpers.generate_probe_model_params_1(),
              "skip_validation_and_conversion": False}
    prismatique.sample.S_matrix_subset_data_size(**kwargs)

    kwargs = {"sample_specification": helpers.generate_sample_specification_1(),
              "probe_model_params": helpers.generate_probe_model_params_1(),
              "skip_validation_and_conversion": True}
    prismatique.sample.S_matrix_subset_data_size(**kwargs)
    
    helpers.remove_output_files()

    return None



def test_1_of_S_matrix_set_data_size(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "probe_model_params": helpers.generate_probe_model_params_1(),
              "skip_validation_and_conversion": False}
    prismatique.sample.S_matrix_set_data_size(**kwargs)

    kwargs = {"sample_specification": helpers.generate_sample_specification_1(),
              "probe_model_params": helpers.generate_probe_model_params_1(),
              "skip_validation_and_conversion": True}
    prismatique.sample.S_matrix_set_data_size(**kwargs)
    
    helpers.remove_output_files()

    return None



def test_1_of_generate_potential_slices(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"sample_model_params": helpers.generate_sample_specification_1(),
              "output_dirname": output_dirname,
              "max_data_size": 2*10**9,
              "worker_params": helpers.generate_worker_params_1(),
              "skip_validation_and_conversion": False}
    prismatique.sample.generate_potential_slices(**kwargs)

    kwargs["skip_validation_and_conversion"] = True
    prismatique.sample.generate_potential_slices(**kwargs)

    with pytest.raises(MemoryError) as err_info:
        kwargs["max_data_size"] = 1
        kwargs["skip_validation_and_conversion"] = False
        prismatique.sample.generate_potential_slices(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_generate_S_matrices(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "probe_model_params": helpers.generate_probe_model_params_1(),
              "output_dirname": helpers.generate_output_dirname_1(),
              "max_data_size": 2*10**9,
              "worker_params": None,
              "skip_validation_and_conversion": False}
    prismatique.sample.generate_S_matrices(**kwargs)

    kwargs["skip_validation_and_conversion"] = True
    kwargs["worker_params"] = helpers.generate_worker_params_1()
    prismatique.sample.generate_S_matrices(**kwargs)

    helpers.remove_output_files()

    return None



###########################
## Define error messages ##
###########################
