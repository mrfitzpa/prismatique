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
r"""Contains tests for the module :mod:`prismatique.load`.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For operations related to unit tests.
import pytest

# For saving objects to HDF5 files.
import h5pywrappers

# For type-checking objects.
import czekitout.isa

# For general array handling.
import numpy as np



# For simulating electron microscopy experiments using Multislice algorithms.
import prismatique



##################################
## Define classes and functions ##
##################################



def test_1_of_scan_pattern_type(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()
    
    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/stem_sim_params.json",
              "skip_validation_and_conversion": True}
    pattern_type = prismatique.load.scan_pattern_type(**kwargs)
    assert (pattern_type == "no underlying rectangular grid")

    kwargs = {"filename": output_dirname + "/stem_sim_intensity_output.h5",
              "skip_validation_and_conversion": False}
    pattern_type = prismatique.load.scan_pattern_type(**kwargs)
    assert (pattern_type == "no underlying rectangular grid")

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        prismatique.load.scan_pattern_type(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_grid_dims_in_units_of_probe_shifts(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()
    
    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/stem_sim_params.json",
              "skip_validation_and_conversion": True}
    grid_dims = prismatique.load.grid_dims_in_units_of_probe_shifts(**kwargs)
    assert (grid_dims == "N/A")

    kwargs = {"filename": output_dirname + "/stem_sim_intensity_output.h5",
              "skip_validation_and_conversion": False}
    grid_dims = prismatique.load.grid_dims_in_units_of_probe_shifts(**kwargs)
    assert (grid_dims == "N/A")

    helpers.generate_atomic_coords_file_2()

    sim_params = helpers.generate_stem_sim_params_2()
    helpers.generate_dummy_sim_output_files(sim_params)

    kwargs = {"filename": output_dirname + "/stem_sim_intensity_output.h5",
              "skip_validation_and_conversion": False}
    grid_dims = prismatique.load.grid_dims_in_units_of_probe_shifts(**kwargs)
    assert (grid_dims == (1, 1))

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        prismatique.load.grid_dims_in_units_of_probe_shifts(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_probe_positions(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()
    
    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/stem_sim_params.json",
              "force_2_col_shape": False,
              "skip_validation_and_conversion": True}
    probe_positions = np.array(prismatique.load.probe_positions(**kwargs))
    assert czekitout.isa.real_two_column_numpy_matrix(probe_positions)

    helpers.generate_atomic_coords_file_2()

    sim_params = helpers.generate_stem_sim_params_2()
    helpers.generate_dummy_sim_output_files(sim_params)

    for force_2_col_shape in (False, True):
        kwargs = {"filename": output_dirname + "/stem_sim_intensity_output.h5",
                  "force_2_col_shape": force_2_col_shape,
                  "skip_validation_and_conversion": False}
        probe_positions = np.array(prismatique.load.probe_positions(**kwargs))
        assert czekitout.isa.real_numpy_array(probe_positions)
        assert (len(probe_positions.shape) == 3-force_2_col_shape)

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        probe_positions = prismatique.load.probe_positions(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_output_layer_depths(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/stem_sim_params.json",
              "skip_validation_and_conversion": True}
    layer_depths = np.array(prismatique.load.output_layer_depths(**kwargs))
    assert czekitout.isa.real_numpy_array_1d(layer_depths)

    kwargs = {"filename": output_dirname + "/stem_sim_intensity_output.h5",
              "skip_validation_and_conversion": False}
    layer_depths = np.array(prismatique.load.output_layer_depths(**kwargs))
    assert czekitout.isa.real_numpy_array_1d(layer_depths)

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        prismatique.load.output_layer_depths(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_defocii(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/stem_sim_params.json",
              "skip_validation_and_conversion": True}
    beam_defocii = np.array(prismatique.load.defocii(**kwargs))
    assert czekitout.isa.real_numpy_array_1d(beam_defocii)

    kwargs = {"filename": \
              output_dirname + "/stem_sim_wavefunction_output_of_subset_0.h5",
              "skip_validation_and_conversion": \
              False}
    beam_defocii = np.array(prismatique.load.defocii(**kwargs))
    assert czekitout.isa.real_numpy_array_1d(beam_defocii)

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        prismatique.load.defocii(**kwargs)

    sim_params = helpers.generate_hrtem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    kwargs = {"filename": output_dirname + "/hrtem_sim_params.json",
              "skip_validation_and_conversion": False}
    beam_defocii = np.array(prismatique.load.defocii(**kwargs))
    assert czekitout.isa.real_numpy_array_1d(beam_defocii)

    helpers.remove_output_files()

    return None



def test_1_of_num_slices(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/sample_model_params.json",
              "skip_validation_and_conversion": True}
    num_sample_supercell_slices = prismatique.load.num_slices(**kwargs)
    assert (isinstance(num_sample_supercell_slices, int))

    kwargs = {"filename": \
              output_dirname + "/potential_slices_of_subset_0.h5",
              "skip_validation_and_conversion": \
              False}
    num_sample_supercell_slices = prismatique.load.num_slices(**kwargs)
    assert (isinstance(num_sample_supercell_slices, int))

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        prismatique.load.num_slices(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_num_frozen_phonon_configs_in_subset(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()
    helpers.generate_invalid_sim_output_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    func_alias = prismatique.load.num_frozen_phonon_configs_in_subset

    filenames = \
        (output_dirname + "/potential_slices_of_subset_0.h5",
         output_dirname + "/S_matrices_of_subset_0.h5",
         output_dirname + "/stem_sim_wavefunction_output_of_subset_0.h5")

    for filename in filenames:
        kwargs = {"filename": filename,
                  "skip_validation_and_conversion": True}
        num_configs_in_subset = func_alias(**kwargs)
        assert (isinstance(num_configs_in_subset, int))

    sim_params = helpers.generate_hrtem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    kwargs = {"filename": \
              output_dirname + "/hrtem_sim_wavefunction_output_of_subset_0.h5",
              "skip_validation_and_conversion": \
              False}
    num_configs_in_subset = func_alias(**kwargs)
    assert (isinstance(num_configs_in_subset, int))

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_invalid_sim_output_filename_1()
        func_alias(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_k_x_coords(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/stem_sim_intensity_output.h5",
              "skip_validation_and_conversion": False}
    k_x = np.array(prismatique.load.cbed_k_x_coords(**kwargs))
    assert czekitout.isa.real_numpy_array_1d(k_x)

    with pytest.raises(IOError) as err_info:
        kwargs = {"filename": helpers.generate_atomic_coords_filename_1(),
                  "skip_validation_and_conversion": True}
        prismatique.load.cbed_k_x_coords(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_k_y_coords(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/stem_sim_intensity_output.h5",
              "skip_validation_and_conversion": False}
    k_y = np.array(prismatique.load.cbed_k_y_coords(**kwargs))
    assert czekitout.isa.real_numpy_array_1d(k_y)

    with pytest.raises(IOError) as err_info:
        kwargs = {"filename": helpers.generate_atomic_coords_filename_1(),
                  "skip_validation_and_conversion": True}
        prismatique.load.cbed_k_y_coords(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_k_xy_coords_of_3d_stem_output(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/stem_sim_intensity_output.h5",
              "skip_validation_and_conversion": False}
    k_xy = np.array(prismatique.load.k_xy_coords_of_3d_stem_output(**kwargs))
    assert czekitout.isa.real_numpy_array_1d(k_xy)

    with pytest.raises(IOError) as err_info:
        kwargs = {"filename": helpers.generate_atomic_coords_filename_1(),
                  "skip_validation_and_conversion": True}
        prismatique.load.k_xy_coords_of_3d_stem_output(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_integration_limits_of_2d_stem_output(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    func_alias = prismatique.load.integration_limits_of_2d_stem_output

    kwargs = {"filename": output_dirname + "/stem_sim_intensity_output.h5",
              "skip_validation_and_conversion": False}
    integration_limits = func_alias(**kwargs)

    kwargs = {"obj": integration_limits, "obj_name": "integration_limits"}
    czekitout.check.if_pair_of_nonnegative_floats(**kwargs)

    with pytest.raises(IOError) as err_info:
        kwargs = {"filename": helpers.generate_atomic_coords_filename_1(),
                  "skip_validation_and_conversion": True}
        func_alias(**kwargs)

    attr_names = ("lower integration limit in mrads",
                  "upper integration limit in mrads")

    kwargs = {"filename": output_dirname + "/stem_sim_intensity_output.h5",
              "path_in_file": "/data/2D_STEM/integrated_intensities"}
    obj_id = h5pywrappers.obj.ID(**kwargs)

    for attr_name_idx, attr_name in enumerate(attr_names):
        kwargs = {"attr": -attr_name_idx+1,
                  "attr_id": h5pywrappers.attr.ID(obj_id, attr_name),
                  "write_mode": "a"}
        h5pywrappers.attr.save(**kwargs)

    with pytest.raises(IOError) as err_info:
        kwargs = {"filename": output_dirname + "/stem_sim_intensity_output.h5",
                  "skip_validation_and_conversion": False}
        func_alias(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_S_matrix_k_xy_vectors(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()
    
    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/stem_sim_params.json",
              "skip_validation_and_conversion": True}
    k_xy_vectors = np.array(prismatique.load.S_matrix_k_xy_vectors(**kwargs))
    assert czekitout.isa.real_two_column_numpy_matrix(k_xy_vectors)

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        prismatique.load.S_matrix_k_xy_vectors(**kwargs)

    helpers.generate_atomic_coords_file_2()

    sim_params = helpers.generate_stem_sim_params_4()
    helpers.generate_dummy_sim_output_files(sim_params)

    with pytest.raises(IOError) as err_info:
        kwargs = {"filename": output_dirname + "/stem_sim_params.json",
                  "skip_validation_and_conversion": False}
        prismatique.load.S_matrix_k_xy_vectors(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_hrtem_beam_tilts(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_hrtem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/hrtem_sim_params.json",
              "skip_validation_and_conversion": True}
    beam_tilts = np.array(prismatique.load.hrtem_beam_tilts(**kwargs))
    assert czekitout.isa.real_two_column_numpy_matrix(beam_tilts)

    kwargs = {"filename": \
              output_dirname + "/hrtem_sim_wavefunction_output_of_subset_0.h5",
              "skip_validation_and_conversion": \
              False}
    beam_tilts = np.array(prismatique.load.hrtem_beam_tilts(**kwargs))
    assert czekitout.isa.real_two_column_numpy_matrix(beam_tilts)

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        prismatique.load.hrtem_beam_tilts(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_hrtem_image_x_coords(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_hrtem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/hrtem_sim_intensity_output.h5",
              "skip_validation_and_conversion": False}
    r_x = np.array(prismatique.load.hrtem_image_x_coords(**kwargs))
    assert czekitout.isa.real_numpy_array_1d(r_x)

    with pytest.raises(IOError) as err_info:
        kwargs = {"filename": helpers.generate_atomic_coords_filename_1(),
                  "skip_validation_and_conversion": True}
        prismatique.load.hrtem_image_x_coords(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_hrtem_image_y_coords(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_hrtem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/hrtem_sim_intensity_output.h5",
              "skip_validation_and_conversion": False}
    r_y = np.array(prismatique.load.hrtem_image_y_coords(**kwargs))
    assert czekitout.isa.real_numpy_array_1d(r_y)

    with pytest.raises(IOError) as err_info:
        kwargs = {"filename": helpers.generate_atomic_coords_filename_1(),
                  "skip_validation_and_conversion": True}
        prismatique.load.hrtem_image_y_coords(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_potential_slices(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    multi_dim_slices = (None, (0, 0), (slice(None), 0), (0, slice(None)))

    for multi_dim_slice in multi_dim_slices:
        for num_superslices in (0, 1):
            for average_thermally in (True, False):
                for average_laterally_in_space in (True, False):
                    filename = (output_dirname
                                + "/potential_slices_of_subset_0.h5")
                    kwargs = {"filename": \
                              filename,
                              "multi_dim_slice": \
                              multi_dim_slice,
                              "num_superslices": \
                              num_superslices,
                              "average_thermally": \
                              average_thermally,
                              "average_laterally_in_space": \
                              average_laterally_in_space,
                              "skip_validation_and_conversion": \
                              False}
                    result = prismatique.load.potential_slices(**kwargs)
                    potential_slice_signal = result[0]
                    navigational_to_original_indices_map = result[1]

    helpers.generate_atomic_coords_file_2()

    sim_params = helpers.generate_stem_sim_params_2()
    helpers.generate_dummy_sim_output_files(sim_params)

    result = prismatique.load.potential_slices(**kwargs)

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        kwargs["skip_validation_and_conversion"] = True
        prismatique.load.potential_slices(**kwargs)

    helpers.remove_output_files()

    return None



def test_2_of_potential_slices(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_2()

    sim_params = helpers.generate_stem_sim_params_2()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/potential_slices_of_subset_0.h5",
              "multi_dim_slice": None,
              "num_superslices": 2,
              "average_thermally": False,
              "average_laterally_in_space": False,
              "skip_validation_and_conversion": False}
    prismatique.load.potential_slices(**kwargs)

    error_types = (IndexError, ValueError, ValueError)
    multi_dim_slices = ((slice(None), 5), (slice(None), [0, 0]), (0, 0, 0))
    zip_obj = zip(error_types, multi_dim_slices)

    for error_type, multi_dim_slice in zip_obj:
        with pytest.raises(error_type) as err_info:
            kwargs["multi_dim_slice"] = multi_dim_slice
            prismatique.load.potential_slices(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_S_matrix_wavefunctions(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_2()

    sim_params = helpers.generate_stem_sim_params_2()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    multi_dim_slices = (None, (0, 0), (slice(None), 0), (0, slice(None)))

    for multi_dim_slice in multi_dim_slices:
        kwargs = {"filename": output_dirname + "/S_matrices_of_subset_0.h5",
                  "multi_dim_slice": multi_dim_slice,
                  "skip_validation_and_conversion": False}
        result = prismatique.load.S_matrix_wavefunctions(**kwargs)
        S_matrix_wavefunction_signal = result[0]
        navigational_to_original_indices_map = result[1]

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        kwargs["skip_validation_and_conversion"] = True
        prismatique.load.S_matrix_wavefunctions(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_cbed_wavefunctions(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    multi_dim_slices = ((0, 0, 0, 0), None)

    for multi_dim_slice in multi_dim_slices:
        for use_two_axes_to_map_probe_position_if_possible in (False, True):
            filename = (output_dirname
                        + "/stem_sim_wavefunction_output_of_subset_0.h5")
            kwargs = {"filename": \
                      filename,
                      "multi_dim_slice": \
                      multi_dim_slice,
                      "use_two_axes_to_map_probe_position_if_possible": \
                      use_two_axes_to_map_probe_position_if_possible,
                      "skip_validation_and_conversion": \
                      False}
            result = prismatique.load.cbed_wavefunctions(**kwargs)
            cbed_wavefunction_signal = result[0]
            navigational_to_original_indices_map = result[1]

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        kwargs["skip_validation_and_conversion"] = True
        prismatique.load.cbed_wavefunctions(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_cbed_intensity_patterns(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_2()

    sim_params = helpers.generate_stem_sim_params_4()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    multi_dim_slices = ((0, 0), None)

    for multi_dim_slice in multi_dim_slices:
        for use_two_axes_to_map_probe_position_if_possible in (False, True):
            kwargs = {"filename": \
                      output_dirname + "/stem_sim_intensity_output.h5",
                      "multi_dim_slice": \
                      multi_dim_slice,
                      "use_two_axes_to_map_probe_position_if_possible": \
                      use_two_axes_to_map_probe_position_if_possible,
                      "skip_validation_and_conversion": \
                      False}
            result = prismatique.load.cbed_intensity_patterns(**kwargs)
            cbed_intensity_pattern_signal = result[0]
            navigational_to_original_indices_map = result[1]

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        kwargs["skip_validation_and_conversion"] = True
        prismatique.load.cbed_intensity_patterns(**kwargs)

    kwargs = {"filename": output_dirname + "/stem_sim_intensity_output.h5",
              "path_in_file": "/data/4D_STEM/intensity_DPs"}
    dataset_id = h5pywrappers.obj.ID(**kwargs)

    kwargs = {"dataset": np.zeros((3, 3)),
              "dataset_id": dataset_id,
              "write_mode": "a"}
    h5pywrappers.dataset.save(**kwargs)

    with pytest.raises(IOError) as err_info:
        kwargs = {"filename": output_dirname + "/stem_sim_intensity_output.h5",
                  "multi_dim_slice": None,
                  "use_two_axes_to_map_probe_position_if_possible": True}
        prismatique.load.cbed_intensity_patterns(**kwargs)
        
    helpers.remove_output_files()

    return None



def test_1_of_com_momenta(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_2()

    sim_params = helpers.generate_stem_sim_params_4()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    multi_dim_slices = ((0, 0), None)

    for multi_dim_slice in multi_dim_slices:
        for use_two_axes_to_map_probe_position_if_possible in (False, True):
            filename = output_dirname + "/stem_sim_intensity_output.h5"
            kwargs = {"filename": \
                      filename,
                      "multi_dim_slice": \
                      multi_dim_slice,
                      "use_two_axes_to_map_probe_position_if_possible": \
                      use_two_axes_to_map_probe_position_if_possible,
                      "skip_validation_and_conversion": \
                      False}
            result = prismatique.load.com_momenta(**kwargs)
            com_momentum_signal = result[0]
            navigational_to_original_indices_map = result[1]

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        kwargs["skip_validation_and_conversion"] = True
        prismatique.load.com_momenta(**kwargs)
        
    helpers.remove_output_files()

    return None



def test_1_of_stem_intensity_images(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_2()

    sim_params = helpers.generate_stem_sim_params_4()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    multi_dim_slices = ((0,), None)

    for multi_dim_slice in multi_dim_slices:
        for use_two_axes_to_map_probe_position_if_possible in (False, True):
            filename = output_dirname + "/stem_sim_intensity_output.h5"
            kwargs = {"filename": \
                      filename,
                      "multi_dim_slice": \
                      multi_dim_slice,
                      "use_two_axes_to_map_probe_position_if_possible": \
                      use_two_axes_to_map_probe_position_if_possible,
                      "skip_validation_and_conversion": \
                      False}
            result = prismatique.load.stem_intensity_images(**kwargs)
            stem_image_signal = result[0]
            navigational_to_original_indices_map = result[1]

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        kwargs["skip_validation_and_conversion"] = True
        prismatique.load.stem_intensity_images(**kwargs)
        
    helpers.remove_output_files()

    return None



def test_1_of_azimuthally_integrated_cbed_intensity_patterns(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_2()

    sim_params = helpers.generate_stem_sim_params_4()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    multi_dim_slices = ((0, 0), None)

    func_alias = prismatique.load.azimuthally_integrated_cbed_intensity_patterns

    for multi_dim_slice in multi_dim_slices:
        for use_two_axes_to_map_probe_position_if_possible in (False, True):
            filename = output_dirname + "/stem_sim_intensity_output.h5"
            kwargs = {"filename": \
                      filename,
                      "multi_dim_slice": \
                      multi_dim_slice,
                      "use_two_axes_to_map_probe_position_if_possible": \
                      use_two_axes_to_map_probe_position_if_possible,
                      "skip_validation_and_conversion": \
                      False}
            result = func_alias(**kwargs)
            stem_image_signal = result[0]
            navigational_to_original_indices_map = result[1]

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        kwargs["skip_validation_and_conversion"] = True
        func_alias(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_hrtem_image_wavefunctions(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_hrtem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    multi_dim_slices = ((0, 0, 0), None)

    for multi_dim_slice in multi_dim_slices:
        filename = (output_dirname
                    + "/hrtem_sim_wavefunction_output_of_subset_0.h5")
        kwargs = {"filename": filename,
                  "multi_dim_slice": multi_dim_slice,
                  "skip_validation_and_conversion": False}
        result = prismatique.load.hrtem_image_wavefunctions(**kwargs)
        hrtem_image_wavefunction_signal = result[0]
        navigational_to_original_indices_map = result[1]

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        kwargs["skip_validation_and_conversion"] = True
        prismatique.load.hrtem_image_wavefunctions(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_hrtem_intensity_image(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_2()

    sim_params = helpers.generate_hrtem_sim_params_4()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    kwargs = {"filename": output_dirname + "/hrtem_sim_intensity_output.h5",
              "skip_validation_and_conversion": False}
    result = prismatique.load.hrtem_intensity_image(**kwargs)
    hrtem_intensity_image_signal = result

    with pytest.raises(IOError) as err_info:
        kwargs["filename"] = helpers.generate_atomic_coords_filename_1()
        kwargs["skip_validation_and_conversion"] = True
        prismatique.load.hrtem_intensity_image(**kwargs)

    kwargs = {"filename": output_dirname + "/hrtem_sim_intensity_output.h5",
              "path_in_file": "/data/intensity_image"}
    dataset_id = h5pywrappers.obj.ID(**kwargs)

    kwargs = {"dataset": np.zeros((3, 3, 3)),
              "dataset_id": dataset_id,
              "write_mode": "a"}
    h5pywrappers.dataset.save(**kwargs)

    with pytest.raises(IOError) as err_info:
        kwargs = {"filename": output_dirname + "/hrtem_sim_intensity_output.h5"}
        prismatique.load.hrtem_intensity_image(**kwargs)
        
    helpers.remove_output_files()

    return None



###########################
## Define error messages ##
###########################
