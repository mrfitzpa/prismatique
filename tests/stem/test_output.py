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
r"""Contains tests for the module :mod:`prismatique.stem.output`.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For operations related to unit tests.
import pytest



# For simulating electron microscopy experiments using Multislice algorithms.
import prismatique



##################################
## Define classes and functions ##
##################################



def test_1_of_Params(helpers):
    prismatique.stem.output.Params()

    with pytest.raises(TypeError) as err_info:
        kwargs = {"alg_specific_params": 0}
        prismatique.stem.output.Params(**kwargs)

    return None



def test_1_of_layer_depths(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename_1 = output_dirname + "/potential_slices_of_subset_0.h5"
    filename_2 = output_dirname + "/S_matrices_of_subset_0.h5"

    kwargs = {"filenames": (filename_1,), "interpolation_factors": (2, 1)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    with pytest.raises(ValueError) as err_info:
        kwargs = {"sample_specification": sample_specification}
        prismatique.stem.output.layer_depths(**kwargs)

    kwargs = {"filenames": (filename_1,), "interpolation_factors": (1, 1)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification}
    prismatique.stem.output.layer_depths(**kwargs)

    kwargs = {"filenames": (filename_2,)}
    sample_specification = prismatique.sample.SMatrixSubsetIDs(**kwargs)

    with pytest.raises(TypeError) as err_info:
        kwargs = {"sample_specification": sample_specification}
        prismatique.stem.output.layer_depths(**kwargs)

    helpers.remove_output_files()

    return None



def test_2_of_layer_depths(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sample_specification = helpers.generate_sample_specification_4()

    kwargs = {"num_slices_per_output": 1, "z_start_output": 0}
    alg_specific_params = prismatique.stem.output.multislice.Params(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "alg_specific_params": alg_specific_params,
              "skip_validation_and_conversion": True}
    prismatique.stem.output.layer_depths(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_data_size(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    stem_system_model_params = helpers.generate_stem_system_model_params_1()

    kwargs = {"stem_system_model_params": stem_system_model_params}
    prismatique.stem.output.data_size(**kwargs)

    for save_final_intensity in (True, False):
        kwargs = {"save_final_intensity": save_final_intensity}
        cbed_params = prismatique.cbed.Params(**kwargs)

        kwargs = {"cbed_params": cbed_params,
                  "radial_step_size_for_3d_stem": 0,
                  "radial_range_for_2d_stem": (0, 0),
                  "save_com": False}
        base_output_params = prismatique.stem.output.base.Params(**kwargs)

        kwargs = {"base_params": base_output_params}
        output_params = prismatique.stem.output.Params(**kwargs)

        kwargs = {"stem_system_model_params": stem_system_model_params,
                  "output_params": output_params,
                  "skip_validation_and_conversion": True}
        prismatique.stem.output.data_size(**kwargs)

    helpers.remove_output_files()

    return None



def test_2_of_data_size(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/S_matrices_of_subset_0.h5"

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.SMatrixSubsetIDs(**kwargs)

    stem_system_model_params = sim_params.core_attrs["stem_system_model_params"]
    
    new_core_attr_subset_candidate = {"sample_specification": \
                                      sample_specification}
    stem_system_model_params.update(new_core_attr_subset_candidate)

    bool_set = (False, True)
    error_types = (TypeError, ValueError)
    zip_obj = zip(bool_set, error_types)

    for save_potential_slices, error_type in zip_obj:
        kwargs = {"save_potential_slices": save_potential_slices}
        base_output_params = prismatique.stem.output.base.Params(**kwargs)

        kwargs = {"base_params": base_output_params}
        output_params = prismatique.stem.output.Params(**kwargs)

        with pytest.raises(error_type) as err_info:
            kwargs = {"stem_system_model_params": stem_system_model_params,
                      "output_params": output_params}
            prismatique.stem.output.data_size(**kwargs)

    helpers.remove_output_files()

    return None



def test_3_of_data_size(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    stem_system_model_params = sim_params.core_attrs["stem_system_model_params"]
            
    kwargs = {"filenames": (filename,), "interpolation_factors": (2, 1)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    new_core_attr_subset_candidate = {"sample_specification": \
                                      sample_specification}
    stem_system_model_params.update(new_core_attr_subset_candidate)

    with pytest.raises(ValueError) as err_info:
        kwargs = {"stem_system_model_params": stem_system_model_params}
        prismatique.stem.output.data_size(**kwargs)

    helpers.remove_output_files()

    return None



###########################
## Define error messages ##
###########################
