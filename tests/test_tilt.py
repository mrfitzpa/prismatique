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
r"""Contains tests for the module :mod:`prismatique.tilt`.

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



def test_1_of_Params():
    windows = ((0, 0.5, 1), (1, 0))
    
    for window in windows:
        with pytest.raises(TypeError) as err_info:
            kwargs = {"window": window}
            prismatique.tilt.Params(**kwargs)

    return None



def test_1_of_step_size(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_hrtem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    sample_specification = helpers.generate_sample_specification_1()

    kwargs = {"sample_specification": sample_specification,
              "mean_beam_energy": 30,
              "skip_validation_and_conversion": True}
    prismatique.tilt.step_size(**kwargs)

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "mean_beam_energy": 30,
              "skip_validation_and_conversion": False}
    prismatique.tilt.step_size(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_series(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_hrtem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()
    filename = output_dirname + "/potential_slices_of_subset_0.h5"

    sample_specification = helpers.generate_sample_specification_1()

    kwargs = {"window": (0, 1, 0, 1)}
    tilt_params = prismatique.tilt.Params(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "mean_beam_energy": 30,
              "tilt_params": tilt_params,
              "skip_validation_and_conversion": True}
    prismatique.tilt.series(**kwargs)

    kwargs = {"filenames": (filename,)}
    sample_specification = prismatique.sample.PotentialSliceSubsetIDs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "mean_beam_energy": 30,
              "tilt_params": tilt_params,
              "skip_validation_and_conversion": False}
    prismatique.tilt.series(**kwargs)

    helpers.remove_output_files()

    return None



###########################
## Define error messages ##
###########################
