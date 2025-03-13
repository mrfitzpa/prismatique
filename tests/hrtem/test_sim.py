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
r"""Contains tests for the module :mod:`prismatique.hrtem.sim`.

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
    sim_params = helpers.generate_hrtem_sim_params_1()

    sim_params.validation_and_conversion_funcs
    sim_params.pre_serialization_funcs
    sim_params.de_pre_serialization_funcs

    kwargs = {"serializable_rep": sim_params.pre_serialize()}
    prismatique.hrtem.sim.Params.de_pre_serialize(**kwargs)

    return None



def test_1_of_run(helpers):
    helpers.remove_output_files()
    
    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_hrtem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)
    prismatique.hrtem.sim.run(sim_params)

    helpers.generate_atomic_coords_file_2()

    sim_params = helpers.generate_hrtem_sim_params_2()
    prismatique.hrtem.sim.run(sim_params)
    prismatique.hrtem.sim.run(sim_params)

    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_hrtem_sim_params_3()
    with pytest.raises(MemoryError) as err_info:
        prismatique.hrtem.sim.run(sim_params)

    sim_params = helpers.generate_hrtem_sim_params_4()
    prismatique.hrtem.sim.run(sim_params)

    helpers.remove_output_files()

    return None



###########################
## Define error messages ##
###########################
