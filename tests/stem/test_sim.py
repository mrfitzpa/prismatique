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
r"""Contains tests for the module :mod:`prismatique.stem.sim`.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For moving files.
import shutil



# For operations related to unit tests.
import pytest



# For simulating electron microscopy experiments using Multislice algorithms.
import prismatique



##################################
## Define classes and functions ##
##################################



def test_1_of_Params(helpers):
    sim_params = helpers.generate_stem_sim_params_1()

    sim_params.validation_and_conversion_funcs
    sim_params.pre_serialization_funcs
    sim_params.de_pre_serialization_funcs

    kwargs = {"serializable_rep": sim_params.pre_serialize()}
    prismatique.stem.sim.Params.de_pre_serialize(**kwargs)

    sim_params = helpers.generate_stem_sim_params_3()

    kwargs = {"serializable_rep": sim_params.pre_serialize()}
    prismatique.stem.sim.Params.de_pre_serialize(**kwargs)

    return None



def test_1_of_run(helpers):
    helpers.remove_output_files()
    
    helpers.generate_atomic_coords_file_1()

    sim_params = helpers.generate_stem_sim_params_1()
    helpers.generate_dummy_sim_output_files(sim_params)
    prismatique.stem.sim.run(sim_params)

    output_dirname_1 = helpers.generate_output_dirname_1()
    output_dirname_2 = helpers.generate_output_dirname_2()

    shutil.move(output_dirname_1, output_dirname_2)

    helpers.generate_atomic_coords_file_2()

    sim_params = helpers.generate_stem_sim_params_2()
    prismatique.stem.sim.run(sim_params)
    prismatique.stem.sim.run(sim_params)

    helpers.generate_atomic_coords_file_1()
    helpers.generate_scan_config_file_1()

    sim_params = helpers.generate_stem_sim_params_3()
    with pytest.raises(MemoryError) as err_info:
        prismatique.stem.sim.run(sim_params)

    sim_params = helpers.generate_stem_sim_params_4()
    prismatique.stem.sim.run(sim_params)

    sim_params = helpers.generate_stem_sim_params_5()
    prismatique.stem.sim.run(sim_params)

    shutil.move(output_dirname_2, output_dirname_1)

    helpers.remove_output_files()

    return None



###########################
## Define error messages ##
###########################
