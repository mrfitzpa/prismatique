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
r"""Contains tests for the module :mod:`prismatique.stem.system`.

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



def test_1_of_ModelParams(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_2()
    helpers.generate_scan_config_file_1()

    sim_params = helpers.generate_stem_sim_params_3()
    helpers.generate_dummy_sim_output_files(sim_params)

    output_dirname = helpers.generate_output_dirname_1()

    filename_1 = output_dirname + "/potential_slices_of_subset_0.h5"
    filename_2 = output_dirname + "/S_matrices_of_subset_0.h5"

    cls_alias_set = (prismatique.sample.PotentialSliceSubsetIDs,
                     prismatique.sample.SMatrixSubsetIDs)
    
    for cls_alias in cls_alias_set:
        if cls_alias == prismatique.sample.PotentialSliceSubsetIDs:
            kwargs = {"filenames": (filename_1,)}
        else:
            kwargs = {"filenames": (filename_2,)}
        sample_specification = cls_alias(**kwargs)

        kwargs = {"sample_specification": sample_specification}
        stem_system_model_params = prismatique.stem.system.ModelParams(**kwargs)
            
        kwargs = {"serializable_rep": stem_system_model_params.pre_serialize()}
        prismatique.stem.system.ModelParams.de_pre_serialize(**kwargs)

    helpers.remove_output_files()

    return None



###########################
## Define error messages ##
###########################
