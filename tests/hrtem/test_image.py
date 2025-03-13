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
r"""Contains tests for the module :mod:`prismatique.hrtem.image`.

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



def test_1_of_blank_unprocessed_image_signal(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    for skip_validation_and_conversion in (False, True):
        kwargs = {"sample_specification": \
                  helpers.generate_sample_specification_1(),
                  "skip_validation_and_conversion": \
                  skip_validation_and_conversion}
        prismatique.hrtem.image.blank_unprocessed_image_signal(**kwargs)

    helpers.remove_output_files()

    return None



###########################
## Define error messages ##
###########################
