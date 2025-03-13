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
r"""Contains tests for the module :mod:`prismatique.scan`.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For operations related to unit tests.
import pytest

# For setting file permissions.
import os

# For performing operations on file and directory paths.
import pathlib

# For removing directories.
import shutil



# For simulating electron microscopy experiments using Multislice algorithms.
import prismatique



##################################
## Define classes and functions ##
##################################



def test_1_of_generate_probe_positions(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()
    helpers.generate_invalid_scan_config_file_1()
    helpers.generate_invalid_scan_config_file_2()

    sample_specification = helpers.generate_sample_specification_1()

    for skip_validation_and_conversion in (False, True):
        kwargs = {"sample_specification": \
                  sample_specification,
                  "scan_config": \
                  ((0, 0),),
                  "skip_validation_and_conversion": \
                  skip_validation_and_conversion}
        prismatique.scan.generate_probe_positions(**kwargs)

    with pytest.raises(TypeError) as err_info:
        kwargs = {"sample_specification": sample_specification,
                  "scan_config": slice(None),
                  "skip_validation_and_conversion": False}
        prismatique.scan.generate_probe_positions(**kwargs)

    filenames = (helpers.generate_scan_config_filename_2(),
                 helpers.generate_scan_config_filename_3())

    for filename in filenames:
        with pytest.raises(IOError) as err_info:
            kwargs = {"sample_specification": sample_specification,
                      "scan_config": filename}
            prismatique.scan.generate_probe_positions(**kwargs)

    os.chmod(filenames[1], 0o111)
    with pytest.raises(PermissionError) as err_info:
        prismatique.scan.generate_probe_positions(**kwargs)
    os.chmod(filenames[1], 0o711)

    with pytest.raises(FileNotFoundError) as err_info:
        kwargs = {"sample_specification": sample_specification,
                  "scan_config": "foobar"}
        prismatique.scan.generate_probe_positions(**kwargs)

    helpers.remove_output_files()

    return None



def test_2_of_generate_probe_positions(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    dirname_1 = helpers.generate_output_dirname_1()
    dirname_2 = dirname_1 + "/probe_position_output"
    filename = dirname_2 + "/scan_config_file.txt"

    kwargs = {"sample_specification": helpers.generate_sample_specification_1(),
              "scan_config": ((0, 0),),
              "filename": filename}
    prismatique.scan.generate_probe_positions(**kwargs)

    os.chmod(filename, 0o111)
    with pytest.raises(PermissionError) as err_info:
        prismatique.scan.generate_probe_positions(**kwargs)
    os.chmod(filename, 0o711)

    pathlib.Path(filename).unlink(missing_ok=True)

    os.chmod(dirname_2, 0o111)
    with pytest.raises(PermissionError) as err_info:
        prismatique.scan.generate_probe_positions(**kwargs)
    os.chmod(dirname_2, 0o711)

    shutil.rmtree(dirname_2, ignore_errors=True)

    os.chmod(dirname_1, 0o111)
    with pytest.raises(PermissionError) as err_info:
        prismatique.scan.generate_probe_positions(**kwargs)
    os.chmod(dirname_1, 0o711)

    prismatique.scan.generate_probe_positions(**kwargs)

    helpers.remove_output_files()

    return None



def test_1_of_pattern_type(helpers):
    for skip_validation_and_conversion in (False, True):
        kwargs = {"scan_config": \
                  ((0, 0),),
                  "skip_validation_and_conversion": \
                  skip_validation_and_conversion}
        prismatique.scan.pattern_type(**kwargs)

    return None



def test_1_of_grid_dims_in_units_of_probe_shifts(helpers):
    helpers.remove_output_files()

    helpers.generate_atomic_coords_file_1()

    sample_specification = helpers.generate_sample_specification_1()

    for skip_validation_and_conversion in (False, True):
        kwargs = {"sample_specification": \
                  sample_specification,
                  "scan_config": \
                  ((0, 0),),
                  "skip_validation_and_conversion": \
                  skip_validation_and_conversion}
        prismatique.scan.grid_dims_in_units_of_probe_shifts(**kwargs)

    helpers.remove_output_files()

    return None



###########################
## Define error messages ##
###########################
