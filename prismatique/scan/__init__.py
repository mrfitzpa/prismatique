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
r"""For specifying simulation parameters and calculating quantities related to
probe scan patterns.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For checking whether a file exists at a given path and making directories.
import pathlib

# For removing directories.
import shutil



# For general array handling.
import numpy as np

# For validating and converting objects.
import czekitout.check
import czekitout.convert



# Import child modules and packages of current package.
import prismatique.scan.rectangular

# For validating instances of the classes
# :class:`prismatique.sample.ModelParams`,
# :class:`prismatique.sample.PotentialSliceSubsetIDs`,
# :class:`prismatique.sample.SMatrixSubsetIDs`, and
# :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`; and for
# determining the dimensions of sample supercells.
import prismatique.sample

# For recycling helper functions and/or constants.
import prismatique.cbed



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["generate_probe_positions",
           "pattern_type",
           "grid_dims_in_units_of_probe_shifts"]



def _check_and_convert_sample_specification(params):
    module_alias = prismatique.cbed
    func_alias = module_alias._check_and_convert_sample_specification
    sample_specification = func_alias(params)

    return sample_specification



def _check_and_convert_scan_config(params):
    obj_name = "scan_config"

    current_func_name = "_check_and_convert_scan_config"

    try:
        kwargs = \
            {"obj": params[obj_name], "obj_name": obj_name}
        matrix = \
            czekitout.convert.to_real_two_column_numpy_matrix(**kwargs)

        scan_config = tuple()
        for idx, row in enumerate(matrix):
            scan_config += (tuple(row.tolist()),)
    except:
        try:
            scan_config = \
                czekitout.convert.to_str_from_str_like(**kwargs)
        except:
            try:
                params = \
                    {"rectangular_scan_params": params["scan_config"]}
                module_alias = \
                    prismatique.scan.rectangular
                func_alias = \
                    module_alias._check_and_convert_rectangular_scan_params
                scan_config = \
                    func_alias(params)

                del params["rectangular_scan_params"]
            except:
                del params["rectangular_scan_params"]
                err_msg = globals()[current_func_name+"_err_msg_1"]
                raise TypeError(err_msg)
    
    return scan_config



def _check_and_convert_filename(params):
    obj_name = "filename"
    obj = params[obj_name]

    try:
        kwargs = {"obj": obj, "obj_name": obj_name}
        filename = czekitout.convert.to_str_from_str_like(**kwargs)
    except:
        kwargs = {"obj": obj,
                  "obj_name": obj_name,
                  "accepted_types": (str, type(None))}
        czekitout.check.if_instance_of_any_accepted_types(**kwargs)

        filename = obj

    return filename



def _check_and_convert_skip_validation_and_conversion(params):
    module_alias = prismatique.sample
    func_alias = module_alias._check_and_convert_skip_validation_and_conversion
    skip_validation_and_conversion = func_alias(params)

    return skip_validation_and_conversion



_module_alias = \
    prismatique.cbed
_default_scan_config = \
    None
_default_filename = \
    None
_default_skip_validation_and_conversion = \
    _module_alias._default_skip_validation_and_conversion



def generate_probe_positions(sample_specification,
                             scan_config=\
                             _default_scan_config,
                             filename=\
                             _default_filename,
                             skip_validation_and_conversion=\
                             _default_skip_validation_and_conversion):
    r"""Generate the probe positions specified by a given scanning configuration
    for a given sample.

    Parameters
    ----------
    sample_specification : :class:`prismatique.sample.ModelParams` | :class:`prismatique.sample.PotentialSliceSubsetIDs` | :class:`prismatique.sample.SMatrixSubsetIDs` | :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`
        The simulation parameters specifying the sample model. 

        If ``sample_specification`` is of the type
        :class:`prismatique.sample.ModelParams`, then ``sample_specifications``
        specifies sample model parameters that are used to construct the model
        from scratch, i.e. the potential slices for each frozen phonon
        configuration subset are calculated from said model parameters. See the
        documentation for the classes :class:`prismatique.discretization.Params`
        and :class:`prismatique.thermal.Params` for discussions on potential
        slices and frozen phonon configuration subsets respectively. Note that
        of parameters stored in ``sample_specification``, only the following are
        used:

        - sample_specification

          * atomic_coords_filename

          * unit_cell_tiling

        Otherwise, if ``sample_specification`` is an instance of the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs`, the class
        :class:`prismatique.sample.SMatrixSubsetIDs`, or the class
        :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`, then
        ``sample_specification`` specifies a set of files, where each file
        stores either the pre-calculated potential slices or :math:`S`-matrices
        for a frozen phonon configuration subset. See the documentation for the
        aforementioned classes for further discussions on specifying
        pre-calculated objects. See the documentation for the subpackage
        :mod:`prismatique.stem` for a discussion on :math:`S`-matrices.
    scan_config : `array_like` (`float`, shape=(``num_positions``, ``2``)) | :class:`prismatique.scan.rectangular.Params` | `str` | `None`, optional
        If ``scan_config`` is a real-valued two-column matrix, then it specifies
        a set of probe positions, where ``scan_config[i][0]`` and
        ``scan_config[i][1]`` specify respectively the :math:`x`- and
        :math:`y`-coordinates of the probe position indexed by ``i``, in units
        of angstroms, with ``0<=i<num_positions`` and ``num_positions`` being
        the number of probe positions. If ``scan_config`` is of the type
        :class:`prismatique.scan.rectangular.Params`, then it specifies a
        rectangular grid-like pattern of probe positions. See the documentation
        for this class for more details. If ``scan_config`` is a string, then
        ``scan_config`` is a path to a file that specifies a set of probe
        positions. The file must be encoded as ASCII text (UTF-8).  The file
        should be formatted as follows: the first line can be whatever header
        the user desires, i.e. the first line is treated as a comment; each
        subsequent line except for the last is of the form "x y", where "x" and
        "y" are the :math:`x`- and :math:`y`-coordinates of a probe position, in
        units of angstroms; and the last line in the file should be "-1".

        Since periodic boundaries conditions (PBCs) are assumed in the
        :math:`x`- and :math:`y`-dimensions, :math:`x`- and
        :math:`y`-coordinates can take on any real numbers. If ``scan_config``
        is set to `None`, [i.e.  the default value], then the probe is to be
        scanned across the entire unit cell of the simulation, in steps of 0.25
        angstroms in both the :math:`x`- and :math:`y`-directions.
    filename : `str` | `None`, optional
        If ``filename`` is set to a valid filename, then the generated probe
        positions are saved to a file with the filename ``filename``. The file
        is formatted as follows: the first line is "probe positions in units of
        Å"; each subsequent line except for the last is of the form "x y", where
        "x" and "y" are the :math:`x`- and :math:`y`-coordinates of a probe
        position, in units of angstroms; and the last line in the file is "-1".
        Otherwise, if ``filename`` is set to `None` [i.e. the default value],
        then the generated probe positions are not saved to file.
    skip_validation_and_conversion : `bool`, optional
        If ``skip_validation_and_conversion`` is set to ``False``, then
        validations and conversions are performed on the above
        parameters. 

        Otherwise, if ``skip_validation_and_conversion`` is set to ``True``, no
        validations and conversions are performed on the above parameters. This
        option is desired primarily when the user wants to avoid potentially
        expensive validation and/or conversion operations.

    Returns
    -------
    probe_positions : `array_like` (`float`, shape=(``num_positions``, ``2``))
        If we let ``num_positions`` be the number of probe positions, then
        ``probe_positions[i][0]`` and ``probe_positions[i][1]`` are the
        :math:`x`- and :math:`y`-coordinates of the ``i`` th probe position in
        units of Å, where ``0<=i<num_positions``.

    """
    params = locals()

    func_alias = _check_and_convert_skip_validation_and_conversion
    skip_validation_and_conversion = func_alias(params)

    if (skip_validation_and_conversion == False):
        global_symbol_table = globals()
        for param_name in params:
            if param_name == "skip_validation_and_conversion":
                continue
            func_name = "_check_and_convert_" + param_name
            func_alias = global_symbol_table[func_name]
            params[param_name] = func_alias(params)

    kwargs = params
    del kwargs["skip_validation_and_conversion"]
    probe_positions = _generate_probe_positions(**kwargs)

    return probe_positions



def _generate_probe_positions(sample_specification, scan_config, filename):
    output_filename = filename

    if isinstance(scan_config, str):
        _check_scan_config_file_format(scan_config)
        probe_positions = _load_probe_positions(scan_config)
    elif isinstance(scan_config, prismatique.scan.rectangular.Params):
        module_alias = prismatique.scan
        func_alias = module_alias.rectangular._generate_probe_positions
        kwargs = {"rectangular_scan_params": scan_config,
                  "sample_specification": sample_specification}
        probe_positions = func_alias(**kwargs)
    else:
        probe_positions = scan_config

    if output_filename is not None:
        _save_probe_positions(probe_positions,
                              sample_specification,
                              output_filename)
    
    return probe_positions



def _check_scan_config_file_format(scan_config):
    input_filename = scan_config
    _pre_load(input_filename)
    
    with open(input_filename, 'rb') as file_obj:
        lines = file_obj.readlines()[1:]

    current_func_name = "_check_scan_config_file_format"

    if ((lines[-1] != b"_1\n")
        and (lines[-1] != b"_1\r\n")
        and (lines[-1] != b"-1")):
        unformatted_err_msg = globals()[current_func_name+"_err_msg_1"]
        err_msg = unformatted_err_msg.format(input_filename)
        raise IOError(err_msg)

    for count, line in enumerate(lines[:-1]):
        try:
            x, y = tuple(float(elem) for elem in line.split())
        except:
            line_num = count + 2
            unformatted_err_msg = globals()[current_func_name+"_err_msg_2"]
            err_msg = unformatted_err_msg.format(line_num, input_filename)
            raise IOError(err_msg)

    return None



def _pre_load(input_filename):
    current_func_name = "_pre_load"

    try:
        with open(input_filename, 'rb') as file_obj:
            lines = file_obj.readline()
    except FileNotFoundError:
        unformatted_err_msg = globals()[current_func_name+"_err_msg_1"]
        err_msg = unformatted_err_msg.format(input_filename)
        raise FileNotFoundError(err_msg)
    except PermissionError:
        unformatted_err_msg = globals()[current_func_name+"_err_msg_2"]
        err_msg = unformatted_err_msg.format(input_filename)
        raise PermissionError(err_msg)

    return None



def _load_probe_positions(input_filename):
    with open(input_filename, 'rb') as file_obj:
        lines = file_obj.readlines()[1:-1]
    num_lines = len(lines)
    
    probe_positions = tuple()
    for idx, line in enumerate(lines):
        x_coord, y_coord = tuple(float(elem) for elem in line.split())
        probe_positions += ((x_coord, y_coord),)

    return probe_positions



def _save_probe_positions(probe_positions,
                          sample_specification,
                          output_filename):
    _pre_save(output_filename)
    _mk_parent_dir(filename=output_filename)
    
    Delta_X, Delta_Y, _ = \
        prismatique.sample._supercell_dims(sample_specification)

    with open(output_filename, 'w') as file_obj:
        file_obj.write("x y\n")  # Header.
        for x_coord, y_coord in probe_positions:
            line = str(x_coord) + " " + str(y_coord) + "\n"
            file_obj.write(line)
        file_obj.write("-1")  # End of file.

    return None



def _pre_save(output_filename):
    first_new_dir_made = _mk_parent_dir(filename=output_filename)

    current_func_name = "_pre_save"

    if not pathlib.Path(output_filename).is_file():
        try:
            with open(output_filename, "w") as file_obj:
                pass
        except PermissionError:
            unformatted_err_msg = globals()[current_func_name+"_err_msg_1"]
            err_msg = unformatted_err_msg.format(output_filename)
            raise PermissionError(err_msg)
        
        pathlib.Path(output_filename).unlink(missing_ok=True)
    else:
        try:
            with open(output_filename, "a") as file_obj:
                pass
        except PermissionError:
            unformatted_err_msg = globals()[current_func_name+"_err_msg_2"]
            err_msg = unformatted_err_msg.format(output_filename)
            raise PermissionError(err_msg)

    if first_new_dir_made is not None:
        shutil.rmtree(first_new_dir_made)

    return None



def _mk_parent_dir(filename):
    current_func_name = "_mk_parent_dir"

    try:
        parent_dir_path = pathlib.Path(filename).resolve().parent
        temp_dir_path = pathlib.Path(parent_dir_path.root)

        parent_dir_did_not_already_exist = False

        for path_part in parent_dir_path.parts[1:]:
            temp_dir_path = pathlib.Path.joinpath(temp_dir_path, path_part)
            if not temp_dir_path.is_dir():
                parent_dir_did_not_already_exist = True
                break

        pathlib.Path(parent_dir_path).mkdir(parents=True, exist_ok=True)
    except PermissionError:
        err_msg = globals()[current_func_name+"_err_msg_1"].format(filename)
        raise PermissionError(err_msg)

    first_new_dir_made = (temp_dir_path
                          if parent_dir_did_not_already_exist
                          else None)

    return first_new_dir_made



def pattern_type(scan_config=\
                 _default_scan_config,
                 skip_validation_and_conversion=\
                 _default_skip_validation_and_conversion):
    r"""The scan pattern type of a given probe scan pattern.

    Parameters
    ----------
    scan_config : `array_like` (`float`, shape=(``num_positions``, ``2``)) | :class:`prismatique.scan.rectangular.Params` | `str` | `None`, optional
        If ``scan_config`` is a real-valued two-column matrix, then it specifies
        a set of probe positions, where ``scan_config[i][0]`` and
        ``scan_config[i][1]`` specify respectively the :math:`x`- and
        :math:`y`-coordinates of the probe position indexed by ``i``, in units
        of angstroms, with ``0<=i<num_positions`` and ``num_positions`` being
        the number of probe positions. If ``scan_config`` is of the type
        :class:`prismatique.scan.rectangular.Params`, then it specifies a
        rectangular grid-like pattern of probe positions. See the documentation
        for this class for more details. If ``scan_config`` is a string, then
        ``scan_config`` is a path to a file that specifies a set of probe
        positions. The file must be encoded as ASCII text (UTF-8).  The file
        should be formatted as follows: the first line can be whatever header
        the user desires, i.e. the first line is treated as a comment; each
        subsequent line except for the last is of the form "x y", where "x" and
        "y" are the :math:`x`- and :math:`y`-coordinates of a probe position, in
        units of angstroms; and the last line in the file should be "-1".

        Since periodic boundaries conditions (PBCs) are assumed in the
        :math:`x`- and :math:`y`-dimensions, :math:`x`- and
        :math:`y`-coordinates can take on any real numbers. If ``scan_config``
        is set to `None`, [i.e.  the default value], then the probe is to be
        scanned across the entire unit cell of the simulation, in steps of 0.25
        angstroms in both the :math:`x`- and :math:`y`-directions.
    skip_validation_and_conversion : `bool`, optional
        If ``skip_validation_and_conversion`` is set to ``False``, then
        validations and conversions are performed on the above
        parameters. 

        Otherwise, if ``skip_validation_and_conversion`` is set to ``True``, no
        validations and conversions are performed on the above parameters. This
        option is desired primarily when the user wants to avoid potentially
        expensive validation and/or conversion operations.

    Returns
    -------
    scan_pattern_type : "rectangular grid" | "jittered rectangular grid" | "no underlying rectangular grid"
        If ``scan_pattern_type=="rectangular_grid"``, then the probe positions
        making up the scan pattern lie exactly on a regular rectangular grid.
        If ``scan_pattern_type=="jittered rectangular grid"``, then the set of
        probe positions making up the scan pattern lie is equal to the set of
        positions obtained by generating an underlying rectangular grid to which
        a random positional deviation is applied to each point. In this case,
        the pattern is irregular but rectangular grid-like. If
        ``scan_pattern_type=="no underlying rectangular grid"``, then the scan
        pattern is irregular and not rectangular grid-like, i.e. this case
        is different from the previous two.

    """
    params = locals()

    func_alias = _check_and_convert_skip_validation_and_conversion
    skip_validation_and_conversion = func_alias(params)

    if (skip_validation_and_conversion == False):
        global_symbol_table = globals()
        for param_name in params:
            if param_name == "skip_validation_and_conversion":
                continue
            func_name = "_check_and_convert_" + param_name
            func_alias = global_symbol_table[func_name]
            params[param_name] = func_alias(params)

    kwargs = params
    del kwargs["skip_validation_and_conversion"]
    if isinstance(kwargs["scan_config"], str):
        _check_scan_config_file_format(**kwargs)
    scan_pattern_type = _pattern_type(**kwargs)
    
    return scan_pattern_type



def _pattern_type(scan_config):
    if not isinstance(scan_config, prismatique.scan.rectangular.Params):
        scan_pattern_type = "no underlying rectangular grid"
    else:
        scan_config_core_attrs = scan_config.get_core_attrs(deep_copy=False)
        jitter = scan_config_core_attrs["jitter"]
        scan_pattern_type = ("rectangular grid"
                             if (jitter == 0)
                             else "jittered rectangular grid")

    return scan_pattern_type



def grid_dims_in_units_of_probe_shifts(sample_specification,
                                       scan_config=\
                                       _default_scan_config,
                                       skip_validation_and_conversion=\
                                       _default_skip_validation_and_conversion):
    r"""The underlying grid dimensions of a given probe scan pattern, in units
    of probe shifts.

    Parameters
    ----------
    sample_specification : :class:`prismatique.sample.ModelParams` | :class:`prismatique.sample.PotentialSliceSubsetIDs` | :class:`prismatique.sample.SMatrixSubsetIDs` | :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`
        The simulation parameters specifying the sample model. 

        If ``sample_specification`` is of the type
        :class:`prismatique.sample.ModelParams`, then ``sample_specifications``
        specifies sample model parameters that are used to construct the model
        from scratch, i.e. the potential slices for each frozen phonon
        configuration subset are calculated from said model parameters. See the
        documentation for the classes :class:`prismatique.discretization.Params`
        and :class:`prismatique.thermal.Params` for discussions on potential
        slices and frozen phonon configuration subsets respectively. Note that
        of parameters stored in ``sample_specification``, only the following are
        used:

        - sample_specification

          * atomic_coords_filename

          * unit_cell_tiling

        Otherwise, if ``sample_specification`` is an instance of the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs`, the class
        :class:`prismatique.sample.SMatrixSubsetIDs`, or the class
        :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`, then
        ``sample_specification`` specifies a set of files, where each file
        stores either the pre-calculated potential slices or :math:`S`-matrices
        for a frozen phonon configuration subset. See the documentation for the
        aforementioned classes for further discussions on specifying
        pre-calculated objects. See the documentation for the subpackage
        :mod:`prismatique.stem` for a discussion on :math:`S`-matrices.
    scan_config : `array_like` (`float`, shape=(``num_positions``, ``2``)) | :class:`prismatique.scan.rectangular.Params` | `str` | `None`, optional
        If ``scan_config`` is a real-valued two-column matrix, then it specifies
        a set of probe positions, where ``scan_config[i][0]`` and
        ``scan_config[i][1]`` specify respectively the :math:`x`- and
        :math:`y`-coordinates of the probe position indexed by ``i``, in units
        of angstroms, with ``0<=i<num_positions`` and ``num_positions`` being
        the number of probe positions. If ``scan_config`` is of the type
        :class:`prismatique.scan.rectangular.Params`, then it specifies a
        rectangular grid-like pattern of probe positions. See the documentation
        for this class for more details. If ``scan_config`` is a string, then
        ``scan_config`` is a path to a file that specifies a set of probe
        positions. The file must be encoded as ASCII text (UTF-8).  The file
        should be formatted as follows: the first line can be whatever header
        the user desires, i.e. the first line is treated as a comment; each
        subsequent line except for the last is of the form "x y", where "x" and
        "y" are the :math:`x`- and :math:`y`-coordinates of a probe position, in
        units of angstroms; and the last line in the file should be "-1".

        Since periodic boundaries conditions (PBCs) are assumed in the
        :math:`x`- and :math:`y`-dimensions, :math:`x`- and
        :math:`y`-coordinates can take on any real numbers. If ``scan_config``
        is set to `None`, [i.e.  the default value], then the probe is to be
        scanned across the entire unit cell of the simulation, in steps of 0.25
        angstroms in both the :math:`x`- and :math:`y`-directions.
    skip_validation_and_conversion : `bool`, optional
        If ``skip_validation_and_conversion`` is set to ``False``, then
        validations and conversions are performed on the above
        parameters. 

        Otherwise, if ``skip_validation_and_conversion`` is set to ``True``, no
        validations and conversions are performed on the above parameters. This
        option is desired primarily when the user wants to avoid potentially
        expensive validation and/or conversion operations.

    Returns
    -------
    grid_dims : "N/A" | `array_like` (`float`, shape=(``2``))
        If ``prismatique.scan.pattern_type(scan_config) == "no underlying
        rectangular grid"``, then ``grid_dimensions_in_units_of_probe_shifts ==
        "N/A"``, indicating that there is no notion of grid dimensions that is
        applicable to the scan pattern used. Otherwise, if
        ``prismatique.scan.pattern_type(scan_config) != "no underlying
        rectangular grid"``, then ``grid_dims[0]`` and ``grid_dims[1]`` are the
        number of probe positions along the :math:`x`- and :math:`y`-dimensions
        respectively of the underlying rectangular grid of the scanning pattern.

    """
    params = locals()

    func_alias = _check_and_convert_skip_validation_and_conversion
    skip_validation_and_conversion = func_alias(params)

    if (skip_validation_and_conversion == False):
        global_symbol_table = globals()
        for param_name in params:
            if param_name == "skip_validation_and_conversion":
                continue
            func_name = "_check_and_convert_" + param_name
            func_alias = global_symbol_table[func_name]
            params[param_name] = func_alias(params)

    kwargs = params
    del kwargs["skip_validation_and_conversion"]
    grid_dims = _grid_dims_in_units_of_probe_shifts(**kwargs)

    return grid_dims



def _grid_dims_in_units_of_probe_shifts(sample_specification, scan_config):
    if pattern_type(scan_config) == "no underlying rectangular grid":
        grid_dims = "N/A"
    else:
        Delta_X, Delta_Y, _ = \
            prismatique.sample._supercell_dims(sample_specification)

        scan_config_core_attrs = scan_config.get_core_attrs(deep_copy=False)
        
        min_x_probe_coord = Delta_X * scan_config_core_attrs["window"][0]
        max_x_probe_coord = Delta_X * scan_config_core_attrs["window"][1]
        min_y_probe_coord = Delta_Y * scan_config_core_attrs["window"][2]
        max_y_probe_coord = Delta_Y * scan_config_core_attrs["window"][3]

        x_probe_coord_step = scan_config_core_attrs["step_size"][0]
        y_probe_coord_step = scan_config_core_attrs["step_size"][1]
        tol = 1e-10

        x_coords = np.arange(min_x_probe_coord,
                             max_x_probe_coord+tol,
                             x_probe_coord_step)
        y_coords = np.arange(min_y_probe_coord,
                             max_y_probe_coord+tol,
                             y_probe_coord_step)

        grid_dims = (len(x_coords), len(y_coords))

    return grid_dims



###########################
## Define error messages ##
###########################

_check_and_convert_scan_config_err_msg_1 = \
    ("The object ``scan_config`` must be either a real two-column matrix, "
     "a string, or an instance of the class "
     ":class:`prismatique.scan.rectangular.Params`.")

_check_scan_config_file_format_err_msg_1 = \
    ("The last line in the probe scan configuration file ``{}`` (i.e. the file "
     "specified by the object ``scan_config``) should be '-1'.")
_check_scan_config_file_format_err_msg_2 = \
    ("Line #{} of the probe scan configuration file ``{}`` (i.e. the file "
     "specified by the object ``scan_config``) is not formatted correctly: it "
     "should be of the form 'x y', where 'x' and 'y' are floating-point "
     "numbers.")

_pre_load_err_msg_1 = \
    ("No file exists at the file path ``'{}'``.")
_pre_load_err_msg_2 = \
    ("Cannot access the file path ``'{}'`` because of insufficient "
     "permissions.")

_pre_save_err_msg_1 = \
    _pre_load_err_msg_2
_pre_save_err_msg_2 = \
    ("Cannot write to the file at the file path ``'{}'`` because of "
     "insufficient permissions.")

_mk_parent_dir_err_msg_1 = \
    _pre_load_err_msg_2
