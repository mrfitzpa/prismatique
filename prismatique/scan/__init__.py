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

# To make directories.
import h5pywrappers.obj

# For checking whether floats are approximately zero.
import emconstants



# Import child modules and packages of current package.
import prismatique.scan.rectangular

# For validating instances of the classes
# :class:`prismatique.sample.ModelParams`,
# :class:`prismatique.sample.PotentialSliceSubsetIDs`,
# :class:`prismatique.sample.SMatrixSubsetIDs`, and
# :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`; and for
# determining the dimensions of sample supercells.
import prismatique.sample



############################
## Authorship information ##
############################

__author__     = "Matthew Fitzpatrick"
__copyright__  = "Copyright 2023"
__credits__    = ["Matthew Fitzpatrick"]
__maintainer__ = "Matthew Fitzpatrick"
__email__      = "mrfitzpa@uvic.ca"
__status__     = "Development"



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["generate_probe_positions",
           "pattern_type",
           "grid_dims_in_units_of_probe_shifts"]



def generate_probe_positions(sample_specification,
                             scan_config=None,
                             filename=None):
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

    Returns
    -------
    probe_positions : `array_like` (`float`, shape=(``num_positions``, ``2``))
        If we let ``num_positions`` be the number of probe positions, then
        ``probe_positions[i][0]`` and ``probe_positions[i][1]`` are the
        :math:`x`- and :math:`y`-coordinates of the ``i`` th probe position in
        units of Å, where ``0<=i<num_positions``.

    """
    accepted_types = (prismatique.sample.ModelParams,
                      prismatique.sample.PotentialSliceSubsetIDs,
                      prismatique.sample.SMatrixSubsetIDs,
                      prismatique.sample.PotentialSliceAndSMatrixSubsetIDs)
    prismatique.sample._check_sample_specification(sample_specification,
                                                   accepted_types)

    temp_ctor_params = {"scan_config": scan_config}
    scan_config = _check_and_convert_scan_config(temp_ctor_params)

    probe_positions = _generate_probe_positions(sample_specification,
                                                scan_config,
                                                filename)
        
    return probe_positions



def _check_and_convert_scan_config(ctor_params):
    # To be used in the module :mod:`prismatique.stem.system` as well, hence the
    # function signature.
    try:
        kwargs = {"obj": ctor_params["scan_config"],
                  "obj_name": "scan_config"}
        scan_config = \
            czekitout.convert.to_real_two_column_numpy_matrix(**kwargs)
    except:
        try:
            scan_config = czekitout.convert.to_str_from_path_like(**kwargs)
        except:
            try:
                temp_ctor_params = \
                    {"rectangular_scan_params": ctor_params["scan_config"]}
                mod_alias = \
                    prismatique.scan.rectangular
                check_and_convert_rectangular_scan_params = \
                    mod_alias._check_and_convert_rectangular_scan_params
                scan_config = \
                    check_and_convert_rectangular_scan_params(temp_ctor_params)
            except:
                raise TypeError(_check_and_convert_scan_config_err_msg_1)
    
    return scan_config



def _generate_probe_positions(sample_specification, scan_config, filename):
    output_filename = filename

    if isinstance(scan_config, str):
        _check_scan_config_file_format(scan_config)
        probe_positions = _load_probe_positions(scan_config)
    elif isinstance(scan_config, prismatique.scan.rectangular.Params):
        kwargs = \
            {"rectangular_scan_params": scan_config,
             "sample_specification": sample_specification}
        probe_positions = \
            prismatique.scan.rectangular._generate_probe_positions(**kwargs)
    else:
        probe_positions = scan_config

    if output_filename is not None:
        _save_probe_positions(probe_positions,
                              sample_specification,
                              output_filename,
                              fold_into_one_supercell=False)
    
    return probe_positions



def _check_scan_config_file_format(scan_config):
    input_filename = scan_config
    _pre_load(input_filename)
    
    try:
        with open(input_filename, 'rb') as file_obj:
            lines = file_obj.readlines()[1:]
    except BaseException as err:
        raise err

    if ((lines[-1] != b"_1\n")
        and (lines[-1] != b"_1\r\n")
        and (lines[-1] != b"-1")):
        err_msg = \
            _check_scan_config_file_format_err_msg_1.format(input_filename)
        raise IOError(err_msg)

    for count, line in enumerate(lines[:-1]):
        try:
            x, y = tuple(float(elem) for elem in line.split())
        except:
            line_num = count + 2
            _check_scan_config_file_format_err_msg_2.format(line_num,
                                                            input_filename)
            raise IOError(err_msg)

    return None



def _pre_load(input_filename):
    try:
        if not pathlib.Path(input_filename).is_file():
            raise FileNotFoundError
    except FileNotFoundError:
        err_msg = h5pywrappers.obj._pre_load_err_msg_1.format(input_filename)
        raise FileNotFoundError(err_msg)
    except PermissionError:
        err_msg = h5pywrappers.obj._pre_load_err_msg_2.format(input_filename)
        raise PermissionError(err_msg)
    except BaseException as err:
        raise err

    return None



def _load_probe_positions(input_filename):
    with open(input_filename, 'rb') as file_obj:
        lines = file_obj.readlines()[1:-1]
    num_lines = len(lines)
    
    probe_positions = np.empty([num_lines, 2])
    for idx, line in enumerate(lines):
        x_coord, y_coord = tuple(float(elem) for elem in line.split())
        probe_positions[idx][0] = x_coord
        probe_positions[idx][1] = y_coord

    return probe_positions



def _save_probe_positions(probe_positions,
                          sample_specification,
                          output_filename,
                          fold_into_one_supercell):
    _pre_save(output_filename)
    
    Delta_X, Delta_Y, _ = \
        prismatique.sample.supercell_dims(sample_specification)

    with open(output_filename, 'w') as file_obj:
        file_obj.write("x y\n")  # Header.
        for x_coord, y_coord in probe_positions:
            if fold_into_one_supercell:
                x_coord = x_coord % Delta_X
                y_coord = y_coord % Delta_Y
            line = str(x_coord) + " " + str(y_coord) + "\n"
            file_obj.write(line)
        file_obj.write("-1")  # End of file.

    return None



def _pre_save(output_filename):
    output_filename = czekitout.convert.to_str_from_path_like(output_filename,
                                                              "filename")
    temp_dir_path, rm_temp_dir = \
        h5pywrappers.obj._mk_parent_dir(output_filename)

    try:
        if not pathlib.Path(output_filename).is_file():
            try:
                with open(output_filename, "w") as file_obj:
                    pass
            except PermissionError:
                err_msg = \
                    h5pywrappers.obj._pre_save_err_msg_1.format(output_filename)
                raise PermissionError(err_msg)
            except BaseException as err:
                raise err
        
            pathlib.Path(output_filename).unlink(missing_ok=True)
        else:
            try:
                with open(output_filename, "a") as file_obj:
                    pass
            except PermissionError:
                err_msg = \
                    h5pywrappers.obj._pre_save_err_msg_2.format(output_filename)
                raise PermissionError(err_msg)
            except BaseException as err:
                raise err
    except BaseException as err:
        if rm_temp_dir:
            shutil.rmtree(temp_dir_path)
        raise err

    if rm_temp_dir:
        shutil.rmtree(temp_dir_path)

    return None



def pattern_type(scan_config=None):
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
    temp_ctor_params = {"scan_config": scan_config}
    scan_config = _check_and_convert_scan_config(temp_ctor_params)
    if isinstance(scan_config, str):
        _check_scan_config_file_format(scan_config)
        
    scan_pattern_type = _pattern_type(scan_config)
    
    return scan_pattern_type



def _pattern_type(scan_config):
    if not isinstance(scan_config, prismatique.scan.rectangular.Params):
        scan_pattern_type = "no underlying rectangular grid"
    else:
        if scan_config.core_attrs["jitter"] == 0:
            scan_pattern_type = "rectangular grid"
        else:
            scan_pattern_type = "jittered rectangular grid"

    return scan_pattern_type



def grid_dims_in_units_of_probe_shifts(sample_specification, scan_config=None):
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
    accepted_types = (prismatique.sample.ModelParams,
                      prismatique.sample.PotentialSliceSubsetIDs,
                      prismatique.sample.SMatrixSubsetIDs,
                      prismatique.sample.PotentialSliceAndSMatrixSubsetIDs)
    prismatique.sample._check_sample_specification(sample_specification,
                                                   accepted_types)

    temp_ctor_params = {"scan_config": scan_config}
    scan_config = _check_and_convert_scan_config(temp_ctor_params)
    grid_dims = _grid_dims_in_units_of_probe_shifts(sample_specification,
                                                    scan_config)

    return grid_dims



def _grid_dims_in_units_of_probe_shifts(sample_specification, scan_config):
    if pattern_type(scan_config) == "no underlying rectangular grid":
        grid_dims = "N/A"
    else:
        Delta_X, Delta_Y, _ = \
            prismatique.sample._supercell_dims(sample_specification)
        
        min_x_probe_coord = Delta_X * scan_config.core_attrs["window"][0]
        max_x_probe_coord = Delta_X * scan_config.core_attrs["window"][1]
        min_y_probe_coord = Delta_Y * scan_config.core_attrs["window"][2]
        max_y_probe_coord = Delta_Y * scan_config.core_attrs["window"][3]

        x_probe_coord_step = scan_config.core_attrs["step_size"][0]
        y_probe_coord_step = scan_config.core_attrs["step_size"][1]
        tol = emconstants.tol

        x_coords = np.arange(min_x_probe_coord,
                             max_x_probe_coord+tol,
                             x_probe_coord_step)
        y_coords = np.arange(min_y_probe_coord,
                             max_y_probe_coord+tol,
                             y_probe_coord_step)

        grid_dims = np.array((len(x_coords), len(y_coords)))

    return grid_dims



###########################
## Define error messages ##
###########################

_check_scan_config_file_format_err_msg_1 = \
    ("The last line in the probe scan configuration file ``{}`` (i.e. the file "
     "specified by the object ``scan_config``) should be '-1'.")
_check_scan_config_file_format_err_msg_2 = \
    ("Line #{} of the probe scan configuration file ``{}`` (i.e. the file "
     "specified by the object ``scan_config``) is not formatted correctly: it "
     "should be of the form 'x y', where 'x' and 'y' are floating-point "
     "numbers.")
