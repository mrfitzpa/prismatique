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
the modelling of the sample.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For performing operations on file and directory paths.
import pathlib

# For pattern matching.
import re

# For explicit garbage collection.
import gc



# For general array handling.
import numpy as np

# For validating and converting objects.
import czekitout.check
import czekitout.convert

# For defining classes that support enforced validation, updatability,
# pre-serialization, and de-serialization.
import fancytypes

# For calculating the electron beam wavelength and validating attributes of the
# class :class:`embeam.stem.probe.ModelParams`.
import embeam

# For loading objects from and saving objects to HDF5 files.
import h5pywrappers

# For creating a :obj:`pyprismatic.Metadata` object that is responsible for
# running the ``prismatic`` simulation.
import pyprismatic



# For validating, pre-serializing, and de-pre-serializing instances of the
# classes :class:`prismatique.discretization.Params`, and
# :class:`prismatique.thermal.Params`.
import prismatique.discretization
import prismatique.thermal

# For validating the class :class:`prismatique.worker.Params`.
import prismatique.worker



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["ModelParams",
           "check_atomic_coords_file_format",
           "PotentialSliceSubsetIDs",
           "SMatrixSubsetIDs",
           "PotentialSliceAndSMatrixSubsetIDs",
           "are_valid_precalculated_subset_ids",
           "unit_cell_dims",
           "supercell_dims",
           "supercell_xy_dims_in_pixels",
           "supercell_lateral_pixel_size",
           "supercell_slice_thickness",
           "num_slices",
           "num_frozen_phonon_config_subsets",
           "num_frozen_phonon_configs_in_subset",
           "total_num_frozen_phonon_configs",
           "S_matrix_k_xy_vectors",
           "potential_slice_subset_data_size",
           "potential_slice_set_data_size",
           "S_matrix_subset_data_size",
           "S_matrix_set_data_size",
           "generate_potential_slices",
           "generate_S_matrices"]



def _check_and_convert_atomic_coords_filename(params):
    obj_name = "atomic_coords_filename"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    atomic_coords_filename = czekitout.convert.to_str_from_str_like(**kwargs)
    
    return atomic_coords_filename



def _pre_serialize_atomic_coords_filename(atomic_coords_filename):
    obj_to_pre_serialize = atomic_coords_filename
    serializable_rep = atomic_coords_filename

    return serializable_rep



def _de_pre_serialize_atomic_coords_filename(serializable_rep):
    atomic_coords_filename = serializable_rep

    return atomic_coords_filename



def _check_and_convert_unit_cell_tiling(params):
    obj_name = "unit_cell_tiling"

    current_func_name = "_check_and_convert_unit_cell_tiling"

    try:
        kwargs = {"obj": params[obj_name], "obj_name": obj_name}
        unit_cell_tiling = czekitout.convert.to_tuple_of_positive_ints(**kwargs)
        if len(unit_cell_tiling) != 3:
            raise
    except:
        err_msg = globals()[current_func_name+"_err_msg_1"]
        raise TypeError(err_msg)
    
    return unit_cell_tiling



def _pre_serialize_unit_cell_tiling(unit_cell_tiling):
    obj_to_pre_serialize = unit_cell_tiling
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_unit_cell_tiling(serializable_rep):
    unit_cell_tiling = serializable_rep

    return unit_cell_tiling



def _check_and_convert_discretization_params(params):
    module_alias = prismatique.discretization
    func_alias = module_alias._check_and_convert_discretization_params
    discretization_params = func_alias(params)
    
    return discretization_params



def _pre_serialize_discretization_params(discretization_params):
    obj_to_pre_serialize = discretization_params
    module_alias = prismatique.discretization
    func_alias = module_alias._pre_serialize_discretization_params
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_discretization_params(serializable_rep):
    module_alias = prismatique.discretization
    func_alias = module_alias._de_pre_serialize_discretization_params
    discretization_params = func_alias(serializable_rep)

    return discretization_params



def _check_and_convert_atomic_potential_extent(params):
    obj_name = "atomic_potential_extent"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    atomic_potential_extent = czekitout.convert.to_positive_float(**kwargs)
    
    return atomic_potential_extent



def _pre_serialize_atomic_potential_extent(atomic_potential_extent):
    obj_to_pre_serialize = atomic_potential_extent
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_atomic_potential_extent(serializable_rep):
    atomic_potential_extent = serializable_rep

    return atomic_potential_extent



def _check_and_convert_thermal_params(params):
    module_alias = prismatique.thermal
    func_alias = module_alias._check_and_convert_thermal_params
    thermal_params = func_alias(params)
    
    return thermal_params



def _pre_serialize_thermal_params(thermal_params):
    obj_to_pre_serialize = thermal_params
    module_alias = prismatique.thermal
    func_alias = module_alias._pre_serialize_thermal_params
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_thermal_params(serializable_rep):
    module_alias = prismatique.thermal
    func_alias = module_alias._de_pre_serialize_thermal_params
    thermal_params = func_alias(serializable_rep)

    return thermal_params



_module_alias_1 = \
    prismatique.discretization
_module_alias_2 = \
    prismatique.thermal
_default_unit_cell_tiling = \
    3*(1,)
_default_discretization_params = \
    _module_alias_1._default_discretization_params
_default_atomic_potential_extent = \
    3
_default_thermal_params = \
    _module_alias_2._default_thermal_params
_default_skip_validation_and_conversion = \
    _module_alias_1._default_skip_validation_and_conversion



class ModelParams(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related to the modelling of the sample.

    Parameters
    ----------
    atomic_coords_filename : `str`
        See the documentation for the classes
        :class:`prismatique.discretization.Params` and
        :class:`prismatique.thermal.Params` for context: the former discusses
        sample unit cells and supercells; and the latter discusses the effective
        root-mean-squared displacements of atoms, which we denote by
        :math:`u_{i,\text{rms}}`.

        ``atomic_coords_filename`` is a relative or absolute path to a file that
        specifies the zero-temperature expectation value of the atomic
        coordinates for each atom in the unit cell of the sample, along with
        :math:`\frac{1}{\sqrt{3}}u_{i,\text{rms}}` for each atom, where
        :math:`u_{i,\text{rms}}` was introduced in Eq. :eq:`u_i_rms`. The file
        must be encoded as ASCII text (UTF-8).  The file should be formatted as
        follows: the first line can be whatever header the user desires,
        i.e. the first line is treated as a comment; the second line is of the
        form "a b c", where "a", "b", and "c" are the :math:`x`-, :math:`y`-,
        and :math:`z`-dimensions of the unit cell in angstroms; each subsequent
        line except for the last specifies attributes of an atom and is of the
        form "Z x y z' occ u_x_rms", where "Z" is the atomic number of said
        atom, "x", "y" and "z'" are the :math:`x`-, :math:`y`-, and
        :math:`z^{\prime}`-coordinates respectively of said atom in angstroms
        [the :math:`z^{\prime}`-coordinate of an atom is related to the
        :math:`z`-coordinate of the same atom by :math:`z=-z^{\prime}+\Delta Z`,
        with :math:`\Delta Z` being the :math:`z`-dimension of the sample's
        supercell], "occ" is the likelihood that said atom exists, and "u_x_rms"
        is :math:`\frac{1}{\sqrt{3}}u_{i,\text{rms}}` for said atom in
        angstroms; and the last line in the file should be "-1". Example files
        can be found in the ``examples`` directory of the ``prismatique``
        repository. Note that "occ" is ignored in ``prismatique``.
    unit_cell_tiling : `array_like` (`int`, shape=(``3``,)), optional
        ``unit_cell_tiling[0]``, ``unit_cell_tiling[1]``, and
        ``unit_cell_tiling[3]`` specify the number of times to tile the unit
        cell in the :math:`x`-, :math:`y`-, and :math:`z`- directions
        respectively to form the sample supercell. See the documentation for the
        class :class:`prismatique.discretization.Params` for a discussion on the
        sample's unit cell and supercell.
    discretization_params : :class:`prismatique.discretization.Params` | `None`, optional
        The simulation parameters related to the discretization of real-space
        and Fourier/:math:`k`-space. See the documentation for the class
        :class:`prismatique.discretization.Params` for a discussion on said
        parameters. If ``discretization_params`` is set to `None` [i.e. the 
        default value], then the aforementioned simulation parameters are set to
        default values.
    atomic_potential_extent : `float`, optional
        See the documentation for the class
        :class:`prismatique.discretization.Params` for
        context. ``atomic_potential_extent`` specifies the atomic potential
        extent :math:`r_{\text{max}}` in angstroms, where

        .. math ::
            V_{j}^{\left(\text{atom}\right)}
            \left(\left|\mathbf{r}-\mathbf{r}_{j}\right|>r_{\text{max}}\right)
            = 0,
            :label: potential_extent_1

        and

        .. math ::
            \tilde{V}_{j}^{\left(\text{atom}\right)}
            \left(\left|\left\{\mathbf{r}-\mathbf{r}_{j}\right\}_{xy}\right|
            >r_{\text{max}}\right)
            = 0.
            :label: potential_extent_2
    thermal_params : :class:`prismatique.thermal.Params` | `None`, optional
        The simulation parameters related to the thermal properties of the
        sample and its environment. See the documentation for the class
        :class:`prismatique.thermal.Params` for a discussion on said parameters.
        If ``thermal_params`` is set to `None` [i.e. the default value], then
        the aforementioned simulation parameters are set to default values.
    skip_validation_and_conversion : `bool`, optional
        Let ``validation_and_conversion_funcs`` and ``core_attrs`` denote the
        attributes :attr:`~fancytypes.Checkable.validation_and_conversion_funcs`
        and :attr:`~fancytypes.Checkable.core_attrs` respectively, both of which
        being `dict` objects.

        Let ``params_to_be_mapped_to_core_attrs`` denote the `dict`
        representation of the constructor parameters excluding the parameter
        ``skip_validation_and_conversion``, where each `dict` key ``key`` is a
        different constructor parameter name, excluding the name
        ``"skip_validation_and_conversion"``, and
        ``params_to_be_mapped_to_core_attrs[key]`` would yield the value of the
        constructor parameter with the name given by ``key``.

        If ``skip_validation_and_conversion`` is set to ``False``, then for each
        key ``key`` in ``params_to_be_mapped_to_core_attrs``,
        ``core_attrs[key]`` is set to ``validation_and_conversion_funcs[key]
        (params_to_be_mapped_to_core_attrs)``.

        Otherwise, if ``skip_validation_and_conversion`` is set to ``True``,
        then ``core_attrs`` is set to
        ``params_to_be_mapped_to_core_attrs.copy()``. This option is desired
        primarily when the user wants to avoid potentially expensive deep copies
        and/or conversions of the `dict` values of
        ``params_to_be_mapped_to_core_attrs``, as it is guaranteed that no
        copies or conversions are made in this case.

    """
    ctor_param_names = ("atomic_coords_filename",
                        "unit_cell_tiling",
                        "discretization_params",
                        "atomic_potential_extent",
                        "thermal_params")
    kwargs = {"namespace_as_dict": globals(),
              "ctor_param_names": ctor_param_names}
    
    _validation_and_conversion_funcs_ = \
        fancytypes.return_validation_and_conversion_funcs(**kwargs)
    _pre_serialization_funcs_ = \
        fancytypes.return_pre_serialization_funcs(**kwargs)
    _de_pre_serialization_funcs_ = \
        fancytypes.return_de_pre_serialization_funcs(**kwargs)

    del ctor_param_names, kwargs

    
    
    def __init__(self,
                 atomic_coords_filename,
                 unit_cell_tiling=\
                 _default_unit_cell_tiling,
                 discretization_params=\
                 _default_discretization_params,
                 atomic_potential_extent=\
                 _default_atomic_potential_extent,
                 thermal_params=\
                 _default_thermal_params,
                 skip_validation_and_conversion=\
                 _default_skip_validation_and_conversion):
        ctor_params = {key: val
                       for key, val in locals().items()
                       if (key not in ("self", "__class__"))}
        kwargs = ctor_params
        kwargs["skip_cls_tests"] = True
        fancytypes.PreSerializableAndUpdatable.__init__(self, **kwargs)

        return None



    @classmethod
    def get_validation_and_conversion_funcs(cls):
        validation_and_conversion_funcs = \
            cls._validation_and_conversion_funcs_.copy()

        return validation_and_conversion_funcs


    
    @classmethod
    def get_pre_serialization_funcs(cls):
        pre_serialization_funcs = \
            cls._pre_serialization_funcs_.copy()

        return pre_serialization_funcs


    
    @classmethod
    def get_de_pre_serialization_funcs(cls):
        de_pre_serialization_funcs = \
            cls._de_pre_serialization_funcs_.copy()

        return de_pre_serialization_funcs



def _check_and_convert_sample_model_params(params):
    obj_name = "sample_model_params"
    obj = params[obj_name]

    accepted_types = (ModelParams,)
    
    kwargs = {"obj": obj,
              "obj_name": obj_name,
              "accepted_types": accepted_types}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    kwargs = obj.get_core_attrs(deep_copy=False)
    sample_model_params = accepted_types[0](**kwargs)

    return sample_model_params



def _pre_serialize_sample_model_params(sample_model_params):
    obj_to_pre_serialize = sample_model_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_sample_model_params(serializable_rep):
    sample_model_params = ModelParams.de_pre_serialize(serializable_rep)
    
    return sample_model_params



def check_atomic_coords_file_format(atomic_coords_filename):
    r"""Check whether a given atomic coordinates file is of the correct format.

    An exception is raised if the file is invalid.

    Parameters
    ----------
    atomic_coords_filename : `str`
        See the documentation for the classes
        :class:`prismatique.discretization.Params` and
        :class:`prismatique.thermal.Params` for context: the former discusses
        sample unit cells and supercells; and the latter discusses the effective
        root-mean-squared displacements of atoms, which we denote by
        :math:`u_{i,\text{rms}}`.

        ``atomic_coords_filename`` is a relative or absolute path to a file that
        specifies the zero-temperature expectation value of the atomic
        coordinates for each atom in the unit cell of the sample, along with
        :math:`\frac{1}{\sqrt{3}}u_{i,\text{rms}}` for each atom, where
        :math:`u_{i,\text{rms}}` was introduced in Eq. :eq:`u_i_rms`. The file
        must be encoded as ASCII text (UTF-8).  The file should be formatted as
        follows: the first line can be whatever header the user desires,
        i.e. the first line is treated as a comment; the second line is of the
        form "a b c", where "a", "b", and "c" are the :math:`x`-, :math:`y`-,
        and :math:`z`-dimensions of the unit cell in angstroms; each subsequent
        line except for the last specifies attributes of an atom and is of the
        form "Z x y z' occ u_x_rms", where "Z" is the atomic number of said
        atom, "x", "y" and "z'" are the :math:`x`-, :math:`y`-, and
        :math:`z^{\prime}`-coordinates respectively of said atom in angstroms
        [the :math:`z^{\prime}`-coordinate of an atom is related to the
        :math:`z`-coordinate of the same atom by :math:`z=-z^{\prime}+\Delta Z`,
        with :math:`\Delta Z` being the :math:`z`-dimension of the sample's
        supercell], "occ" is the likelihood that said atom exists, and "u_x_rms"
        is :math:`\frac{1}{\sqrt{3}}u_{i,\text{rms}}` for said atom in
        angstroms; and the last line in the file should be "-1". Example files
        can be found in the ``examples`` directory of the ``prismatique``
        repository. Note that "occ" is ignored in ``prismatique``.

    Returns
    -------

    """
    kwargs = {"obj": atomic_coords_filename,
              "obj_name": "atomic_coords_filename"}
    filename = czekitout.convert.to_str_from_str_like(**kwargs)

    current_func_name = "check_atomic_coords_file_format"

    try:
        with open(filename, 'rb') as file_obj:
            lines = file_obj.readlines()[1:]
    except FileNotFoundError:
        unformatted_err_msg = globals()["_"+current_func_name+"_err_msg_1"]
        err_msg = unformatted_err_msg.format(filename)
        raise FileNotFoundError(err_msg)

    if ((lines[-1] != b"_1\n")
        and (lines[-1] != b"_1\r\n")
        and (lines[-1] != b"-1")):
        unformatted_err_msg = globals()["_"+current_func_name+"_err_msg_2"]
        err_msg = unformatted_err_msg.format(filename)
        raise IOError(err_msg)

    try:
        abc = tuple(float(elem) for elem in lines[0].split())
        czekitout.check.if_positive_float_seq(abc, "abc")
    except:
        unformatted_err_msg = globals()["_"+current_func_name+"_err_msg_3"]
        err_msg = unformatted_err_msg.format(filename)
        raise IOError(err_msg)

    try:
        for count, line in enumerate(lines[1:-1]):
            Z, x, y, z, occ, u_x_rms = \
                tuple(float(elem) for elem in line.split())

            # Atomic potentials have only been calculate for 1<=Z<=103.
            if ((Z > 103) or (not Z.is_integer())
                or (u_x_rms < 0) or (occ < 0) or (occ > 1)):
                raise
    except:
        line_num = count + 2
        unformatted_err_msg = globals()["_"+current_func_name+"_err_msg_4"]
        err_msg = unformatted_err_msg.format(line_num, filename)
        raise IOError(err_msg)

    return None
    


def _check_and_convert_filenames(params):
    obj_name = "filenames"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    filenames = czekitout.convert.to_tuple_of_strs(**kwargs)

    current_func_name = "_check_and_convert_filenames"

    if len(filenames) == 0:
        err_msg = globals()[current_func_name+"_err_msg_1"]
        raise ValueError(err_msg)
    
    return filenames



def _pre_serialize_filenames(filenames):
    obj_to_pre_serialize = filenames
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_filenames(serializable_rep):
    filenames = serializable_rep

    return filenames



def _check_and_convert_interpolation_factors(params):
    module_alias = prismatique.discretization
    func_alias = module_alias._check_and_convert_interpolation_factors
    interpolation_factors = func_alias(params)
    
    return interpolation_factors



def _pre_serialize_interpolation_factors(interpolation_factors):
    obj_to_pre_serialize = interpolation_factors
    module_alias = prismatique.discretization
    func_alias = module_alias._pre_serialize_interpolation_factors
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_interpolation_factors(serializable_rep):
    module_alias = prismatique.discretization
    func_alias = module_alias._de_pre_serialize_interpolation_factors
    interpolation_factors = func_alias(serializable_rep)

    return interpolation_factors



def _check_and_convert_max_num_frozen_phonon_configs_per_subset(params):
    obj_name = "max_num_frozen_phonon_configs_per_subset"

    kwargs = \
        {"obj": params[obj_name], "obj_name": obj_name}
    max_num_frozen_phonon_configs_per_subset = \
        czekitout.convert.to_positive_int(**kwargs)

    current_func_name = ("_check_and_convert"
                         "_max_num_frozen_phonon_configs_per_subset")

    if max_num_frozen_phonon_configs_per_subset < 2:
        err_msg = globals()[current_func_name+"_err_msg_1"]
        raise ValueError(err_msg)
    
    return max_num_frozen_phonon_configs_per_subset



def _pre_serialize_max_num_frozen_phonon_configs_per_subset(
        max_num_frozen_phonon_configs_per_subset):
    obj_to_pre_serialize = max_num_frozen_phonon_configs_per_subset
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_max_num_frozen_phonon_configs_per_subset(
        serializable_rep):
    max_num_frozen_phonon_configs_per_subset = serializable_rep

    return max_num_frozen_phonon_configs_per_subset



_module_alias = \
    prismatique.discretization
_default_interpolation_factors = \
    _module_alias._default_interpolation_factors
_default_max_num_frozen_phonon_configs_per_subset = \
    1000



class PotentialSliceSubsetIDs(fancytypes.PreSerializableAndUpdatable):
    r"""A parameter set specifying a set of files, where each file stores the
    pre-calculated potential slices for a frozen phonon configuration subset.

    See the documentation for the classes
    :class:`prismatique.discretization.Params` and
    :class:`prismatique.thermal.Params` for discussions on potential slices and
    frozen phonon configuration subsets respectively.

    HDF5 files storing pre-calculated potential slices can be generated using
    the functions :func:`prismatique.sample.generate_potential_slices`,
    :func:`prismatique.stem.sim.run`, and :func:`prismatique.hrtem.sim.run`.

    Parameters
    ----------
    filenames : `array_like` (`str`, ndim=1)
        For each integer ``i`` satisfying ``0<=i<len(filenames)``,
        ``filenames[i]`` specifies the filename of the pre-calculated potential
        slices for the ``i`` th frozen phonon configuration subset. The
        pre-calculated objects for each frozen phonon configuration subset are
        assumed to correspond to the same sample. In particular, each
        pre-calculated potential slice subset must have been calculated using
        the same sample supercell dimensions :math:`\left(\Delta X, \Delta Y,
        \Delta Z\right)`, slice thickness :math:`\delta z`, and potential slice
        pixel sizes :math:`\left(\Delta x, \Delta y\right)` in the :math:`x`-
        and :math:`y`-directions. See the documentation for the class
        :class:`prismatique.discretization.Params` for discussions on the
        aforementioned quantities. Lastly, note that the total number of frozen
        phonon configuration subsets is ``len(filenames)``, and must be greater
        than or equal to 1.
    interpolation_factors : `array_like` (`int`, shape=(``2``,)), optional
        ``interpolation_factors[0]`` and ``interpolation_factors[1]`` are the
        interpolation factors :math:`f_x` and :math:`f_y` respectively, which
        are discussed in the documentation for the class
        :class:`prismatique.discretization.Params`. Note that
        ``interpolation_factors`` must be set to ``(1, 1)`` if using the
        multislice algorithm for STEM simulations. Otherwise,
        ``interpolation_factors`` must be set to a pair of positive integers.
    max_num_frozen_phonon_configs_per_subset : `int`, optional
        The maximum number of frozen phonon configurations to be loaded, per
        frozen phonon configuration subset. Due to a bug in ``prismatic``, the
        minimum possible value of ``max_num_frozen_phonon_configs_per_subset``
        is 2.
    skip_validation_and_conversion : `bool`, optional
        Let ``validation_and_conversion_funcs`` and ``core_attrs`` denote the
        attributes :attr:`~fancytypes.Checkable.validation_and_conversion_funcs`
        and :attr:`~fancytypes.Checkable.core_attrs` respectively, both of which
        being `dict` objects.

        Let ``params_to_be_mapped_to_core_attrs`` denote the `dict`
        representation of the constructor parameters excluding the parameter
        ``skip_validation_and_conversion``, where each `dict` key ``key`` is a
        different constructor parameter name, excluding the name
        ``"skip_validation_and_conversion"``, and
        ``params_to_be_mapped_to_core_attrs[key]`` would yield the value of the
        constructor parameter with the name given by ``key``.

        If ``skip_validation_and_conversion`` is set to ``False``, then for each
        key ``key`` in ``params_to_be_mapped_to_core_attrs``,
        ``core_attrs[key]`` is set to ``validation_and_conversion_funcs[key]
        (params_to_be_mapped_to_core_attrs)``.

        Otherwise, if ``skip_validation_and_conversion`` is set to ``True``,
        then ``core_attrs`` is set to
        ``params_to_be_mapped_to_core_attrs.copy()``. This option is desired
        primarily when the user wants to avoid potentially expensive deep copies
        and/or conversions of the `dict` values of
        ``params_to_be_mapped_to_core_attrs``, as it is guaranteed that no
        copies or conversions are made in this case.

    """
    ctor_param_names = ("filenames",
                        "interpolation_factors",
                        "max_num_frozen_phonon_configs_per_subset")
    kwargs = {"namespace_as_dict": globals(),
              "ctor_param_names": ctor_param_names}
    
    _validation_and_conversion_funcs_ = \
        fancytypes.return_validation_and_conversion_funcs(**kwargs)
    _pre_serialization_funcs_ = \
        fancytypes.return_pre_serialization_funcs(**kwargs)
    _de_pre_serialization_funcs_ = \
        fancytypes.return_de_pre_serialization_funcs(**kwargs)

    del ctor_param_names, kwargs

    
    
    def __init__(self,
                 filenames,
                 interpolation_factors=\
                 _default_interpolation_factors,
                 max_num_frozen_phonon_configs_per_subset=\
                 _default_max_num_frozen_phonon_configs_per_subset,
                 skip_validation_and_conversion=\
                 _default_skip_validation_and_conversion):
        ctor_params = {key: val
                       for key, val in locals().items()
                       if (key not in ("self", "__class__"))}
        kwargs = ctor_params
        kwargs["skip_cls_tests"] = True
        fancytypes.PreSerializableAndUpdatable.__init__(self, **kwargs)

        return None



    @classmethod
    def get_validation_and_conversion_funcs(cls):
        validation_and_conversion_funcs = \
            cls._validation_and_conversion_funcs_.copy()

        return validation_and_conversion_funcs


    
    @classmethod
    def get_pre_serialization_funcs(cls):
        pre_serialization_funcs = \
            cls._pre_serialization_funcs_.copy()

        return pre_serialization_funcs


    
    @classmethod
    def get_de_pre_serialization_funcs(cls):
        de_pre_serialization_funcs = \
            cls._de_pre_serialization_funcs_.copy()

        return de_pre_serialization_funcs



def _pre_serialize_potential_slice_subset_ids(potential_slice_subset_ids):
    obj_to_pre_serialize = potential_slice_subset_ids
    serializable_rep = potential_slice_subset_ids.pre_serialize()

    return serializable_rep



def _de_pre_serialize_potential_slice_subset_ids(serializable_rep):
    sample_specification = \
        PotentialSliceSubsetIDs.de_pre_serialize(serializable_rep)

    return sample_specification



class SMatrixSubsetIDs(fancytypes.PreSerializableAndUpdatable):
    r"""A parameter set specifying a set of files, where each file stores the
    pre-calculated :math:`S`-matrices for a frozen phonon configuration subset.

    See the documentation for the subpackage :mod:`prismatique.stem` for a
    discussion on :math:`S`-matrices, and the class
    :class:`prismatique.thermal.Params` for a discussion on frozen phonon
    configuration subsets.

    HDF5 files storing pre-calculated :math:`S`-matrices can be generated using
    the functions :func:`prismatique.sample.generate_S_matrices`, and
    :func:`prismatique.stem.sim.run`.

    Parameters
    ----------
    filenames : `array_like` (`str`, ndim=1)
        For each integer ``i`` satisfying ``0<=i<len(filenames)``,
        ``filenames[i]`` specifies the filename of the pre-calculated
        :math:`S`-matrices for the ``i`` th frozen phonon configuration
        subset. The pre-calculated objects for each frozen phonon configuration
        subset are assumed to correspond to the same sample. In particular, each
        pre-calculated :math:`S`-matrix subset must have been calculated using
        the same sample supercell dimensions :math:`\left(\Delta X, \Delta Y,
        \Delta Z\right)`, potential slice pixel sizes :math:`\left(\Delta x,
        \Delta y\right)` in the :math:`x`- and :math:`y`-directions,
        interpolation factors :math:`\left(f_x, f_y\right)`, mean beam energy
        :math:`E`, and convergence semiangle :math:`\alpha_{\max}`. See the
        documentation for the class :class:`prismatique.discretization.Params`
        for discussions on :math:`\left(\Delta X, \Delta Y, \Delta Z\right)`,
        :math:`\left(\Delta x, \Delta y\right)`, and :math:`\left(f_x,
        f_y\right)`. See the documentation for the class
        :class:`embeam.stem.probe.ModelParams` for discussions on the beam
        energy and the convergence semiangle. Lastly, note that the total number
        of frozen phonon configuration subsets is ``len(filenames)``, and must
        be greater than or equal to 1.
    max_num_frozen_phonon_configs_per_subset : `int`, optional
        The maximum number of frozen phonon configurations to be loaded, per
        frozen phonon configuration subset. Due to a bug in ``prismatic``, the
        minimum possible value of ``max_num_frozen_phonon_configs_per_subset``
        is 2.
    skip_validation_and_conversion : `bool`, optional
        Let ``validation_and_conversion_funcs`` and ``core_attrs`` denote the
        attributes :attr:`~fancytypes.Checkable.validation_and_conversion_funcs`
        and :attr:`~fancytypes.Checkable.core_attrs` respectively, both of which
        being `dict` objects.

        Let ``params_to_be_mapped_to_core_attrs`` denote the `dict`
        representation of the constructor parameters excluding the parameter
        ``skip_validation_and_conversion``, where each `dict` key ``key`` is a
        different constructor parameter name, excluding the name
        ``"skip_validation_and_conversion"``, and
        ``params_to_be_mapped_to_core_attrs[key]`` would yield the value of the
        constructor parameter with the name given by ``key``.

        If ``skip_validation_and_conversion`` is set to ``False``, then for each
        key ``key`` in ``params_to_be_mapped_to_core_attrs``,
        ``core_attrs[key]`` is set to ``validation_and_conversion_funcs[key]
        (params_to_be_mapped_to_core_attrs)``.

        Otherwise, if ``skip_validation_and_conversion`` is set to ``True``,
        then ``core_attrs`` is set to
        ``params_to_be_mapped_to_core_attrs.copy()``. This option is desired
        primarily when the user wants to avoid potentially expensive deep copies
        and/or conversions of the `dict` values of
        ``params_to_be_mapped_to_core_attrs``, as it is guaranteed that no
        copies or conversions are made in this case.

    """
    ctor_param_names = ("filenames",
                        "max_num_frozen_phonon_configs_per_subset")
    kwargs = {"namespace_as_dict": globals(),
              "ctor_param_names": ctor_param_names}
    
    _validation_and_conversion_funcs_ = \
        fancytypes.return_validation_and_conversion_funcs(**kwargs)
    _pre_serialization_funcs_ = \
        fancytypes.return_pre_serialization_funcs(**kwargs)
    _de_pre_serialization_funcs_ = \
        fancytypes.return_de_pre_serialization_funcs(**kwargs)

    del ctor_param_names, kwargs

    
    
    def __init__(self,
                 filenames,
                 max_num_frozen_phonon_configs_per_subset=\
                 _default_max_num_frozen_phonon_configs_per_subset,
                 skip_validation_and_conversion=\
                 _default_skip_validation_and_conversion):
        ctor_params = {key: val
                       for key, val in locals().items()
                       if (key not in ("self", "__class__"))}
        kwargs = ctor_params
        kwargs["skip_cls_tests"] = True
        fancytypes.PreSerializableAndUpdatable.__init__(self, **kwargs)

        return None



    @classmethod
    def get_validation_and_conversion_funcs(cls):
        validation_and_conversion_funcs = \
            cls._validation_and_conversion_funcs_.copy()

        return validation_and_conversion_funcs


    
    @classmethod
    def get_pre_serialization_funcs(cls):
        pre_serialization_funcs = \
            cls._pre_serialization_funcs_.copy()

        return pre_serialization_funcs


    
    @classmethod
    def get_de_pre_serialization_funcs(cls):
        de_pre_serialization_funcs = \
            cls._de_pre_serialization_funcs_.copy()

        return de_pre_serialization_funcs



def _pre_serialize_S_matrix_subset_ids(S_matrix_subset_ids):
    obj_to_pre_serialize = S_matrix_subset_ids
    serializable_rep = S_matrix_subset_ids.pre_serialize()

    return serializable_rep



def _de_pre_serialize_S_matrix_subset_ids(serializable_rep):
    sample_specification = SMatrixSubsetIDs.de_pre_serialize(serializable_rep)

    return sample_specification



def _check_and_convert_potential_slice_subset_filenames(params):
    obj_name = "potential_slice_subset_filenames"

    kwargs = \
        {"obj": params[obj_name], "obj_name": obj_name}
    potential_slice_subset_filenames = \
        czekitout.convert.to_tuple_of_strs(**kwargs)

    current_func_name = "_check_and_convert_potential_slice_subset_filenames"

    if len(potential_slice_subset_filenames) == 0:
        err_msg = globals()[current_func_name+"_err_msg_1"]
        raise ValueError(err_msg)
    
    return potential_slice_subset_filenames



def _pre_serialize_potential_slice_subset_filenames(
        potential_slice_subset_filenames):
    obj_to_pre_serialize = potential_slice_subset_filenames
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_potential_slice_subset_filenames(serializable_rep):
    potential_slice_subset_filenames = serializable_rep

    return potential_slice_subset_filenames



def _check_and_convert_S_matrix_subset_filenames(params):
    obj_name = "S_matrix_subset_filenames"

    kwargs = \
        {"obj": params[obj_name], "obj_name": obj_name}
    S_matrix_subset_filenames = \
        czekitout.convert.to_tuple_of_strs(**kwargs)

    current_func_name = "_check_and_convert_S_matrix_subset_filenames"

    if len(S_matrix_subset_filenames) == 0:
        err_msg = globals()[current_func_name+"_err_msg_1"]
        raise ValueError(err_msg)
    
    return S_matrix_subset_filenames



def _pre_serialize_S_matrix_subset_filenames(S_matrix_subset_filenames):
    obj_to_pre_serialize = S_matrix_subset_filenames
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_S_matrix_subset_filenames(serializable_rep):
    S_matrix_subset_filenames = serializable_rep

    return S_matrix_subset_filenames



class PotentialSliceAndSMatrixSubsetIDs(fancytypes.PreSerializableAndUpdatable):
    r"""A parameter set specifying a set of files, where each file stores either
    the pre-calculated potential slices or :math:`S`-matrices for a frozen 
    phonon configuration subset.

    See the documentation for the subpackage :mod:`prismatique.stem` for a
    discussion on :math:`S`-matrices, the class
    :class:`prismatique.thermal.Params` for a discussion on frozen phonon
    configuration subsets, and the class
    :class:`prismatique.discretization.Params` for a discussion on potential
    slices.

    HDF5 files storing pre-calculated potential slices can be generated using
    the functions :func:`prismatique.sample.generate_potential_slices`,
    :func:`prismatique.stem.sim.run`, and :func:`prismatique.hrtem.sim.run`.

    HDF5 files storing pre-calculated :math:`S`-matrices can be generated using
    the functions :func:`prismatique.sample.generate_S_matrices`, and
    :func:`prismatique.stem.sim.run`.

    Parameters
    ----------
    potential_slice_subset_filenames : `array_like` (`str`, ndim=1)
        For each integer ``i`` satisfying
        ``0<=i<len(potential_slice_subset_filenames)``,
        ``potential_slice_subset_filenames[i]`` specifies the filename of the
        pre-calculated potential slices for the ``i`` th frozen phonon
        configuration subset. The pre-calculated objects for each frozen phonon
        configuration subset are assumed to correspond to the same sample. In
        particular, each pre-calculated potential slice subset must have been
        calculated using the same sample supercell dimensions
        :math:`\left(\Delta X, \Delta Y, \Delta Z\right)`, slice thickness
        :math:`\delta z`, and potential slice pixel sizes :math:`\left(\Delta x,
        \Delta y\right)` in the :math:`x`- and :math:`y`-directions. See the
        documentation for the class :class:`prismatique.discretization.Params`
        for discussions on the aforementioned quantities. Note that the total
        number of frozen phonon configurations is
        ``len(potential_slice_subset_filenames)+len(S_matrix_subset_filenames)``.
        Moreover, ``len(potential_slice_subset_filenames)`` must be greater than
        or equal to 1.
    S_matrix_subset_filenames : `array_like` (`str`, ndim=1)
        For each integer ``i`` satisfying
        ``0<=i<len(S_matrix_subset_filenames)``,
        ``S_matrix_subset_filenames[i]`` specifies the filename of the
        pre-calculated :math:`S`-matrices for the ``N+i`` th frozen phonon
        configuration subset, where
        ``N==len(potential_slice_subset_filenames)``. The pre-calculated objects
        for each frozen phonon configuration subset are assumed to correspond to
        the same sample. In particular, each pre-calculated :math:`S`-matrix
        subset must have been calculated using the same sample supercell
        dimensions :math:`\left(\Delta X, \Delta Y, \Delta Z\right)`, potential
        slice pixel sizes :math:`\left(\Delta x, \Delta y\right)` in the
        :math:`x`- and :math:`y`-directions, interpolation factors
        :math:`\left(f_x, f_y\right)`, mean beam energy :math:`E`, and
        convergence semiangle :math:`\alpha_{\max}`. See the documentation for
        the class :class:`prismatique.discretization.Params` for discussions on
        :math:`\left(\Delta X, \Delta Y, \Delta Z\right)`, :math:`\left(\Delta
        x, \Delta y\right)`, and :math:`\left(f_x, f_y\right)`. See the
        documentation for the class :class:`embeam.stem.probe.ModelParams` for
        discussions on the beam energy and the convergence semiangle. Note that
        the total number of frozen phonon configuration subsets is
        ``len(potential_slice_subset_filenames)+len(S_matrix_subset_filenames)``.
        Moreover, ``len(S_matrix_subset_filenames)`` must be greater than or
        equal to 1.
    max_num_frozen_phonon_configs_per_subset : `int`, optional
        The maximum number of frozen phonon configurations to be loaded, per
        frozen phonon configuration subset. Due to a bug in ``prismatic``, the
        minimum possible value of ``max_num_frozen_phonon_configs_per_subset``
        is 2.
    skip_validation_and_conversion : `bool`, optional
        Let ``validation_and_conversion_funcs`` and ``core_attrs`` denote the
        attributes :attr:`~fancytypes.Checkable.validation_and_conversion_funcs`
        and :attr:`~fancytypes.Checkable.core_attrs` respectively, both of which
        being `dict` objects.

        Let ``params_to_be_mapped_to_core_attrs`` denote the `dict`
        representation of the constructor parameters excluding the parameter
        ``skip_validation_and_conversion``, where each `dict` key ``key`` is a
        different constructor parameter name, excluding the name
        ``"skip_validation_and_conversion"``, and
        ``params_to_be_mapped_to_core_attrs[key]`` would yield the value of the
        constructor parameter with the name given by ``key``.

        If ``skip_validation_and_conversion`` is set to ``False``, then for each
        key ``key`` in ``params_to_be_mapped_to_core_attrs``,
        ``core_attrs[key]`` is set to ``validation_and_conversion_funcs[key]
        (params_to_be_mapped_to_core_attrs)``.

        Otherwise, if ``skip_validation_and_conversion`` is set to ``True``,
        then ``core_attrs`` is set to
        ``params_to_be_mapped_to_core_attrs.copy()``. This option is desired
        primarily when the user wants to avoid potentially expensive deep copies
        and/or conversions of the `dict` values of
        ``params_to_be_mapped_to_core_attrs``, as it is guaranteed that no
        copies or conversions are made in this case.

    """
    ctor_param_names = ("potential_slice_subset_filenames",
                        "S_matrix_subset_filenames",
                        "max_num_frozen_phonon_configs_per_subset")
    kwargs = {"namespace_as_dict": globals(),
              "ctor_param_names": ctor_param_names}
    
    _validation_and_conversion_funcs_ = \
        fancytypes.return_validation_and_conversion_funcs(**kwargs)
    _pre_serialization_funcs_ = \
        fancytypes.return_pre_serialization_funcs(**kwargs)
    _de_pre_serialization_funcs_ = \
        fancytypes.return_de_pre_serialization_funcs(**kwargs)

    del ctor_param_names, kwargs

    
    
    def __init__(self,
                 potential_slice_subset_filenames,
                 S_matrix_subset_filenames,
                 max_num_frozen_phonon_configs_per_subset=\
                 _default_max_num_frozen_phonon_configs_per_subset,
                 skip_validation_and_conversion=\
                 _default_skip_validation_and_conversion):
        ctor_params = {key: val
                       for key, val in locals().items()
                       if (key not in ("self", "__class__"))}
        kwargs = ctor_params
        kwargs["skip_cls_tests"] = True
        fancytypes.PreSerializableAndUpdatable.__init__(self, **kwargs)

        return None



    @classmethod
    def get_validation_and_conversion_funcs(cls):
        validation_and_conversion_funcs = \
            cls._validation_and_conversion_funcs_.copy()

        return validation_and_conversion_funcs


    
    @classmethod
    def get_pre_serialization_funcs(cls):
        pre_serialization_funcs = \
            cls._pre_serialization_funcs_.copy()

        return pre_serialization_funcs


    
    @classmethod
    def get_de_pre_serialization_funcs(cls):
        de_pre_serialization_funcs = \
            cls._de_pre_serialization_funcs_.copy()

        return de_pre_serialization_funcs



def _check_and_convert_sample_specification(params):
    obj_name = "sample_specification"
    obj = params[obj_name]

    accepted_types = params.get("accepted_types", None)
    if accepted_types is None:
        accepted_types = (prismatique.sample.ModelParams,
                          prismatique.sample.PotentialSliceSubsetIDs,
                          prismatique.sample.SMatrixSubsetIDs,
                          prismatique.sample.PotentialSliceAndSMatrixSubsetIDs)
    
    kwargs = {"obj": obj,
              "obj_name": obj_name,
              "accepted_types": accepted_types}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    obj_core_attrs = obj.get_core_attrs(deep_copy=False)

    if (("thermal_params" not in obj_core_attrs)
        and ("specimen_tilt" not in params)
        and ("lens_model_params" not in params)):
        _check_precalculated_subsets(sample_specification=obj)

    kwargs = obj_core_attrs
    sample_specification = obj.__class__(**kwargs)

    return sample_specification



def _check_precalculated_subsets(sample_specification):
    _check_that_hdf5_files_exist(sample_specification)
    _check_attr_subsets_in_hdf5_files(sample_specification)
    _check_datasets_in_hdf5_files(sample_specification)

    return None



def _check_that_hdf5_files_exist(sample_specification):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    core_attr_names_to_skip = ("interpolation_factors",
                               "max_num_frozen_phonon_configs_per_subset")

    current_func_name = "_check_that_hdf5_files_exist"
    
    for core_attr_name in sample_specification_core_attrs:
        if core_attr_name in core_attr_names_to_skip:
            continue
        filenames = sample_specification_core_attrs[core_attr_name]
        for idx, filename in enumerate(filenames):
            path_in_file = "/"
            obj_id = h5pywrappers.obj.ID(filename, path_in_file)
            try:
                _ = h5pywrappers.obj.load(obj_id, read_only=True)
            except:
                unformatted_err_msg = globals()[current_func_name+"_err_msg_1"]
                err_msg = unformatted_err_msg.format(core_attr_name,
                                                     idx,
                                                     filename)
                raise IOError(err_msg)
    
    return None



def _check_attr_subsets_in_hdf5_files(sample_specification):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    hdf5_attrs = dict()
    hdf5_attr_names = ("c", "t", "px", "py", "s", "fx", "fy")
    core_attr_names_to_skip = ("interpolation_factors",
                               "max_num_frozen_phonon_configs_per_subset")

    for hdf5_attr_name in hdf5_attr_names:
        for core_attr_name in sample_specification_core_attrs:
            if core_attr_name in core_attr_names_to_skip:
                continue
            filenames = sample_specification_core_attrs[core_attr_name]
            for idx, filename in enumerate(filenames):
                if ((hdf5_attr_name == "s")
                    and (not isinstance(sample_specification,
                                        PotentialSliceSubsetIDs))
                    and (core_attr_name != "potential_slice_subset_filenames")):
                    continue
                if ((hdf5_attr_name in ("fx", "fy"))
                    and (not isinstance(sample_specification,
                                        SMatrixSubsetIDs))
                    and (core_attr_name != "S_matrix_subset_filenames")):
                    continue

                kwargs = {"hdf5_attr_name": hdf5_attr_name,
                          "core_attr_name": core_attr_name,
                          "idx": idx,
                          "filename": filename}
                hdf5_attr = np.array(_load_hdf5_attr(**kwargs))

                kwargs = {"hdf5_attr": hdf5_attr,
                          "hdf5_attr_name": hdf5_attr_name,
                          "hdf5_attrs": hdf5_attrs,
                          "sample_specification": sample_specification,
                          "core_attr_name": core_attr_name}
                _check_hdf5_attr(**kwargs)

    return None



def _load_hdf5_attr(hdf5_attr_name, core_attr_name, idx, filename):
    path_in_file = ("4DSTEM_simulation/metadata/metadata_0/original"
                    "/simulation_parameters")
    obj_id = h5pywrappers.obj.ID(filename, path_in_file)
    attr_id = h5pywrappers.attr.ID(obj_id, hdf5_attr_name)

    current_func_name = "_load_hdf5_attr"

    try:
        try:
            hdf5_attr = h5pywrappers.attr.load(attr_id)
        except:
            end_of_err_msg = "."
            raise

        kwargs = {"obj": hdf5_attr, "obj_name": "hdf5_attr"}
        
        end_of_err_msg = " that stores a "
        if hdf5_attr_name == "c":
            end_of_err_msg += "3-element array of positive real numbers."
            hdf5_attr = czekitout.convert.to_tuple_of_positive_floats(**kwargs)
            if len(hdf5_attr) != 3:
                raise
        elif hdf5_attr_name == "t":
            end_of_err_msg += "3-element array of positive integers."
            hdf5_attr = czekitout.convert.to_tuple_of_positive_ints(**kwargs)
            if len(hdf5_attr) != 3:
                raise
        elif hdf5_attr_name in ("s", "px", "py"):
            end_of_err_msg += "positive real number."
            hdf5_attr = czekitout.convert.to_positive_float(**kwargs)
        else:
            end_of_err_msg += "positive integer."
            hdf5_attr = czekitout.convert.to_positive_int(**kwargs)
    except:
        unformatted_err_msg = globals()[current_func_name+"_err_msg_1"]
        err_msg = unformatted_err_msg.format(filename,
                                             core_attr_name,
                                             idx,
                                             hdf5_attr_name,
                                             end_of_err_msg)
        error_type = IOError if (end_of_err_msg == ".") else TypeError
        raise error_type(err_msg)

    return hdf5_attr



def _check_hdf5_attr(hdf5_attr,
                     hdf5_attr_name,
                     hdf5_attrs,
                     sample_specification,
                     core_attr_name):
    if hdf5_attr_name not in hdf5_attrs:
        hdf5_attrs[hdf5_attr_name] = hdf5_attr
    abs_diff = np.linalg.norm(hdf5_attrs[hdf5_attr_name]-hdf5_attr)
    rel_diff = abs_diff / np.linalg.norm(hdf5_attrs[hdf5_attr_name])

    tol = 1e-6

    current_func_name = "_check_hdf5_attr"
                
    if rel_diff >= tol:
        if (isinstance(sample_specification, PotentialSliceAndSMatrixSubsetIDs)
            and (hdf5_attr_name not in ("s", "fx", "fy"))):
            unformatted_err_msg = globals()[current_func_name+"_err_msg_1"]
            err_msg = unformatted_err_msg.format(hdf5_attr_name)
        else:
            unformatted_err_msg = globals()[current_func_name+"_err_msg_2"]
            err_msg = unformatted_err_msg.format(core_attr_name,
                                                 hdf5_attr_name,
                                                 core_attr_name)
                        
        raise ValueError(err_msg)

    return None



def _check_datasets_in_hdf5_files(sample_specification):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    core_attr_names_to_skip = ("interpolation_factors",
                               "max_num_frozen_phonon_configs_per_subset")

    current_func_name = "_check_datasets_in_hdf5_files"
    
    for core_attr_name in sample_specification_core_attrs:
        if core_attr_name in core_attr_names_to_skip:
            continue
        filenames = sample_specification_core_attrs[core_attr_name]
        
        for idx, filename in enumerate(filenames):
            if (("interpolation_factors" in sample_specification_core_attrs)
                or (core_attr_name == "potential_slice_subset_filenames")):
                num_valid_potential_slice_datasets_in_file = \
                    _num_valid_potential_slice_datasets_in_file(idx,
                                                                filename,
                                                                core_attr_name)
                
                if num_valid_potential_slice_datasets_in_file == 0:
                    unformatted_err_msg = \
                        globals()[current_func_name+"_err_msg_1"]
                    err_msg = \
                        unformatted_err_msg.format(filename,
                                                   core_attr_name,
                                                   idx)
                    
                    raise IOError(err_msg)
            else:
                num_valid_S_matrix_datasets_in_file = \
                    _num_valid_S_matrix_datasets_in_file(idx,
                                                         filename,
                                                         core_attr_name)

                if num_valid_S_matrix_datasets_in_file == 0:
                    unformatted_err_msg = \
                        globals()[current_func_name+"_err_msg_2"]
                    err_msg = \
                        unformatted_err_msg.format(filename,
                                                   core_attr_name,
                                                   idx)
                    
                    raise IOError(err_msg)

    return None



def _num_valid_potential_slice_datasets_in_file(idx, filename, core_attr_name):
    num_valid_potential_slice_datasets_in_file = 0
    
    ppotential_fp_indices = _ppotential_fp_indices(idx,
                                                   filename,
                                                   core_attr_name)
    num_ppotential_fp_indices = len(ppotential_fp_indices)

    sample_supercell_dims, sample_supercell_voxel_size = \
        _load_sample_supercell_dims_and_voxel_size_from_file(idx, filename)
    sample_supercell_dims = \
        np.array(sample_supercell_dims)
    sample_supercell_voxel_size = \
        np.array(sample_supercell_voxel_size)

    for ppotential_fp_idx in ppotential_fp_indices:
        path_in_file = ("4DSTEM_simulation/data/realslices"
                        "/ppotential_fp"+str(ppotential_fp_idx).rjust(4, "0")
                        +"/data")
        dataset_id = h5pywrappers.obj.ID(filename, path_in_file)
        
        try:
            dataset = h5pywrappers.dataset.load(dataset_id)
            dataset_shape = np.array(dataset.shape)
            dataset_dtype = dataset.dtype
            dataset.file.close()

            abs_diff = np.linalg.norm(dataset_shape*sample_supercell_voxel_size
                                      - sample_supercell_dims)
            rel_diff = abs_diff / np.linalg.norm(sample_supercell_dims)

            tol = 1e-6

            if ((num_valid_potential_slice_datasets_in_file
                 != ppotential_fp_idx)
                or (dataset_dtype != "float32")
                or (rel_diff >= tol)):
                break
            
        except:
            break

        num_valid_potential_slice_datasets_in_file += 1
        
    return num_valid_potential_slice_datasets_in_file



def _ppotential_fp_indices(idx, filename, core_attr_name):
    path_in_file = "4DSTEM_simulation/data/realslices"
    obj_id = h5pywrappers.obj.ID(filename, path_in_file)

    current_func_name = "_ppotential_fp_indices"
    
    try:
        obj = h5pywrappers.obj.load(obj_id, read_only=True)
        obj_keys = tuple(obj.keys())
        obj.file.close()
    except:
        unformatted_err_msg = globals()[current_func_name+"_err_msg_1"]
        err_msg = unformatted_err_msg.format(filename, core_attr_name, idx)
        raise IOError(err_msg)

    ppotential_fp_indices = tuple()
    for key in obj_keys:
        match = re.fullmatch(r"ppotential_fp[0-9]{4}", key)
        if match is not None:
            ppotential_fp_idx = int(match.string[13:])
            ppotential_fp_indices += (ppotential_fp_idx,)
    ppotential_fp_indices = sorted(ppotential_fp_indices)

    return ppotential_fp_indices



def _load_sample_supercell_dims_and_voxel_size_from_file(idx, filename):
    sample_supercell_dims = \
        _load_sample_supercell_dims_from_file(filename)
    sample_supercell_lateral_pixel_size = \
        _load_sample_supercell_lateral_pixel_size_from_file(filename)
    slice_thickness = \
        _load_hdf5_attr("s", "", idx, filename)  # 2nd arg not used here.

    sample_supercell_voxel_size = \
        (sample_supercell_lateral_pixel_size[0],
         sample_supercell_lateral_pixel_size[1],
         slice_thickness)
    
    return sample_supercell_dims, sample_supercell_voxel_size



def _load_sample_supercell_dims_from_file(filename):
    path_in_file = ("4DSTEM_simulation/metadata/metadata_0"
                    "/original/simulation_parameters")
    obj_id = h5pywrappers.obj.ID(filename, path_in_file)
    obj = h5pywrappers.obj.load(obj_id, read_only=True)

    sample_supercell_dims = tuple((obj.attrs["c"]*obj.attrs["t"]).tolist())
    
    obj.file.close()

    return sample_supercell_dims



def _load_sample_supercell_lateral_pixel_size_from_file(filename):
    path_in_file = ("4DSTEM_simulation/metadata/metadata_0"
                    "/original/simulation_parameters")
    obj_id = h5pywrappers.obj.ID(filename, path_in_file)
    obj = h5pywrappers.obj.load(obj_id, read_only=True)

    sample_supercell_lateral_pixel_size = (float(obj.attrs["px"]),
                                           float(obj.attrs["py"]))
    
    obj.file.close()

    return sample_supercell_lateral_pixel_size



def _num_valid_S_matrix_datasets_in_file(idx, filename, core_attr_name):
    num_valid_S_matrix_datasets_in_file = 0
    
    smatrix_fp_indices = _smatrix_fp_indices(idx, filename, core_attr_name)
    num_smatrix_fp_indices = len(smatrix_fp_indices)

    interpolation_factors = \
        np.array(_load_interpolation_factors_from_file(idx, filename))
    sample_supercell_xy_dims_in_pixels = \
        np.array(_load_sample_supercell_xy_dims_in_pixels_from_file(filename))

    for smatrix_fp_idx in smatrix_fp_indices:
        path_in_file = ("4DSTEM_simulation/data/realslices"
                        "/smatrix_fp"+str(smatrix_fp_idx).rjust(4, "0")+"/data")
        dataset_id = h5pywrappers.obj.ID(filename, path_in_file)
        
        try:
            dataset = h5pywrappers.dataset.load(dataset_id)
            dataset_shape = np.array(dataset.shape)
            dataset_dtype = dataset.dtype
            dataset.file.close()

            if num_valid_S_matrix_datasets_in_file == 0:
                num_S_matrix_k_xy_vectors = dataset_shape[2]

            abs_diff_1 = \
                np.linalg.norm(2*dataset_shape[0:2]*interpolation_factors
                               - sample_supercell_xy_dims_in_pixels)
            rel_diff_1 = \
                abs_diff_1 / np.linalg.norm(sample_supercell_xy_dims_in_pixels)

            abs_diff_2 = np.linalg.norm(num_S_matrix_k_xy_vectors
                                        - dataset_shape[2])
            rel_diff_2 = abs_diff_2 / np.linalg.norm(num_S_matrix_k_xy_vectors)

            tol = 1e-6

            if ((num_valid_S_matrix_datasets_in_file != smatrix_fp_idx)
                or (dataset_dtype != "complex64")
                or (rel_diff_1 >= tol)
                or (rel_diff_2 >= tol)
                or (len(dataset_shape) != 3)):
                break
            
        except:
            break

        num_valid_S_matrix_datasets_in_file += 1
        
    return num_valid_S_matrix_datasets_in_file



def _load_sample_supercell_xy_dims_in_pixels_from_file(filename):
    sample_supercell_dims = \
        np.array(_load_sample_supercell_dims_from_file(filename))
    sample_supercell_xy_dims = \
        sample_supercell_dims[0:2]
    
    sample_supercell_lateral_pixel_size = \
        np.array(_load_sample_supercell_lateral_pixel_size_from_file(filename))
    
    sample_supercell_xy_dims_in_pixels = \
        sample_supercell_xy_dims / sample_supercell_lateral_pixel_size
    sample_supercell_xy_dims_in_pixels = \
        np.round(sample_supercell_xy_dims_in_pixels).astype(int)
    sample_supercell_xy_dims_in_pixels = \
        tuple(sample_supercell_xy_dims_in_pixels.tolist())

    return sample_supercell_xy_dims_in_pixels



def _load_interpolation_factors_from_file(idx, filename):
    fx = _load_hdf5_attr("fx", "", idx, filename)  # 2nd arg not used here.
    fy = _load_hdf5_attr("fy", "", idx, filename)  # 2nd arg not used here.
    interpolation_factors = (fx, fy)

    return interpolation_factors



def _smatrix_fp_indices(idx, filename, core_attr_name):
    path_in_file = "4DSTEM_simulation/data/realslices"
    obj_id = h5pywrappers.obj.ID(filename, path_in_file)

    current_func_name = "_smatrix_fp_indices"
    
    try:
        obj = h5pywrappers.obj.load(obj_id, read_only=True)
        obj_keys = tuple(obj.keys())
        obj.file.close()
    except:
        unformatted_err_msg = globals()[current_func_name+"_err_msg_1"]
        err_msg = unformatted_err_msg.format(filename, core_attr_name, idx)
        raise IOError(err_msg)

    smatrix_fp_indices = tuple()
    for key in obj_keys:
        match = re.fullmatch(r"smatrix_fp[0-9]{4}", key)
        if match is not None:
            smatrix_fp_idx = int(match.string[13:])
            smatrix_fp_indices += (smatrix_fp_idx,)
    smatrix_fp_indices = sorted(smatrix_fp_indices)

    return smatrix_fp_indices



def _pre_serialize_sample_specification(sample_specification):
    obj_to_pre_serialize = sample_specification

    obj_to_pre_serialize_core_attrs = \
        obj_to_pre_serialize.get_core_attrs(deep_copy=False)

    if "thermal_params" in obj_to_pre_serialize_core_attrs:
        func_name = "_pre_serialize_sample_model_params"
    elif "interpolation_factors" in obj_to_pre_serialize_core_attrs:
        func_name = "_pre_serialize_potential_slice_subset_ids"
    else:
        func_name = "_pre_serialize_S_matrix_subset_ids"
    
    func_alias = globals()[func_name]
    serializable_rep = func_alias(sample_specification)

    return serializable_rep



def _de_pre_serialize_sample_specification(serializable_rep):
    if "thermal_params" in serializable_rep:
        func_name = "_de_pre_serialize_sample_model_params"
    elif "interpolation_factors" in serializable_rep:
        func_name = "_de_pre_serialize_potential_slice_subset_ids"
    else:
        func_name = "_de_pre_serialize_S_matrix_subset_ids"
        
    func_alias = globals()[func_name]
    sample_specification = func_alias(serializable_rep)

    return sample_specification



def _check_and_convert_probe_model_params(params):
    obj_name = "probe_model_params"

    module_alias = embeam.stem.probe.kspace
    cls_alias = module_alias.Intensity
    func_alias = cls_alias.get_validation_and_conversion_funcs()[obj_name]
    probe_model_params = func_alias(params)

    if "specimen_tilt" not in params:
        sample_specification = \
            _check_and_convert_sample_specification(params)
        sample_specification_core_attrs = \
            sample_specification.get_core_attrs(deep_copy=False)

        if "thermal_params" not in sample_specification_core_attrs:
            kwargs = {"sample_specification": sample_specification,
                      "probe_model_params": probe_model_params}
            _check_that_probe_model_is_compatible_with_S_matrices(**kwargs)

    return probe_model_params



def _check_that_probe_model_is_compatible_with_S_matrices(sample_specification,
                                                          probe_model_params):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    if "interpolation_factors" in sample_specification_core_attrs:
        return None

    core_attr_names_to_not_skip = ("filenames", "S_matrix_subset_filenames")

    core_attr_name = ("filenames"
                      if ("filenames" in sample_specification_core_attrs)
                      else "S_matrix_subset_filenames")
    filenames = sample_specification_core_attrs[core_attr_name]

    kwargs = {"sample_specification": sample_specification,
              "probe_model_params": probe_model_params}
    S_matrix_k_xy_vectors = _S_matrix_k_xy_vectors(**kwargs)
    num_S_matrix_k_xy_vectors = len(S_matrix_k_xy_vectors)

    current_func_name = "_check_that_probe_model_is_compatible_with_S_matrices"
    
    for idx, filename in enumerate(filenames):
        path_in_file = ("4DSTEM_simulation/data/realslices"
                        "/smatrix_fp0000/data")
        dataset_id = h5pywrappers.obj.ID(filename, path_in_file)
        dataset = h5pywrappers.dataset.load(dataset_id)
        dataset_shape = np.array(dataset.shape)
        dataset.file.close()

        if num_S_matrix_k_xy_vectors != dataset_shape[2]:
            unformatted_err_msg = globals()[current_func_name+"_err_msg_1"]
            err_msg = unformatted_err_msg.format(filename,
                                                 core_attr_name,
                                                 idx)
            raise IOError(err_msg)

    return None



def _S_matrix_k_xy_vectors(sample_specification, probe_model_params):
    probe_model_params_core_attrs = \
        probe_model_params.get_core_attrs(deep_copy=False)
    gun_model_params = \
        probe_model_params_core_attrs["gun_model_params"]

    gun_model_params_core_attrs = \
        gun_model_params.get_core_attrs(deep_copy=False)
    mean_beam_energy = \
        gun_model_params_core_attrs["mean_beam_energy"]
    
    electron_beam_wavelength = embeam.wavelength(mean_beam_energy)
    
    convergence_semiangle = \
        probe_model_params_core_attrs["convergence_semiangle"]

    max_scattering_angle = convergence_semiangle + 2.5  # In mrads.
    tilt_offset = (0, 0)
    tilt_window = (0, max_scattering_angle)
    for_a_hrtem_calculation = False

    angular_mesh, beam_mask = \
        _angular_mesh_and_beam_mask(sample_specification,
                                    mean_beam_energy,
                                    tilt_offset,
                                    tilt_window,
                                    for_a_hrtem_calculation)

    # The for-loop ordering has been reversed on purpose, i.e. we start with the
    # ``k_y_idx`` for loop rather than the ``k_x_idx``. This is done in order to
    # be consistent with ``prismatic``.
    k_xy_vectors = tuple()
    for k_x_idx in range(beam_mask.shape[0]):
        for k_y_idx in range(beam_mask.shape[1]):
            if beam_mask[k_x_idx][k_y_idx]:
                k_x_coord = \
                    angular_mesh[0][k_x_idx][k_y_idx]/electron_beam_wavelength
                k_y_coord = \
                    angular_mesh[1][k_x_idx][k_y_idx]/electron_beam_wavelength
                k_xy_vectors += ((k_x_coord, k_y_coord),)
    
    return k_xy_vectors



def _angular_mesh_and_beam_mask(sample_specification,
                                mean_beam_energy,
                                tilt_offset,
                                tilt_window,
                                for_a_hrtem_calculation):
    # This code is essentially following the prismatic code.
    angular_mesh = _angular_mesh(sample_specification, mean_beam_energy)
    beam_idx_mesh = _beam_idx_mesh(sample_specification)
    k_mask = _k_mask(sample_specification)
    
    f_x, f_y = \
        _interpolation_factors_from_sample_specification(sample_specification)

    temp_mask = ((beam_idx_mesh[0] % f_x == 0)
                 * (beam_idx_mesh[1] % f_y == 0)
                 * (k_mask == 1))

    if len(tilt_window) == 4:
        rel_x_tilts, rel_y_tilts = \
            _rel_x_and_rel_y_tilts(angular_mesh, tilt_offset)
        min_x_tilt, max_x_tilt, min_y_tilt, max_y_tilt = \
            np.array(tilt_window) / 1000  # In rads.
        beam_mask = \
            (temp_mask
             * (rel_x_tilts <= max_x_tilt)
             * (rel_x_tilts >= min_x_tilt)
             * (rel_y_tilts <= max_y_tilt)
             * (rel_y_tilts >= min_y_tilt))
    else:
        rel_tilts = _rel_tilts(angular_mesh, tilt_offset)
        min_r_tilt, max_r_tilt = np.array(tilt_window) / 1000  # In rads.
        if for_a_hrtem_calculation:
            beam_mask = (temp_mask
                         * (rel_tilts <= max_r_tilt)
                         * (rel_tilts >= min_r_tilt))
        else:
            beam_mask = (temp_mask
                         * (rel_tilts < max_r_tilt)
                         * (rel_tilts >= min_r_tilt))

    return angular_mesh, beam_mask



def _angular_mesh(sample_specification, mean_beam_energy):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    if "atomic_coords_filename" in sample_specification_core_attrs:
        sample_supercell_xy_dims_in_pixels = \
            _supercell_xy_dims_in_pixels(sample_specification)
        sample_supercell_lateral_pixel_size = \
            _supercell_lateral_pixel_size(sample_specification)
    else:
        core_attr_name = ("filenames"
                          if ("filenames" in sample_specification_core_attrs)
                          else "potential_slice_subset_filenames")
        filename = sample_specification_core_attrs[core_attr_name][0]
        
        sample_supercell_xy_dims_in_pixels = \
            _load_sample_supercell_xy_dims_in_pixels_from_file(filename)
        sample_supercell_lateral_pixel_size = \
            _load_sample_supercell_lateral_pixel_size_from_file(filename)
        
    electron_beam_wavelength = embeam.wavelength(mean_beam_energy)

    k_x = _FFT_1D_freqs(sample_supercell_xy_dims_in_pixels[0],
                        sample_supercell_lateral_pixel_size[0])
    k_y = _FFT_1D_freqs(sample_supercell_xy_dims_in_pixels[1],
                        sample_supercell_lateral_pixel_size[1])
    k_mesh = np.meshgrid(k_x, k_y, indexing="ij")

    angular_mesh = (k_mesh[0] * electron_beam_wavelength,
                    k_mesh[1] * electron_beam_wavelength)

    return angular_mesh



def _FFT_1D_freqs(num_samples, sample_spacing):
    freqs = np.zeros([num_samples])
    N_c = num_samples // 2

    df = 1 / ((num_samples-(num_samples%2)) * sample_spacing)
    for idx in range(num_samples):
        freqs[(N_c + idx + (num_samples%2)) % num_samples] = (idx - N_c) * df

    return freqs



def _beam_idx_mesh(sample_specification):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    if "atomic_coords_filename" in sample_specification_core_attrs:
        sample_supercell_xy_dims_in_pixels = \
            _supercell_xy_dims_in_pixels(sample_specification)
    else:
        core_attr_name = ("filenames"
                          if ("filenames" in sample_specification_core_attrs)
                          else "potential_slice_subset_filenames")
        filename = sample_specification_core_attrs[core_attr_name][0]
        
        sample_supercell_xy_dims_in_pixels = \
            _load_sample_supercell_xy_dims_in_pixels_from_file(filename)

    xv = _FFT_1D_freqs(sample_supercell_xy_dims_in_pixels[0],
                       1/sample_supercell_xy_dims_in_pixels[0])
    yv = _FFT_1D_freqs(sample_supercell_xy_dims_in_pixels[1],
                       1/sample_supercell_xy_dims_in_pixels[1])
    beam_idx_mesh = np.round(np.meshgrid(xv, yv, indexing="ij"))
    beam_idx_mesh = beam_idx_mesh.astype(int)

    return beam_idx_mesh



def _k_mask(sample_specification):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    if "atomic_coords_filename" in sample_specification_core_attrs:
        sample_supercell_xy_dims_in_pixels = \
            _supercell_xy_dims_in_pixels(sample_specification)
    else:
        core_attr_name = ("filenames"
                          if ("filenames" in sample_specification_core_attrs)
                          else "potential_slice_subset_filenames")
        filename = sample_specification_core_attrs[core_attr_name][0]
        
        sample_supercell_xy_dims_in_pixels = \
            _load_sample_supercell_xy_dims_in_pixels_from_file(filename)
        
    k_mask = np.zeros(sample_supercell_xy_dims_in_pixels)

    k_x_idx_offset = k_mask.shape[0] // 4
    k_y_idx_offset = k_mask.shape[1] // 4

    N_c_x = k_mask.shape[0] // 2
    N_c_y = k_mask.shape[1] // 2

    for k_x_idx in range(N_c_x):
        adjust_k_x_idx = (k_x_idx - k_x_idx_offset) % k_mask.shape[0]
        for k_y_idx in range(N_c_y):
            adjust_k_y_idx = (k_y_idx - k_y_idx_offset) % k_mask.shape[1]
            k_mask[adjust_k_x_idx, adjust_k_y_idx] = 1

    return k_mask



def _interpolation_factors_from_sample_specification(sample_specification):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    if "atomic_coords_filename" in sample_specification_core_attrs:
        discretization_params = \
            sample_specification_core_attrs["discretization_params"]
        discretization_params_core_attrs = \
            discretization_params.get_core_attrs(deep_copy=False)
        interpolation_factors = \
            discretization_params_core_attrs["interpolation_factors"]
    else:
        if "interpolation_factors" in sample_specification_core_attrs:
            interpolation_factors = \
                sample_specification_core_attrs["interpolation_factors"]
        else:
            core_attr_name = \
                ("filenames"
                 if ("filenames" in sample_specification_core_attrs)
                 else "S_matrix_subset_filenames")
            filename = \
                sample_specification_core_attrs[core_attr_name][0]
            interpolation_factors = \
                _load_interpolation_factors_from_file(0, filename)

    return interpolation_factors



def _rel_x_and_rel_y_tilts(angular_mesh, tilt_offset):
    x_tilts = angular_mesh[0]  # In rads.
    y_tilts = angular_mesh[1]  # In rads.
    rel_x_tilts = np.abs(x_tilts - tilt_offset[0] / 1000)  # In rads.
    rel_y_tilts = np.abs(y_tilts - tilt_offset[1] / 1000)  # In rads.

    return rel_x_tilts, rel_y_tilts



def _rel_tilts(angular_mesh, tilt_offset):
    rel_x_tilts, rel_y_tilts = _rel_x_and_rel_y_tilts(angular_mesh, tilt_offset)
    rel_tilts = np.sqrt(rel_x_tilts*rel_x_tilts + rel_y_tilts*rel_y_tilts)

    return rel_tilts



def _pre_serialize_probe_model_params(probe_model_params):
    obj_to_pre_serialize = probe_model_params

    obj_name = "probe_model_params"

    module_alias = embeam.stem.probe.kspace
    cls_alias = module_alias.Intensity
    func_alias = cls_alias.get_pre_serialization_funcs()[obj_name]
    serializable_rep = func_alias(obj_to_pre_serialize)
    
    return serializable_rep



def _de_pre_serialize_probe_model_params(serializable_rep):
    obj_name = "probe_model_params"

    module_alias = embeam.stem.probe.kspace
    cls_alias = module_alias.Intensity
    func_alias = cls_alias.get_de_pre_serialization_funcs()[obj_name]
    probe_model_params = func_alias(serializable_rep)

    return probe_model_params



def _check_and_convert_skip_validation_and_conversion(params):
    obj_name = "skip_validation_and_conversion"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    skip_validation_and_conversion = czekitout.convert.to_bool(**kwargs)

    return skip_validation_and_conversion



def unit_cell_dims(sample_specification,
                   skip_validation_and_conversion=\
                   _default_skip_validation_and_conversion):
    r"""Calculate the :math:`x`-, :math:`y`-, and :math:`z`-dimensions of the 
    sample's unit cell.

    See the documentation for class :class:`prismatique.discretization.Params`
    for a discussion on sample unit cells.

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
    sample_unit_cell_dims : `array_like` (`int`, shape=(``2``,))
        ``sample_unit_cell_dims[0]``, ``sample_unit_cell_dims[1]``, and
        ``sample_unit_cell_dims[2]`` are respectively the :math:`x`-,
        :math:`y`-, and :math:`z`-dimensions of the sample's unit cell in units
        of angstroms, and are denoted respectively as :math:`\Delta X`,
        :math:`\Delta Y`, and :math:`\Delta Z` in the documentation for the
        class :class:`prismatique.discretization.Params`.

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
    sample_unit_cell_dims = _unit_cell_dims(**kwargs)
    
    return sample_unit_cell_dims



def _unit_cell_dims(sample_specification):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    if "atomic_coords_filename" in sample_specification_core_attrs:
        atomic_coords_filename = \
            sample_specification_core_attrs["atomic_coords_filename"]
        
        check_atomic_coords_file_format(atomic_coords_filename)
        with open(atomic_coords_filename, "r") as file_obj:
            throw_away_line = file_obj.readline()
            strings = file_obj.readline().split()
        sample_unit_cell_dims = tuple(float(s) for s in strings)
    else:
        core_attr_name = ("filenames"
                          if ("filenames" in sample_specification_core_attrs)
                          else "potential_slice_subset_filenames")
        filename = sample_specification_core_attrs[core_attr_name][0]
        
        sample_unit_cell_dims = \
            _load_sample_unit_cell_dims_from_file(filename)

    return sample_unit_cell_dims



def _load_sample_unit_cell_dims_from_file(filename):
    path_in_file = ("4DSTEM_simulation/metadata/metadata_0"
                    "/original/simulation_parameters")
    obj_id = h5pywrappers.obj.ID(filename, path_in_file)
    obj = h5pywrappers.obj.load(obj_id, read_only=True)

    sample_unit_cell_dims = tuple(obj.attrs["c"].tolist())
    
    obj.file.close()

    return sample_unit_cell_dims



def supercell_dims(sample_specification,
                   skip_validation_and_conversion=\
                   _default_skip_validation_and_conversion):
    r"""Calculate the :math:`x`-, :math:`y`-, and :math:`z`-dimensions of the 
    sample's supercell.

    See the documentation for class :class:`prismatique.discretization.Params`
    for a discussion on sample supercells.

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
    sample_supercell_dims : `array_like` (`int`, shape=(``3``,))
        ``sample_supercell_dims[0]``, ``sample_supercell_dims[1]``, and
        ``sample_supercell_dims[2]`` are respectively the :math:`x`-,
        :math:`y`-, and :math:`z`-dimensions of the sample's supercell in units
        of angstroms, and are denoted respectively as :math:`\Delta X`,
        :math:`\Delta Y`, and :math:`\Delta Z` in the documentation for the
        class :class:`prismatique.discretization.Params`.

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
    sample_supercell_dims = _supercell_dims(**kwargs)
    
    return sample_supercell_dims



def _supercell_dims(sample_specification):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    if "atomic_coords_filename" in sample_specification_core_attrs:
        atomic_coords_filename = \
            sample_specification_core_attrs["atomic_coords_filename"]
        
        check_atomic_coords_file_format(atomic_coords_filename)
        with open(atomic_coords_filename, "r") as file_obj:
            throw_away_line = file_obj.readline()
            strings = file_obj.readline().split()
        sample_unit_cell_dims = tuple(float(s) for s in strings)

        unit_cell_tiling = sample_specification_core_attrs["unit_cell_tiling"]

        sample_supercell_dims = (sample_unit_cell_dims[0]*unit_cell_tiling[0],
                                 sample_unit_cell_dims[1]*unit_cell_tiling[1],
                                 sample_unit_cell_dims[2]*unit_cell_tiling[2])
    else:
        core_attr_name = \
            ("filenames"
             if ("filenames" in sample_specification_core_attrs)
             else "potential_slice_subset_filenames")
        filename = \
            sample_specification_core_attrs[core_attr_name][0]
        
        sample_supercell_dims = \
            _load_sample_supercell_dims_from_file(filename)

    return sample_supercell_dims



def supercell_xy_dims_in_pixels(sample_specification,
                                skip_validation_and_conversion=\
                                _default_skip_validation_and_conversion):
    r"""Calculate the :math:`x`- and :math:`y`-dimensions of the sample's 
    supercell in units of pixels.

    See the documentation for class :class:`prismatique.discretization.Params`
    for a discussion on sample supercells and how they are discretized in
    real-space.

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

          * discretization_params

            + sample_supercell_reduced_xy_dims_in_pixels

            + interpolation_factors

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
    sample_supercell_xy_dims_in_pixels : `array_like` (`int`, shape=(``2``,))
        ``sample_supercell_xy_dims_in_pixels[0]`` and
        ``sample_supercell_xy_dims_in_pixels[1]`` are respectively the
        :math:`x`- and :math:`y`-dimensions of the sample's supercell in units
        of pixels, and are denoted respectively as :math:`N_x` and :math:`N_y`
        in the documentation for the class
        :class:`prismatique.discretization.Params`.

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
    sample_supercell_xy_dims_in_pixels = _supercell_xy_dims_in_pixels(**kwargs)

    return sample_supercell_xy_dims_in_pixels



def _supercell_xy_dims_in_pixels(sample_specification):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    if "atomic_coords_filename" in sample_specification_core_attrs:
        discretization_params = \
            sample_specification_core_attrs["discretization_params"]
        discretization_params_core_attrs = \
            discretization_params.get_core_attrs(deep_copy=False)

        core_attr_name = "sample_supercell_reduced_xy_dims_in_pixels"
        tilde_N_x, tilde_N_y = discretization_params_core_attrs[core_attr_name]

        core_attr_name = "interpolation_factors"
        f_x, f_y = discretization_params_core_attrs[core_attr_name]

        sample_supercell_xy_dims_in_pixels = (4*f_x*tilde_N_x, 4*f_y*tilde_N_y)
    else:
        core_attr_name = \
            ("filenames"
             if ("filenames" in sample_specification_core_attrs)
             else "potential_slice_subset_filenames")
        filename = \
            sample_specification_core_attrs[core_attr_name][0]
        
        sample_supercell_xy_dims_in_pixels = \
            _load_sample_supercell_xy_dims_in_pixels_from_file(filename)

    return sample_supercell_xy_dims_in_pixels



def supercell_lateral_pixel_size(sample_specification,
                                 skip_validation_and_conversion=\
                                 _default_skip_validation_and_conversion):
    r"""Calculate the :math:`x`- and :math:`y` dimensions of the real-space 
    pixels used to discretize the sample's supercell.

    See the documentation for class :class:`prismatique.discretization.Params`
    for a discussion on sample supercells and how they are discretized in
    real-space.

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

          * discretization_params

            + sample_supercell_reduced_xy_dims_in_pixels

            + interpolation_factors

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
    sample_supercell_lateral_pixel_size : `array_like` (`float`, shape=(``2``,))
        ``sample_supercell_xy_dims_in_pixels[0]`` and
        ``sample_supercell_xy_dims_in_pixels[1]`` are respectively the
        :math:`x`- and :math:`y` dimensions of the real-space pixels used to
        discretize the sample's supercell, in units of angstroms, and are
        denoted respectively as :math:`\Delta x` and :math:`\Delta y` in the
        documentation for the class :class:`prismatique.discretization.Params`.

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

    func_alias = _supercell_lateral_pixel_size
    kwargs = params
    del kwargs["skip_validation_and_conversion"]    
    sample_supercell_lateral_pixel_size = func_alias(**kwargs)

    return sample_supercell_lateral_pixel_size



def _supercell_lateral_pixel_size(sample_specification):
    Delta_X, Delta_Y, _ = _supercell_dims(sample_specification)
    N_x, N_y = _supercell_xy_dims_in_pixels(sample_specification)

    sample_supercell_lateral_pixel_size = (Delta_X/N_x, Delta_Y/N_y)

    return sample_supercell_lateral_pixel_size



def supercell_slice_thickness(sample_specification,
                              skip_validation_and_conversion=\
                              _default_skip_validation_and_conversion):
    r"""Calculate the sample supercell's slice thickness.

    See the documentation for class :class:`prismatique.discretization.Params`
    for a discussion on sample supercells and how they are discretized in
    real-space.

    Parameters
    ----------
    sample_specification : :class:`prismatique.sample.ModelParams` | :class:`prismatique.sample.PotentialSliceSubsetIDs` | :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`
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

          * discretization_params

            + num_slices

        Otherwise, if ``sample_specification`` is an instance of the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs`, or the class
        :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`, then
        ``sample_specification`` specifies a set of files, where each file
        stores either the pre-calculated potential slices or :math:`S`-matrices
        for a frozen phonon configuration subset. See the documentation for the
        aforementioned classes for further discussions on specifying
        pre-calculated objects. See the documentation for the subpackage
        :mod:`prismatique.stem` for a discussion on :math:`S`-matrices.
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
    sample_supercell_slice_thickness : `float`
        The sample supercell's slice thickness in units of angstroms, denoted as
        :math:`\delta z` in the documentation for the class
        :class:`prismatique.discretization.Params`.

    """
    params = locals()

    func_alias = _check_and_convert_skip_validation_and_conversion
    skip_validation_and_conversion = func_alias(params)

    if (skip_validation_and_conversion == False):
        params["accepted_types"] = (ModelParams,
                                    PotentialSliceSubsetIDs,
                                    PotentialSliceAndSMatrixSubsetIDs)
        global_symbol_table = globals()
        for param_name in params:
            if ((param_name in "skip_validation_and_conversion")
                or (param_name in "accepted_types")):
                continue
            func_name = "_check_and_convert_" + param_name
            func_alias = global_symbol_table[func_name]
            params[param_name] = func_alias(params)
        del params["accepted_types"]

    kwargs = params
    del kwargs["skip_validation_and_conversion"]
    sample_supercell_slice_thickness = _supercell_slice_thickness(**kwargs)

    return sample_supercell_slice_thickness



def _supercell_slice_thickness(sample_specification):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    if "atomic_coords_filename" in sample_specification_core_attrs:
        _, _, Delta_Z = \
            _supercell_dims(sample_specification)
        discretization_params = \
            sample_specification_core_attrs["discretization_params"]

        discretization_params_core_attrs = \
            discretization_params.get_core_attrs(deep_copy=False)
        num_slices = \
            discretization_params_core_attrs["num_slices"]

        eps = 1e-6
        sample_supercell_slice_thickness = Delta_Z/num_slices + eps
    else:
        core_attr_name = \
            ("filenames"
             if ("filenames" in sample_specification_core_attrs)
             else "potential_slice_subset_filenames")
        filename = \
            sample_specification_core_attrs[core_attr_name][0]
        
        sample_supercell_slice_thickness = \
            _load_hdf5_attr("s", "", 0, filename)  # 2nd arg not used here.

    return sample_supercell_slice_thickness



def num_slices(sample_specification,
               skip_validation_and_conversion=\
               _default_skip_validation_and_conversion):
    r"""Calculate the number of slices :math:`N_{\text{slices}}` used to 
    partition the sample.

    See the documentation for class :class:`prismatique.discretization.Params`
    for a discussion on sample supercells and how they are discretized/sliced in
    real-space.

    Parameters
    ----------
    sample_specification : :class:`prismatique.sample.ModelParams` | :class:`prismatique.sample.PotentialSliceSubsetIDs` | :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`
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

          * discretization_params

            + num_slices

        Otherwise, if ``sample_specification`` is an instance of the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs`, or the class
        :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`, then
        ``sample_specification`` specifies a set of files, where each file
        stores either the pre-calculated potential slices or :math:`S`-matrices
        for a frozen phonon configuration subset. See the documentation for the
        aforementioned classes for further discussions on specifying
        pre-calculated objects. See the documentation for the subpackage
        :mod:`prismatique.stem` for a discussion on :math:`S`-matrices.
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
    num_sample_supercell_slices : `int`        
        The number of slices :math:`N_{\text{slices}}` used to partition the
        sample.

    """
    params = locals()

    func_alias = _check_and_convert_skip_validation_and_conversion
    skip_validation_and_conversion = func_alias(params)

    if (skip_validation_and_conversion == False):
        params["accepted_types"] = (ModelParams,
                                    PotentialSliceSubsetIDs,
                                    PotentialSliceAndSMatrixSubsetIDs)
        global_symbol_table = globals()
        for param_name in params:
            if ((param_name in "skip_validation_and_conversion")
                or (param_name in "accepted_types")):
                continue
            func_name = "_check_and_convert_" + param_name
            func_alias = global_symbol_table[func_name]
            params[param_name] = func_alias(params)
        del params["accepted_types"]

    kwargs = params
    del kwargs["skip_validation_and_conversion"]
    num_sample_supercell_slices = _num_slices(**kwargs)

    return num_sample_supercell_slices



def _num_slices(sample_specification):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    if "atomic_coords_filename" in sample_specification_core_attrs:
        discretization_params = \
            sample_specification_core_attrs["discretization_params"]
        discretization_params_core_attrs = \
            discretization_params.get_core_attrs(deep_copy=False)
        num_sample_supercell_slices = \
            discretization_params_core_attrs["num_slices"]
    else:
        _, _, Delta_Z = _supercell_dims(sample_specification)

        sample_supercell_slice_thickness = \
            (Delta_Z
             if isinstance(sample_specification, SMatrixSubsetIDs)
             else _supercell_slice_thickness(sample_specification))

        num_sample_supercell_slices = \
            int(np.ceil(Delta_Z / sample_supercell_slice_thickness))

    return num_sample_supercell_slices



def num_frozen_phonon_config_subsets(sample_specification,
                                     skip_validation_and_conversion=\
                                     _default_skip_validation_and_conversion):
    r"""Calculate the number of frozen phonon configurations subsets in a
    given sample model.

    See the documentation for the class :class:`prismatique.thermal.Params` for
    a discussion on frozen phonon configuration subsets.

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

          * thermal_params

            + num_subsets

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
    num_subsets : `int`
        The number of frozen phonon configurations subsets in the given sample
        model.

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
    num_subsets = _num_frozen_phonon_config_subsets(**kwargs)

    return num_subsets



def _num_frozen_phonon_config_subsets(sample_specification):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    if "atomic_coords_filename" in sample_specification_core_attrs:
        thermal_params = sample_specification_core_attrs["thermal_params"]
        
        thermal_params_core_attrs = \
            thermal_params.get_core_attrs(deep_copy=False)
        
        num_subsets = thermal_params_core_attrs["num_subsets"]
    else:
        num_subsets = 0
        core_attr_names_to_skip = ("interpolation_factors",
                                   "max_num_frozen_phonon_configs_per_subset")
        for core_attr_name in sample_specification_core_attrs:
            if core_attr_name in core_attr_names_to_skip:
                continue
            num_subsets += len(sample_specification_core_attrs[core_attr_name])
        
    return num_subsets



def _check_and_convert_subset_idx(params):
    obj_name = "subset_idx"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    subset_idx = czekitout.convert.to_nonnegative_int(**kwargs)

    kwargs = {"sample_specification": params["sample_specification"],
              "skip_validation_and_conversion": False}
    num_subsets = num_frozen_phonon_config_subsets(**kwargs)

    current_func_name = "_check_and_convert_subset_idx"

    if subset_idx >= num_subsets:
        err_msg = globals()[current_func_name+"_err_msg_1"]
        raise ValueError(err_msg)
    
    return subset_idx



_default_subset_idx = 0



def num_frozen_phonon_configs_in_subset(
        sample_specification,
        subset_idx=\
        _default_subset_idx,
        skip_validation_and_conversion=\
        _default_skip_validation_and_conversion):
    r"""Calculate the number of frozen phonon configurations in a given subset 
    of a given sample model.

    See the documentation for the class :class:`prismatique.thermal.Params` for
    a discussion on frozen phonon configurations and subsets thereof.

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

          * thermal_params

            + num_frozen_phonon_configs_per_subset

            + num_subsets

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
    subset_idx : `int`, optional
        The index specifying the frozen phonon configuration subset of interest.
        Must satisfy ``0<=subset_idx<prismatique.sample.num_frozen_phonon_config_subsets(sample_specification)``.
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
    num_configs_in_subset : `int`
        The number of frozen phonon configurations in the given subset of the
        given sample model.

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
    num_configs_in_subset = _num_frozen_phonon_configs_in_subset(**kwargs)

    return num_configs_in_subset



def _num_frozen_phonon_configs_in_subset(sample_specification, subset_idx):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)
    
    if "atomic_coords_filename" in sample_specification_core_attrs:
        thermal_params = \
            sample_specification_core_attrs["thermal_params"]
        thermal_params_core_attrs = \
            thermal_params.get_core_attrs(deep_copy=False)

        core_attr_name = "num_frozen_phonon_configs_per_subset"
        num_configs_in_subset = thermal_params_core_attrs[core_attr_name]
    else:
        core_attr_name = "potential_slice_subset_filenames"
        if core_attr_name in sample_specification_core_attrs:
            filenames = sample_specification_core_attrs[core_attr_name]
            if subset_idx < len(filenames):
                idx = subset_idx
                func_to_eval = _num_valid_potential_slice_datasets_in_file
            else:
                core_attr_name = "S_matrix_subset_filenames"
                idx = subset_idx - len(filenames)
                func_to_eval = _num_valid_S_matrix_datasets_in_file
        else:
            core_attr_name = "filenames"
            idx = subset_idx
            if "interpolation_factors" in sample_specification_core_attrs:
                func_to_eval = _num_valid_potential_slice_datasets_in_file
            else:
                func_to_eval = _num_valid_S_matrix_datasets_in_file

        kwargs = {"idx": \
                  idx,
                  "filename": \
                  sample_specification_core_attrs[core_attr_name][idx],
                  "core_attr_name": \
                  core_attr_name}
        num_configs_in_subset = func_to_eval(**kwargs)

        core_attr_name = \
            "max_num_frozen_phonon_configs_per_subset"
        max_num_frozen_phonon_configs_per_subset = \
            sample_specification_core_attrs[core_attr_name]

        num_configs_in_subset = min(num_configs_in_subset,
                                    max_num_frozen_phonon_configs_per_subset)

    return num_configs_in_subset



def total_num_frozen_phonon_configs(sample_specification,
                                    skip_validation_and_conversion=\
                                    _default_skip_validation_and_conversion):
    r"""Calculate the total number of frozen phonon configurations a given 
    sample model.

    See the documentation for the class :class:`prismatique.thermal.Params` for
    a discussion on frozen phonon configurations.

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

          * thermal_params

            + num_frozen_phonon_configs_per_subset

            + num_subsets

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
    total_num_configs : `int`
        The total number of frozen phonon configurations a given sample model.

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
    total_num_configs = _total_num_frozen_phonon_configs(**kwargs)

    return total_num_configs



def _total_num_frozen_phonon_configs(sample_specification):
    total_num_configs = 0
    num_subsets = _num_frozen_phonon_config_subsets(sample_specification)
    for subset_idx in range(num_subsets):
        kwargs = {"sample_specification": sample_specification,
                  "subset_idx": subset_idx}
        total_num_configs += _num_frozen_phonon_configs_in_subset(**kwargs)

    return total_num_configs



_default_probe_model_params = None



def S_matrix_k_xy_vectors(sample_specification,
                          probe_model_params=\
                          _default_probe_model_params,
                          skip_validation_and_conversion=\
                          _default_skip_validation_and_conversion):
    r"""Determine the :math:`\mathbf{k}_{xy}`-momentum vectors of the plane
    waves used to calculate the :math:`S`-matrix 
    :math:`S_{m_{x},m_{y}}\left(x,y\right)`.

    See the documentation for the subpackage :mod:`prismatique.stem` for a
    discussion on :math:`S`-matrices. As discussed therein, each pair
    :math:`\left(m_x, m_y\right)` corresponds to a different
    :math:`\mathbf{k}_{xy}`-momentum vector. In ``prismatic``, each pair
    :math:`\left(m_x, m_y\right)` is essentially mapped to a unique integer
    :math:`i`.

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

          * discretization_params

            + sample_supercell_reduced_xy_dims_in_pixels

            + interpolation_factors

        Otherwise, if ``sample_specification`` is an instance of the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs`, the class
        :class:`prismatique.sample.SMatrixSubsetIDs`, or the class
        :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`, then
        ``sample_specification`` specifies a set of files, where each file
        stores either the pre-calculated potential slices or :math:`S`-matrices
        for a frozen phonon configuration subset. See the documentation for the
        aforementioned classes for further discussions on specifying
        pre-calculated objects. See the documentation for the subpackage
        :mod:`prismatique.stem` for a discussion on :math:`S`-matrices. Note
        that any :math:`S`-matrices specified in ``sample_specification`` must
        be pre-calculated using the same convergence semiangle and mean beam
        energy as that specified by the model parameters of the probe
        ``probe_model_params`` below.
    probe_model_params : :class:`embeam.stem.probe.ModelParams` | `None`, optional
        The model parameters of the probe. See the documentation for the class
        :class:`embeam.stem.probe.ModelParams` for a discussion on said
        parameters. If ``probe_model_params`` is set to ``None`` [i.e. the
        default value], then the parameter will be reassigned to the value
        ``embeam.stem.probe.ModelParams()``.

        Note that of parameters stored in ``probe_model_params``, only the
        following are used:

        - probe_model_params

          * convergence_semiangle

          * gun_model_params

            + mean_beam_energy

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
    k_xy_vectors : `array_like` (`float`, shape=(``num_vectors``, ``2``))
        If we let ``num_vectors`` be the number of
        :math:`\mathbf{k}_{xy}`-momentum vectors, then ``k_xy_vectors[i][0]``
        and ``k_xy_vectors[i][1]`` are the :math:`x`- and :math:`y`-components
        of the ``i`` th :math:`\mathbf{k}_{xy}`-momentum vector in units of 1/Å,
        where ``0<=i<num_vectors``.

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
    k_xy_vectors = _S_matrix_k_xy_vectors(**kwargs)
    
    return k_xy_vectors



def potential_slice_subset_data_size(sample_specification,
                                     subset_idx=\
                                     _default_subset_idx,
                                     skip_validation_and_conversion=\
                                     _default_skip_validation_and_conversion):
    r"""Calculate the data size of a specified subset of potential slices that 
    one could generate according to a given sample model.

    See the documentation for the classes
    :class:`prismatique.discretization.Params` and
    :class:`prismatique.thermal.Params` for discussions on potential slices and
    frozen phonon configuration subsets respectively.

    Note that data size due to HDF5 file overhead and metadata are not taken
    into account.

    Parameters
    ----------
    sample_specification : :class:`prismatique.sample.ModelParams` | :class:`prismatique.sample.PotentialSliceSubsetIDs`
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

          * discretization_params

            + sample_supercell_reduced_xy_dims_in_pixels

            + interpolation_factors

            + num_slices

          * thermal_params

            + num_frozen_phonon_configs_per_subset

            + num_subsets

        Otherwise, if ``sample_specification`` is an instance of the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs` then
        ``sample_specification`` specifies a set of files, where each file
        stores the pre-calculated potential slices for a frozen phonon
        configuration subset. See the documentation for the aforementioned
        class for a further discussion on specifying pre-calculated
        potential slices. 
    subset_idx : `int`, optional
        The index specifying the frozen phonon configuration subset of which we
        want to calculate the size of the potential slice file.  Must satisfy
        ``0<=subset_idx<prismatique.sample.num_frozen_phonon_config_subsets(sample_specification)``.
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
    data_size : `int`
        The data size in units of bytes.

    """
    params = locals()

    func_alias = _check_and_convert_skip_validation_and_conversion
    skip_validation_and_conversion = func_alias(params)

    if (skip_validation_and_conversion == False):
        params["accepted_types"] = (ModelParams, PotentialSliceSubsetIDs)
        global_symbol_table = globals()
        for param_name in params:
            if ((param_name in "skip_validation_and_conversion")
                or (param_name in "accepted_types")):
                continue
            func_name = "_check_and_convert_" + param_name
            func_alias = global_symbol_table[func_name]
            params[param_name] = func_alias(params)
        del params["accepted_types"]

    kwargs = params
    del kwargs["skip_validation_and_conversion"]
    data_size = _potential_slice_subset_data_size(**kwargs)
        
    return data_size



def _potential_slice_subset_data_size(sample_specification, subset_idx):
    sample_supercell_xy_dims_in_pixels = \
        _supercell_xy_dims_in_pixels(sample_specification)
    num_sample_supercell_slices = \
        _num_slices(sample_specification)
    num_configs_in_subset = \
        _num_frozen_phonon_configs_in_subset(sample_specification, subset_idx)

    num_elems = (sample_supercell_xy_dims_in_pixels[0]
                 * sample_supercell_xy_dims_in_pixels[1]
                 * num_sample_supercell_slices
                 * num_configs_in_subset)

    size_of_single = 4  # In bytes.
    
    data_size = num_elems * size_of_single
        
    return data_size



def potential_slice_set_data_size(sample_specification,
                                  skip_validation_and_conversion=\
                                  _default_skip_validation_and_conversion):
    r"""Calculate the data size of a set of potential slices that one could 
    generate according to a given sample model.

    See the documentation for the class
    :class:`prismatique.discretization.Params` for a discussion on potential
    slices.

    Note that data size due to HDF5 file overhead and metadata are not taken
    into account.

    Parameters
    ----------
    sample_specification : :class:`prismatique.sample.ModelParams` | :class:`prismatique.sample.PotentialSliceSubsetIDs`
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

          * discretization_params

            + sample_supercell_reduced_xy_dims_in_pixels

            + interpolation_factors

            + num_slices

          * thermal_params

            + num_frozen_phonon_configs_per_subset

            + num_subsets

        Otherwise, if ``sample_specification`` is an instance of the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs` then
        ``sample_specification`` specifies a set of files, where each file
        stores the pre-calculated potential slices for a frozen phonon
        configuration subset. See the documentation for the aforementioned
        class for a further discussion on specifying pre-calculated
        potential slices. 
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
    data_size : `int`
        The data size in units of bytes.

    """
    params = locals()

    func_alias = _check_and_convert_skip_validation_and_conversion
    skip_validation_and_conversion = func_alias(params)

    if (skip_validation_and_conversion == False):
        params["accepted_types"] = (ModelParams,
                                    PotentialSliceSubsetIDs)
        global_symbol_table = globals()
        for param_name in params:
            if ((param_name in "skip_validation_and_conversion")
                or (param_name in "accepted_types")):
                continue
            func_name = "_check_and_convert_" + param_name
            func_alias = global_symbol_table[func_name]
            params[param_name] = func_alias(params)
        del params["accepted_types"]

    kwargs = params
    del kwargs["skip_validation_and_conversion"]
    data_size = _potential_slice_set_data_size(**kwargs)
        
    return data_size



def _potential_slice_set_data_size(sample_specification):
    data_size = 0
    num_subsets = _num_frozen_phonon_config_subsets(sample_specification)
    for subset_idx in range(num_subsets):
        kwargs = {"sample_specification": sample_specification,
                  "subset_idx": subset_idx}
        data_size += _potential_slice_subset_data_size(**kwargs)

    data_size = int(data_size)
        
    return data_size



def S_matrix_subset_data_size(sample_specification,
                              probe_model_params=\
                              _default_probe_model_params,
                              subset_idx=\
                              _default_subset_idx,
                              skip_validation_and_conversion=\
                              _default_skip_validation_and_conversion):
    r"""Calculate the data size of a specified subset of :math:`S`-matrices that
    one could generate according to a given sample model and probe model.

    See the documentation for the subpackage :mod:`prismatique.stem` and the
    class :class:`prismatique.thermal.Params` for discussions on
    :math:`S`-matrices and frozen phonon configurations subsets respectively.

    Note that data size due to HDF5 file overhead and metadata are not taken
    into account.

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

          * discretization_params

            + sample_supercell_reduced_xy_dims_in_pixels

            + interpolation_factors

          * thermal_params

            + num_frozen_phonon_configs_per_subset

            + num_subsets

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

    probe_model_params : :class:`embeam.stem.probe.ModelParams` | `None`, optional 
        The model parameters of the probe. See the documentation for the class
        :class:`embeam.stem.probe.ModelParams` for a discussion on said
        parameters. If ``probe_model_params`` is set to ``None`` [i.e. the
        default value], then the parameter will be reassigned to the value
        ``embeam.stem.probe.ModelParams()``.

        Note that of parameters stored in ``probe_model_params``, only the
        following are used:

        - probe_model_params

          * convergence_semiangle

          * gun_model_params

            + mean_beam_energy

    subset_idx : `int`, optional
        The index specifying the frozen phonon configuration subset of which we
        want to calculate the size of the :math:`S`-matrix file.  Must satisfy
        ``0<=subset_idx<prismatique.sample.num_frozen_phonon_config_subsets(sample_specification)``.
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
    data_size : `int`
        The data size in units of bytes.

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
    data_size = _S_matrix_subset_data_size(**kwargs)
        
    return data_size



def _S_matrix_subset_data_size(sample_specification,
                               probe_model_params,
                               subset_idx):
    sample_supercell_xy_dims_in_pixels = \
        _supercell_xy_dims_in_pixels(sample_specification)
    interpolation_factors = \
        _interpolation_factors_from_sample_specification(sample_specification)
    num_S_matrix_k_xy_vectors = \
        len(_S_matrix_k_xy_vectors(sample_specification, probe_model_params))
    num_configs_in_subset = \
        _num_frozen_phonon_configs_in_subset(sample_specification, subset_idx)

    num_elems = (sample_supercell_xy_dims_in_pixels[0]
                 * sample_supercell_xy_dims_in_pixels[1]
                 * num_S_matrix_k_xy_vectors
                 * num_configs_in_subset
                 / interpolation_factors[0] / 2
                 / interpolation_factors[1] / 2)

    size_of_single = 4  # In bytes.
    size_of_complex_single = 2 * size_of_single  # In bytes.
    
    data_size = num_elems * size_of_complex_single
        
    return data_size



def S_matrix_set_data_size(sample_specification,
                           probe_model_params=\
                           _default_probe_model_params,
                           skip_validation_and_conversion=\
                           _default_skip_validation_and_conversion):
    r"""Calculate the data size of a set of :math:`S`-matrices that one could 
    generate according to a given sample model and probe model.

    See the documentation for the subpackage :mod:`prismatique.stem` for a
    discussion on :math:`S`-matrices.

    Note that data size due to HDF5 file overhead and metadata are not taken
    into account.

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

          * discretization_params

            + sample_supercell_reduced_xy_dims_in_pixels

            + interpolation_factors

          * thermal_params

            + num_frozen_phonon_configs_per_subset

            + num_subsets

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

    probe_model_params : :class:`embeam.stem.probe.ModelParams` | `None`, optional 
        The model parameters of the probe. See the documentation for the class
        :class:`embeam.stem.probe.ModelParams` for a discussion on said
        parameters. If ``probe_model_params`` is set to ``None`` [i.e. the
        default value], then the parameter will be reassigned to the value
        ``embeam.stem.probe.ModelParams()``.

        Note that of parameters stored in ``probe_model_params``, only the
        following are used:

        - probe_model_params

          * convergence_semiangle

          * gun_model_params

            + mean_beam_energy

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
    data_size : `int`
        The data size in units of bytes.

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
    data_size = _S_matrix_set_data_size(**kwargs)
        
    return data_size



def _S_matrix_set_data_size(sample_specification, probe_model_params):
    data_size = 0
    num_subsets = _num_frozen_phonon_config_subsets(sample_specification)
    for subset_idx in range(num_subsets):
        data_size += _S_matrix_subset_data_size(sample_specification,
                                                probe_model_params,
                                                subset_idx)

    data_size = int(data_size)
        
    return data_size



def _check_and_convert_output_dirname(params):
    obj_name = "output_dirname"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    output_dirname = czekitout.convert.to_str_from_str_like(**kwargs)

    if "worker_params" in params:
        _ = _check_and_convert_max_data_size(params)
        _ = _check_and_convert_worker_params(params)

        sample_specification = (params["sample_model_params"]
                                if ("sample_model_params" in params)
                                else params["sample_specification"])

        unformatted_basename = ("potential_slices_of_subset_{}.h5"
                                if ("sample_model_params" in params)
                                else "S_matrices_of_subset_{}.h5")
        
        kwargs = {"sample_specification": sample_specification,
                  "output_dirname": output_dirname,
                  "unformatted_basename": unformatted_basename}
        _pre_save(**kwargs)
    
    return output_dirname



def _pre_save(sample_specification, output_dirname, unformatted_basename):
    num_frozen_phonon_config_subsets = \
        _num_frozen_phonon_config_subsets(sample_specification)
    
    filenames = tuple()
    for subset_idx in range(num_frozen_phonon_config_subsets):
        basename = unformatted_basename.format(subset_idx)
        filenames += (output_dirname + "/" + basename,)
    _check_hdf5_filenames(filenames)

    pathlib.Path(output_dirname).mkdir(parents=True, exist_ok=True)

    return None



def _check_hdf5_filenames(filenames):
    current_func_name = "_check_hdf5_filenames"
    
    for filename in filenames:
        kwargs = {"group": None,
                  "group_id": h5pywrappers.obj.ID(filename, "/"),
                  "write_mode": "w"}
        h5pywrappers.group.save(**kwargs)

    return None



def _check_and_convert_max_data_size(params):
    obj_name = "max_data_size"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    max_data_size = czekitout.convert.to_positive_int(**kwargs)

    if "worker_params" in params:
        if "probe_model_params" in params:
            probe_model_params = _check_and_convert_probe_model_params(params)
            
            sample_specification = (params["sample_model_params"]
                                    if ("sample_model_params" in params)
                                    else params["sample_specification"])
        else:
            probe_model_params = None

            unformatted_func_name = "_check_and_convert_sample_{}"
            func_name = (unformatted_func_name.format("model_params")
                         if ("sample_model_params" in params)
                         else unformatted_func_name.format("specification"))
            func_alias = globals()[func_name]
            sample_specification = func_alias(params)
        
        kwargs = {"sample_specification": sample_specification,
                  "probe_model_params": probe_model_params,
                  "max_data_size": max_data_size}
        _check_data_size(**kwargs)
    
    return max_data_size



def _check_data_size(sample_specification, probe_model_params, max_data_size):
    kwargs = {"obj": max_data_size, "obj_name": "max_data_size"}
    max_data_size = czekitout.convert.to_positive_int(**kwargs)
    
    if probe_model_params is None:
        objs_to_generate = "potential slices"
        data_size = _potential_slice_set_data_size(sample_specification)
    else:
        objs_to_generate = "S-matrices"
        data_size = _S_matrix_set_data_size(sample_specification,
                                            probe_model_params)

    current_func_name = "_check_data_size"
        
    if max_data_size < data_size:
        unformatted_err_msg = globals()[current_func_name+"_err_msg_1"]
        err_msg = unformatted_err_msg.format(objs_to_generate,
                                             data_size,
                                             max_data_size)
        raise MemoryError(err_msg)

    return None



def _check_and_convert_worker_params(params):
    module_alias = prismatique.worker
    func_alias = module_alias._check_and_convert_worker_params
    worker_params = func_alias(params)

    return worker_params



_module_alias = prismatique.worker
_default_output_dirname = "potential_slice_generator_output"
_default_max_data_size = 2*10**9
_default_worker_params = _module_alias._default_worker_params



def generate_potential_slices(sample_model_params,
                              output_dirname=\
                              _default_output_dirname,
                              max_data_size=\
                              _default_max_data_size,
                              worker_params=\
                              _default_worker_params,
                              skip_validation_and_conversion=\
                              _default_skip_validation_and_conversion):
    r"""Generate the potential slices for a given sample model.

    For each frozen phonon configuration subset, this Python function generates
    the corresponding potential slices and saves them to a set of files. See the
    documentation for the classes :class:`prismatique.discretization.Params` and
    :class:`prismatique.thermal.Params` for discussions on potential slices and
    frozen phonon configuration subsets respectively.

    For each frozen phonon configuration subset, the corresponding potential
    slice data is written to a file with the basename
    ``"potential_slices_of_subset_"+str(i)+".h5"``, with ``i`` being the subset
    index.

    In addition to the potential slice files, two other files are generated as
    well: The first is a JSON file with the basename
    ``"sample_model_params.json"`` which contains, in a serialized format, the
    simulation parameters related to the modelling of the sample; the second
    file is a JSON file with the basename ``"worker_params.json"`` which
    contains, in a serialized format, the simulation parameters related to GPU
    and CPU workers.

    Parameters
    ----------
    sample_model_params : :class:`prismatique.sample.ModelParams`
        The simulation parameters related to the modelling of the sample.
    output_dirname : `str`, optional
        The relative or absolute path to the directory in which all output files
        are to be saved. If the directory doesn't exist upon saving the output 
        files, it will be created if possible.
    max_data_size : `int`, optional
        The data size limit, in bytes, of the potential slices to be
        generated. If the potential slices to be generated would require a data
        size larger than the aforementioned limit, then an exception is raised
        and the potential slices are not generated. Note that data size due to
        HDF5 file overhead and metadata are not taken into account.
    worker_params : :class:`prismatique.worker.Params` | `None`, optional
        The simulation parameters related to GPU and CPU workers. See the
        documentation for the class :class:`prismatique.worker.Params` for a
        discussion on said parameters. If ``worker_params`` is set to `None`
        [i.e. the default value], then the aforementioned simulation parameters 
        are set to default values.
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
    _generate_potential_slices(**kwargs)
    
    return None



def _generate_potential_slices(sample_model_params,
                               output_dirname,
                               worker_params,
                               max_data_size):
    kwargs = {"sample_specification": sample_model_params}
    rng_seeds = _generate_rng_seeds_from_sample_specification(**kwargs)

    output_param_subset = {"output_dirname": output_dirname,
                           "save_potential_slices": True,
                           "save_S_matrices": False}

    kwargs = {"output_param_subset": output_param_subset,
              "subset_idx": 0,
              "first_or_last_call": True}
    _remove_temp_files(**kwargs)
    
    func_alias = \
        _initialize_pyprismatic_sim_obj_for_potential_slice_generation
    kwargs = \
        {"sample_model_params": sample_model_params,
         "output_dirname": output_dirname,
         "worker_params": worker_params}
    pyprismatic_sim_obj = \
        func_alias(**kwargs)

    for subset_idx, rng_seed in enumerate(rng_seeds):
        kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
                  "sample_model_params": sample_model_params,
                  "output_param_subset": output_param_subset,
                  "subset_idx": subset_idx,
                  "rng_seeds": rng_seeds}
        _generate_potential_slice_subset(**kwargs)

    kwargs = {"filename": output_dirname + "/sample_model_params.json",
              "overwrite": True}
    sample_model_params.dump(**kwargs)
    
    kwargs["filename"] = output_dirname + "/worker_params.json"
    worker_params.dump(**kwargs)
    
    return None



def _remove_temp_files(output_param_subset, subset_idx, first_or_last_call):
    output_dirname = output_param_subset["output_dirname"]

    filenames = ("scratch_param.txt",
                 "prismatic_scratch.h5",
                 output_dirname + "/prismatic_output.h5")

    temp_scan_config_filename = \
        _generate_temp_scan_config_filename(output_dirname)
    temp_aberration_filename = \
        _generate_temp_aberration_filename(output_dirname)
    temp_atomic_coords_filename = \
        _generate_temp_atomic_coords_filename(output_dirname)
    
    if first_or_last_call:
        filenames += (temp_scan_config_filename,)
        filenames += (temp_aberration_filename,)
        filenames += (temp_atomic_coords_filename,)
        save_S_matrices = False
        save_potential_slices = False
    else:
        save_S_matrices = output_param_subset["save_S_matrices"]
        save_potential_slices = output_param_subset["save_potential_slices"]

    for filename in filenames:
        if pathlib.Path(filename).is_file():
            pathlib.Path(filename).unlink(missing_ok=True)
    
    precalculated_types = ("S_matrices", "potential_slices")
    save = (save_S_matrices, save_potential_slices)
    
    for idx in range(2):
        temp_filename = (output_dirname + "/"
                         + precalculated_types[idx]  + "_temp.h5")
        if pathlib.Path(temp_filename).is_file():
            if save[idx]:
                new_filename = (output_dirname + "/"
                                + precalculated_types[idx]
                                + "_of_subset_{}.h5".format(subset_idx))
                pathlib.Path(new_filename).unlink(missing_ok=True)
                pathlib.Path(temp_filename).rename(new_filename)
            else:
                pathlib.Path(temp_filename).unlink(missing_ok=True)
        
    return None



def _generate_temp_scan_config_filename(output_dirname):
    temp_scan_config_filename = (output_dirname
                                 + "/_temp_scan_config_filename.txt")

    return temp_scan_config_filename



def _generate_temp_aberration_filename(output_dirname):
    temp_aberration_filename = (output_dirname
                                + "/_temp_aberration_filename.txt")

    return temp_aberration_filename



def _generate_temp_atomic_coords_filename(output_dirname):
    temp_atomic_coords_filename = (output_dirname
                                   + "/_temp_atomic_coords_filename.txt")

    return temp_atomic_coords_filename



def _initialize_pyprismatic_sim_obj_for_potential_slice_generation(
        sample_model_params, output_dirname, worker_params):
    pyprismatic_sim_obj = pyprismatic.Metadata()

    _set_pyprismatic_sim_obj_attrs_to_default_values(pyprismatic_sim_obj,
                                                     output_dirname)

    sample_specification = sample_model_params
    _unpack_sample_specification_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                          output_dirname,
                                                          sample_specification)

    pyprismatic_sim_obj.filenameOutput = output_dirname + "/prismatic_output.h5"
    
    _unpack_worker_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                   worker_params)

    return pyprismatic_sim_obj



def _set_pyprismatic_sim_obj_attrs_to_default_values(pyprismatic_sim_obj,
                                                     output_dirname):
    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "output_dirname": output_dirname,
              "sample_specification": None}
    _unpack_sample_specification_into_pyprismatic_sim_obj(**kwargs)

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "probe_model_params": embeam.stem.probe.ModelParams(),
              "output_dirname": output_dirname}
    _unpack_probe_model_params_into_pyprismatic_sim_obj(**kwargs)

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "tilt_params": prismatique.tilt.Params()}
    _unpack_tilt_params_into_pyprismatic_sim_obj(**kwargs)

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "scan_config": "",
              "output_dirname": output_dirname}
    _unpack_scan_config_into_pyprismatic_sim_obj(**kwargs)

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "output_params": None,
              "sample_specification": None}
    _unpack_output_params_into_pyprismatic_sim_obj(**kwargs)

    pyprismatic_sim_obj.xTiltStep = 0.0
    pyprismatic_sim_obj.yTiltStep = 0.0

    pyprismatic_sim_obj.probeXtilt = 0.0
    pyprismatic_sim_obj.probeYtilt = 0.0

    pyprismatic_sim_obj.numFP = 1  # Updated elsewhere.
    pyprismatic_sim_obj.randomSeed = 1  # Updated elsewhere.

    return None



def _unpack_sample_specification_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                          output_dirname,
                                                          sample_specification):
    temp_atomic_coords_filename = \
        _generate_temp_atomic_coords_filename(output_dirname)

    if sample_specification is not None:
        sample_specification_core_attrs = \
            sample_specification.get_core_attrs(deep_copy=False)
        
        if "atomic_coords_filename" in sample_specification_core_attrs:
            sample_model_params = sample_specification
        else:
            func_alias = _interpolation_factors_from_sample_specification
            f_x, f_y = func_alias(sample_specification)
            
            N_x, N_y = _supercell_xy_dims_in_pixels(sample_specification)
            
            kwargs = {"interpolation_factors": \
                      (f_x, f_y),
                      "num_slices": \
                      _num_slices(sample_specification),
                      "sample_supercell_reduced_xy_dims_in_pixels": \
                      (N_x / 4 / f_x, N_y / 4 / f_y)}
            discretization_params = prismatique.discretization.Params(**kwargs)
            
            _generate_trivial_atomic_coords_file(sample_specification,
                                                 output_dirname)
            
            kwargs = {"atomic_coords_filename": temp_atomic_coords_filename,
                      "discretization_params": discretization_params}
            sample_model_params = ModelParams(**kwargs)
    else:
        _generate_trivial_atomic_coords_file(sample_specification,
                                             output_dirname)

        kwargs = {"atomic_coords_filename": temp_atomic_coords_filename}
        sample_model_params = ModelParams(**kwargs)

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "sample_model_params": sample_model_params}
    _unpack_sample_model_params_into_pyprismatic_sim_obj(**kwargs)

    pyprismatic_sim_obj.includeOccupancy = False
    pyprismatic_sim_obj.importPath = ""  # Ignored.    
    pyprismatic_sim_obj.importPotential = False  # Updated elsewhere, if at all.
    pyprismatic_sim_obj.importSMatrix = False  # Updated elsewhere, if at all.
    pyprismatic_sim_obj.importFile = ""  # Updated elsewhere, if at all.
    
    return None



def _generate_trivial_atomic_coords_file(sample_specification, output_dirname):
    try:
        sample_unit_cell_dims = _unit_cell_dims(sample_specification)
    except:
        sample_unit_cell_dims = (5, 5, 5)  # Default.

    temp_atomic_coords_filename = \
        _generate_temp_atomic_coords_filename(output_dirname)

    with open(temp_atomic_coords_filename, "w") as file_obj:
        line = "Trivial atomic coordinates file\n"
        file_obj.write(line)

        unformatted_line = "\t{:18.14f}\t{:18.14f}\t{:18.14f}\n"
        formatted_line = unformatted_line.format(*sample_unit_cell_dims)
        file_obj.write(formatted_line)

        occ = 1
        Z = 8
        x = 0
        y = 0
        z = 0
        u_x_rms = 0.1

        unformatted_line = ("{}\t{:18.14f}\t{:18.14f}"
                            "\t{:18.14f}\t{:18.14f}\t{:18.14f}\n")
        formatted_line = unformatted_line.format(Z, x, y, z, occ, u_x_rms)
        file_obj.write(formatted_line)

        file_obj.write("-1")

    return None



def _unpack_sample_model_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                         sample_model_params):
    sample_model_params_core_attrs = \
        sample_model_params.get_core_attrs(deep_copy=False)
    
    atomic_coords_filename = \
        sample_model_params_core_attrs["atomic_coords_filename"]

    check_atomic_coords_file_format(atomic_coords_filename)
    
    pyprismatic_sim_obj.realspacePixelSizeX = \
        _supercell_lateral_pixel_size(sample_model_params)[0]
    pyprismatic_sim_obj.realspacePixelSizeY = \
        _supercell_lateral_pixel_size(sample_model_params)[1]
    pyprismatic_sim_obj.sliceThickness = \
        _supercell_slice_thickness(sample_model_params)
        
    pyprismatic_sim_obj.filenameAtoms = atomic_coords_filename

    unit_cell_tiling = sample_model_params_core_attrs["unit_cell_tiling"]
    pyprismatic_sim_obj.tileX = unit_cell_tiling[0]
    pyprismatic_sim_obj.tileY = unit_cell_tiling[1]
    pyprismatic_sim_obj.tileZ = unit_cell_tiling[2]

    discretization_params = \
        sample_model_params_core_attrs["discretization_params"]
    discretization_params_core_attrs = \
        discretization_params.get_core_attrs(deep_copy=False)
    
    z_supersampling = discretization_params_core_attrs["z_supersampling"]
    pyprismatic_sim_obj.potential3D = (z_supersampling > 0)
    pyprismatic_sim_obj.zSampling = z_supersampling + 16*(z_supersampling == 0)

    f_x, f_y = discretization_params_core_attrs["interpolation_factors"]
    pyprismatic_sim_obj.interpolationFactorX = f_x
    pyprismatic_sim_obj.interpolationFactorY = f_y

    pyprismatic_sim_obj.potBound = \
        sample_model_params_core_attrs["atomic_potential_extent"]

    thermal_params = \
        sample_model_params_core_attrs["thermal_params"]
    thermal_params_core_attrs = \
        thermal_params.get_core_attrs(deep_copy=False)
    
    pyprismatic_sim_obj.includeThermalEffects = \
        thermal_params_core_attrs["enable_thermal_effects"]

    return None



def _unpack_probe_model_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                        probe_model_params,
                                                        output_dirname):
    probe_model_params_core_attrs = \
        probe_model_params.get_core_attrs(deep_copy=False)

    gun_model_params = probe_model_params_core_attrs["gun_model_params"]
    _unpack_gun_model_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                      gun_model_params)

    lens_model_params = probe_model_params_core_attrs["lens_model_params"]
    _unpack_lens_model_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                       lens_model_params,
                                                       output_dirname)

    pyprismatic_sim_obj.probeSemiangle = \
        probe_model_params_core_attrs["convergence_semiangle"]

    pyprismatic_sim_obj.alphaBeamMax = 24  # Ignored.

    return None



def _unpack_gun_model_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                      gun_model_params):
    gun_model_params_core_attrs = \
        gun_model_params.get_core_attrs(deep_copy=False)

    pyprismatic_sim_obj.E0 = gun_model_params.core_attrs["mean_beam_energy"]

    return None



def _unpack_lens_model_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                       lens_model_params,
                                                       output_dirname):
    pyprismatic_sim_obj.probeDefocus = float("NaN")
    pyprismatic_sim_obj.aberrations_file = ""
    
    pyprismatic_sim_obj.probeDefocus_min = 0.0  # Ignored.
    pyprismatic_sim_obj.probeDefocus_max = 0.0  # Ignored.
    pyprismatic_sim_obj.probeDefocus_step = 0.0  # Ignored.
    pyprismatic_sim_obj.probeDefocus_sigma = 0.0  # Ignored.

    pyprismatic_sim_obj.C3 = float("NaN")  # Ignored.
    pyprismatic_sim_obj.C5 = float("NaN")  # Ignored.

    lens_model_params_core_attrs = \
        lens_model_params.get_core_attrs(deep_copy=False)
    
    coherent_aberrations = lens_model_params_core_attrs["coherent_aberrations"]
    
    lines = ("m n C_mag C_ang\n",)
    
    for coherent_aberration in coherent_aberrations:
        coherent_aberration_core_attrs = \
            coherent_aberration.get_core_attrs(deep_copy=False)
        
        m = coherent_aberration_core_attrs["m"]
        n = coherent_aberration_core_attrs["n"]
        C_mag = coherent_aberration_core_attrs["C_mag"]
        C_ang = coherent_aberration_core_attrs["C_ang"]
        
        if (m == 2) and (n == 0):
            wavelength = embeam.wavelength(pyprismatic_sim_obj.E0)
            pyprismatic_sim_obj.probeDefocus = wavelength*C_mag/np.pi
        else:
            line = "{} {} {} {}\n".format(m, n, C_mag, C_ang*(180/np.pi))
            lines += (line,)
            
    lines += ("-1",)

    if len(lines) > 2:
        temp_aberration_filename = \
            _generate_temp_aberration_filename(output_dirname)
        pyprismatic_sim_obj.aberrations_file = \
            temp_aberration_filename
        
        with open(temp_aberration_filename, 'w') as file_obj:
            for line in lines:
                file_obj.write(line)

    return None



def _unpack_tilt_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                 tilt_params):
    tilt_params_core_attrs = \
        tilt_params.get_core_attrs(deep_copy=False)

    pyprismatic_sim_obj.xTiltOffset = tilt_params_core_attrs["offset"][0]
    pyprismatic_sim_obj.yTiltOffset = tilt_params_core_attrs["offset"][1]

    if len(tilt_params_core_attrs["window"]) == 4:
        pyprismatic_sim_obj.minXtilt = tilt_params_core_attrs["window"][0]
        pyprismatic_sim_obj.maxXtilt = tilt_params_core_attrs["window"][1]
        pyprismatic_sim_obj.minYtilt = tilt_params_core_attrs["window"][2]
        pyprismatic_sim_obj.maxYtilt = tilt_params_core_attrs["window"][3]
    else:
        pyprismatic_sim_obj.minRtilt = tilt_params_core_attrs["window"][0]
        pyprismatic_sim_obj.maxRtilt = tilt_params_core_attrs["window"][1]
    
    return None



def _unpack_scan_config_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                 scan_config,
                                                 output_dirname):
    pyprismatic_sim_obj.probes_file = \
        (scan_config
         if isinstance(scan_config, str)
         else _generate_temp_scan_config_filename(output_dirname))

    pyprismatic_sim_obj.scanWindowXMin = 0.5
    pyprismatic_sim_obj.scanWindowXMax = 0.5+1.0e-14
    pyprismatic_sim_obj.scanWindowYMin = 0.5
    pyprismatic_sim_obj.scanWindowYMax = 0.5+1.0e-14
    pyprismatic_sim_obj.scanWindowXMin_r = 0.0  # Ignored.
    pyprismatic_sim_obj.scanWindowXMax_r = 0.0  # Ignored.
    pyprismatic_sim_obj.scanWindowYMin_r = 0.0  # Ignored.
    pyprismatic_sim_obj.scanWindowYMax_r = 0.0  # Ignored.

    pyprismatic_sim_obj.probeStepX = 0.25
    pyprismatic_sim_obj.probeStepY = 0.25
    
    pyprismatic_sim_obj.nyquistSampling = False  # Ignored.

    return None



def _unpack_output_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                   output_params,
                                                   sample_specification):
    if output_params is None:
        pyprismatic_sim_obj.numSlices = 0
        pyprismatic_sim_obj.zStart = 0.0

        pyprismatic_sim_obj.algorithm = "multislice"
        pyprismatic_sim_obj.matrixRefocus = False

        pyprismatic_sim_obj.save2DOutput = False  # Ignored.
        pyprismatic_sim_obj.save3DOutput = False  # Ignored.
        pyprismatic_sim_obj.save4DOutput = False
        pyprismatic_sim_obj.saveDPC_CoM = False  # Ignored.
        pyprismatic_sim_obj.savePotentialSlices = False
        pyprismatic_sim_obj.saveSMatrix = False
        pyprismatic_sim_obj.saveComplexOutputWave = True
        pyprismatic_sim_obj.saveProbe = 0  # Ignored.

        pyprismatic_sim_obj.integrationAngleMin = 0  # Ignored.
        pyprismatic_sim_obj.integrationAngleMax = 1.0  # Ignored.
        pyprismatic_sim_obj.detectorAngleStep = 1.0  # Ignored.
    
        pyprismatic_sim_obj.crop4DOutput = False  # Ignored.
        pyprismatic_sim_obj.crop4Damax = 100.0  # Ignored.

        pyprismatic_sim_obj.maxFileSize = 300*10**9
        pyprismatic_sim_obj.filenameOutput = "prismatic_output.h5"
    else:
        output_params_core_attrs = output_params.get_core_attrs(deep_copy=False)

        if "image_params" in output_params_core_attrs:
            kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
                      "hrtem_output_params": output_params}
            _unpack_hrtem_output_params_into_pyprismatic_sim_obj(**kwargs)
        else:
            kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
                      "stem_output_params": output_params,
                      "sample_specification": sample_specification}
            _unpack_stem_output_params_into_pyprismatic_sim_obj(**kwargs)
    
    return None



def _unpack_hrtem_output_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                         hrtem_output_params):
    pyprismatic_sim_obj.algorithm = "hrtem"

    hrtem_output_params_core_attrs = \
        hrtem_output_params.get_core_attrs(deep_copy=False)

    output_dirname = hrtem_output_params_core_attrs["output_dirname"]
    filename = output_dirname + "/prismatic_output.h5"

    current_func_name = "_unpack_hrtem_output_params_into_pyprismatic_sim_obj"
    
    kwargs = {"group": None,
              "group_id": h5pywrappers.obj.ID(filename, "/"),
              "write_mode": "w"}
    h5pywrappers.group.save(**kwargs)
    
    pyprismatic_sim_obj.filenameOutput = filename

    pyprismatic_sim_obj.savePotentialSlices = \
        hrtem_output_params_core_attrs["save_potential_slices"]

    return None



def _unpack_stem_output_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                        stem_output_params,
                                                        sample_specification):
    stem_output_params_core_attrs = \
        stem_output_params.get_core_attrs(deep_copy=False)

    base_stem_output_params = stem_output_params_core_attrs["base_params"]
    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "base_stem_output_params": base_stem_output_params}
    _unpack_base_stem_output_params_into_pyprismatic_sim_obj(**kwargs)

    alg_specific_params = stem_output_params_core_attrs["alg_specific_params"]
    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "alg_specific_params": alg_specific_params,
              "sample_specification": sample_specification}
    _unpack_alg_specific_params_into_pyprismatic_sim_obj(**kwargs)

    return None



def _unpack_base_stem_output_params_into_pyprismatic_sim_obj(
        pyprismatic_sim_obj, base_stem_output_params):
    base_stem_output_params_core_attrs = \
        base_stem_output_params.get_core_attrs(deep_copy=False)

    output_dirname = base_stem_output_params_core_attrs["output_dirname"]
    filename = output_dirname + "/prismatic_output.h5"
    
    kwargs = {"group": None,
              "group_id": h5pywrappers.obj.ID(filename, "/"),
              "write_mode": "w"}
    h5pywrappers.group.save(**kwargs)
    
    pyprismatic_sim_obj.filenameOutput = filename
        
    pyprismatic_sim_obj.savePotentialSlices = \
        base_stem_output_params_core_attrs["save_potential_slices"]
    
    radial_step_size_for_3d_stem = \
        base_stem_output_params_core_attrs["radial_step_size_for_3d_stem"]
    radial_range_for_2d_stem = \
        base_stem_output_params_core_attrs["radial_range_for_2d_stem"]
    cbed_params = \
        base_stem_output_params_core_attrs["cbed_params"]

    cbed_params_core_attrs = cbed_params.get_core_attrs(deep_copy=False)

    if (cbed_params_core_attrs["save_final_intensity"]
        or cbed_params_core_attrs["save_wavefunctions"]
        or base_stem_output_params_core_attrs["save_com"]
        or (radial_step_size_for_3d_stem > 0)
        or (radial_range_for_2d_stem[0] != radial_range_for_2d_stem[1])):
        pyprismatic_sim_obj.save4DOutput = True

    return None



def _unpack_alg_specific_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                         alg_specific_params,
                                                         sample_specification):
    alg_specific_params_core_attrs = \
        alg_specific_params.get_core_attrs(deep_copy=False)

    if "num_slices_per_output" in alg_specific_params_core_attrs:
        kwargs = \
            {"sample_specification": sample_specification,
             "alg_specific_params": alg_specific_params}
        z_start_output = \
            _z_start_output_passed_to_pyprismatic(**kwargs)

        pyprismatic_sim_obj.numSlices = \
            alg_specific_params_core_attrs["num_slices_per_output"]
        pyprismatic_sim_obj.zStart = \
            z_start_output
        pyprismatic_sim_obj.algorithm = \
            "multislice"
    else:
        pyprismatic_sim_obj.matrixRefocus = \
            alg_specific_params_core_attrs["enable_S_matrix_refocus"]
        pyprismatic_sim_obj.saveSMatrix = \
            alg_specific_params_core_attrs["save_S_matrices"]
        pyprismatic_sim_obj.algorithm = \
            "prism"

    return None



def _z_start_output_passed_to_pyprismatic(sample_specification,
                                          alg_specific_params):
    alg_specific_params_core_attrs = \
        alg_specific_params.get_core_attrs(deep_copy=False)

    _, _, Delta_Z = _supercell_dims(sample_specification)
        
    sample_supercell_slice_thickness = \
        _supercell_slice_thickness(sample_specification)
        
    z_start_output = alg_specific_params_core_attrs["z_start_output"]
    z_start_output = z_start_output if (z_start_output < Delta_Z) else Delta_Z
    
    tol = 1e-10
    temp_1 = z_start_output / sample_supercell_slice_thickness
    temp_2 = temp_1 - round(temp_1)

    z_start_output = ((round(temp_1)-tol) * sample_supercell_slice_thickness
                      if ((abs(temp_2) < tol) and (temp_2 > 0))
                      else z_start_output)

    return z_start_output



def _unpack_worker_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                   worker_params):
    worker_params_core_attrs = worker_params.get_core_attrs(deep_copy=False)

    cpu_params = worker_params_core_attrs["cpu_params"]
    _unpack_cpu_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj, cpu_params)

    gpu_params = worker_params_core_attrs["gpu_params"]
    _unpack_gpu_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj, gpu_params)

    return None



def _unpack_cpu_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                cpu_params):
    cpu_params_core_attrs = cpu_params.get_core_attrs(deep_copy=False)

    pyprismatic_sim_obj.alsoDoCPUWork = \
        cpu_params_core_attrs["enable_workers"]
    pyprismatic_sim_obj.numThreads = \
        cpu_params_core_attrs["num_worker_threads"]
    pyprismatic_sim_obj.batchSizeTargetCPU = \
        cpu_params_core_attrs["batch_size"]
    pyprismatic_sim_obj.earlyCPUStopCount = \
        cpu_params_core_attrs["early_stop_count"]

    pyprismatic_sim_obj.batchSizeCPU = 1  # Ignored.
    
    return None



def _unpack_gpu_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                gpu_params):
    gpu_params_core_attrs = gpu_params.get_core_attrs(deep_copy=False)

    pyprismatic_sim_obj.numGPUs = \
        gpu_params_core_attrs["num_gpus"]
    pyprismatic_sim_obj.batchSizeTargetGPU = \
        gpu_params_core_attrs["batch_size"]
    pyprismatic_sim_obj.numStreamsPerGPU = \
        gpu_params_core_attrs["num_streams_per_gpu"]

    data_transfer_mode = \
        gpu_params_core_attrs["data_transfer_mode"]
    pyprismatic_sim_obj.transferMode = \
        (data_transfer_mode
         if (data_transfer_mode in ("streaming", "auto"))
         else "singlexfer")

    pyprismatic_sim_obj.batchSizeGPU = 1  # Ignored.
    
    return None



def _generate_rng_seeds_from_sample_specification(sample_specification):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)

    if "thermal_params" in sample_specification_core_attrs:
        thermal_params = \
            sample_specification_core_attrs["thermal_params"]
        thermal_params_core_attrs = \
            thermal_params.get_core_attrs(deep_copy=False)
        
        np.random.seed(thermal_params_core_attrs["rng_seed"])

    num_subsets = _num_frozen_phonon_config_subsets(sample_specification)
    rng_seeds = np.random.randint(low=0, high=999999, size=num_subsets)
    np.random.seed()

    return rng_seeds



def _generate_potential_slice_subset(pyprismatic_sim_obj,
                                     sample_model_params,
                                     output_param_subset,
                                     subset_idx,
                                     rng_seeds):
    output_dirname = \
        output_param_subset["output_dirname"]
    new_output_filename = \
        output_dirname + "/potential_slices_of_subset_{}.h5".format(subset_idx)
        
    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "sample_specification": sample_model_params,
              "output_param_subset": output_param_subset,
              "subset_idx": subset_idx,
              "defocus_idx": 0,
              "defocii": (0,),
              "rng_seeds": rng_seeds}
    _update_pyprismatic_sim_obj_for_next_prismatic_sim(**kwargs)

    unformatted_msg = ("Generating the potential slices for the frozen phonon "
                       "configuration subset #{}...\n")
    msg = unformatted_msg.format(subset_idx)
    print(msg)

    # Run prismatic simulation.
    _call_pyprismatic_sim_obj_go(pyprismatic_sim_obj)
    gc.collect()
    
    func_alias = _move_sample_specification_output_to_separate_file
    kwargs = {"output_dirname": output_dirname,
              "sample_specification_type": "potential_slice_subset",
              "new_output_filename": new_output_filename}
    func_alias(**kwargs)
    gc.collect()

    num_subsets = len(rng_seeds)
    first_or_last_call = True if (subset_idx == num_subsets-1) else False
    
    kwargs = {"output_param_subset": output_param_subset,
              "subset_idx": subset_idx,
              "first_or_last_call": first_or_last_call}
    _remove_temp_files(**kwargs)
    
    unformatted_msg = ("Finished generating the potential slices for the "
                       "frozen phonon configuration subset #{}.\n\n\n")
    msg = unformatted_msg.format(subset_idx)
    print()
    print(msg)

    return None



def _update_pyprismatic_sim_obj_for_next_prismatic_sim(pyprismatic_sim_obj,
                                                       sample_specification,
                                                       output_param_subset,
                                                       subset_idx,
                                                       defocus_idx,
                                                       defocii,
                                                       rng_seeds):
    pyprismatic_sim_obj.numFP = \
        _num_frozen_phonon_configs_in_subset(sample_specification, subset_idx)

    pyprismatic_sim_obj.probeDefocus = defocii[defocus_idx]

    pyprismatic_sim_obj.randomSeed = rng_seeds[subset_idx]

    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "sample_specification": sample_specification,
              "output_param_subset": output_param_subset,
              "subset_idx": subset_idx,
              "defocus_idx": defocus_idx,
              "defocii": defocii}
    _update_save_and_import_params_for_potential_slices(**kwargs)
    _update_save_and_import_params_for_S_matrices(**kwargs)

    return None



def _update_save_and_import_params_for_potential_slices(pyprismatic_sim_obj,
                                                        sample_specification,
                                                        output_param_subset,
                                                        subset_idx,
                                                        defocus_idx,
                                                        defocii):
    save_potential_slices = output_param_subset["save_potential_slices"]

    kwargs = \
        {"sample_specification": sample_specification, "subset_idx": subset_idx}
    precalculated_potential_slice_subset_filename = \
        _precalculated_potential_slice_subset_filename(**kwargs)
    
    if precalculated_potential_slice_subset_filename is None:
        pyprismatic_sim_obj.importFile = ""
        pyprismatic_sim_obj.importPotential = False
    else:
        filename = precalculated_potential_slice_subset_filename
        pyprismatic_sim_obj.importFile = filename
        pyprismatic_sim_obj.importPotential = True
    
    if len(defocii) == 1:
        pyprismatic_sim_obj.savePotentialSlices = save_potential_slices
    else:
        if pyprismatic_sim_obj.algorithm in ("hrtem", "multislice"):
            if defocus_idx == 0:
                pyprismatic_sim_obj.savePotentialSlices = True
            else:
                pyprismatic_sim_obj.savePotentialSlices = False
                pyprismatic_sim_obj.importPotential = True
                output_dirname = output_param_subset["output_dirname"]
                file_to_import = output_dirname + "/potential_slices_temp.h5"
                pyprismatic_sim_obj.importFile = file_to_import
        else:
            if defocus_idx == 0:
                pyprismatic_sim_obj.savePotentialSlices = save_potential_slices
            else:
                pyprismatic_sim_obj.savePotentialSlices = False
                pyprismatic_sim_obj.importPotential = False
                pyprismatic_sim_obj.importFile = ""

    return None



def _precalculated_potential_slice_subset_filename(sample_specification,
                                                   subset_idx):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)
    
    if "atomic_coords_filename" in sample_specification_core_attrs:
        precalculated_potential_slice_subset_filename = None
    else:
        core_attr_name_1 = "interpolation_factors"
        core_attr_name_2 = "filenames"
        
        if core_attr_name_1 in sample_specification_core_attrs:
            filenames = \
                sample_specification_core_attrs[core_attr_name_2]
            precalculated_potential_slice_subset_filename = \
                filenames[subset_idx]
        else:
            precalculated_potential_slice_subset_filename = \
                None

    return precalculated_potential_slice_subset_filename



def _update_save_and_import_params_for_S_matrices(pyprismatic_sim_obj,
                                                  sample_specification,
                                                  output_param_subset,
                                                  subset_idx,
                                                  defocus_idx,
                                                  defocii):
    if pyprismatic_sim_obj.algorithm != "prism":
        pyprismatic_sim_obj.importSMatrix = False
        pyprismatic_sim_obj.saveSMatrix = False
        return None
        
    save_S_matrices = output_param_subset["save_S_matrices"]

    kwargs = \
        {"sample_specification": sample_specification, "subset_idx": subset_idx}
    precalculated_S_matrix_subset_filename = \
        _precalculated_S_matrix_subset_filename(**kwargs)
    
    if precalculated_S_matrix_subset_filename is None:
        pyprismatic_sim_obj.importSMatrix = False
    else:
        pyprismatic_sim_obj.importFile = precalculated_S_matrix_subset_filename
        pyprismatic_sim_obj.importSMatrix = True
    
    if len(defocii) == 1:
        pyprismatic_sim_obj.saveSMatrix = save_S_matrices
    else:
        if defocus_idx == 0:
            pyprismatic_sim_obj.saveSMatrix = True
        else:
            pyprismatic_sim_obj.saveSMatrix = False
            pyprismatic_sim_obj.importSMatrix = True
            output_dirname = output_param_subset["output_dirname"]
            file_to_import = output_dirname + "/S_matrices_temp.h5"
            pyprismatic_sim_obj.importFile = file_to_import

    return None



def _precalculated_S_matrix_subset_filename(sample_specification, subset_idx):
    sample_specification_core_attrs = \
        sample_specification.get_core_attrs(deep_copy=False)
    
    if "atomic_coords_filename" in sample_specification_core_attrs:
        precalculated_S_matrix_subset_filename = None
    else:
        core_attr_name_1 = "interpolation_factors"
        core_attr_name_2 = "filenames"

        if core_attr_name_1 in sample_specification_core_attrs:
            precalculated_S_matrix_subset_filename = \
                None
        else:
            filenames = \
                sample_specification_core_attrs[core_attr_name_2]
            precalculated_S_matrix_subset_filename = \
                filenames[subset_idx]

    return precalculated_S_matrix_subset_filename



def _call_pyprismatic_sim_obj_go(pyprismatic_sim_obj):
    pyprismatic_sim_obj.go()

    return None



def _move_sample_specification_output_to_separate_file(
        output_dirname, sample_specification_type, new_output_filename):
    new_output_filename_candidate = \
        (output_dirname + "/potential_slices_temp.h5"
         if (sample_specification_type == "potential_slice_subset")
         else output_dirname + "/S_matrices_temp.h5")    
    new_output_filename += \
        (new_output_filename == "")*new_output_filename_candidate

    old_output_filename = output_dirname + "/prismatic_output.h5"

    kwargs = {"old_output_filename": old_output_filename,
              "new_output_filename": new_output_filename,
              "sample_specification_type": sample_specification_type}
    _move_sample_specification_output_metadata_to_separate_file(**kwargs)
    _move_sample_specification_output_data_to_separate_file(**kwargs)
    
    return None



def _move_sample_specification_output_metadata_to_separate_file(
        old_output_filename, new_output_filename, sample_specification_type):
    hdf5_attr_names = (("c", "t", "px", "py", "s")
                       if ("potential_slice" in sample_specification_type)
                       else ("c", "t", "px", "py", "fx", "fy"))

    path_in_old_output_file = ("/4DSTEM_simulation/metadata/metadata_0/original"
                               "/simulation_parameters")

    path_in_new_output_file = path_in_old_output_file

    kwargs = {"filename": new_output_filename,
              "path_in_file": path_in_new_output_file}
    group_id = h5pywrappers.obj.ID(**kwargs)

    kwargs = {"group": None, "group_id": group_id, "write_mode": "a"}
    h5pywrappers.group.save(**kwargs)

    for hdf5_attr_name in hdf5_attr_names:
        kwargs = {"filename": old_output_filename,
                  "path_in_file": path_in_old_output_file}
        obj_id = h5pywrappers.obj.ID(**kwargs)

        kwargs = {"obj_id": obj_id, "attr_name": hdf5_attr_name}
        attr_id = h5pywrappers.attr.ID(**kwargs)
        
        hdf5_attr = h5pywrappers.attr.load(attr_id)

        kwargs = {"filename": new_output_filename,
                  "path_in_file": path_in_new_output_file}
        obj_id = h5pywrappers.obj.ID(**kwargs)

        kwargs = {"obj_id": obj_id, "attr_name": hdf5_attr_name}
        attr_id = h5pywrappers.attr.ID(**kwargs)

        kwargs = {"attr": hdf5_attr, "attr_id": attr_id, "write_mode": "a"}
        h5pywrappers.attr.save(**kwargs)

    return None



def _move_sample_specification_output_data_to_separate_file(
        old_output_filename, new_output_filename, sample_specification_type):
    partial_key = ("ppotential"
                   if ("potential_slice" in sample_specification_type)
                   else "smatrix")
        
    path_in_old_output_file = "/4DSTEM_simulation/data/realslices"

    kwargs = {"filename": old_output_filename,
              "path_in_file": path_in_old_output_file}
    group_id = h5pywrappers.obj.ID(**kwargs)

    kwargs = {"group_id": group_id, "read_only": False}
    group_in_old_output_file = h5pywrappers.group.load(**kwargs)

    for key in group_in_old_output_file.keys():
        if partial_key in key:
            dataset = group_in_old_output_file[key]["data"]
            path_in_new_output_file = (path_in_old_output_file
                                       + "/" + key + "/data")

            kwargs = {"filename": new_output_filename,
                      "path_in_file": path_in_new_output_file}
            dataset_id = h5pywrappers.obj.ID(**kwargs)

            kwargs = {"dataset": dataset,
                      "dataset_id": dataset_id,
                      "write_mode": "a"}
            h5pywrappers.dataset.save(**kwargs)
            
            del group_in_old_output_file[key]

    group_in_old_output_file.file.close()

    return None



_default_output_dirname = "S_matrix_generator_output"



def generate_S_matrices(sample_specification,
                        probe_model_params=\
                        _default_probe_model_params,
                        output_dirname=\
                        _default_output_dirname,
                        max_data_size=\
                        _default_max_data_size,
                        worker_params=\
                        _default_worker_params,
                        skip_validation_and_conversion=\
                        _default_skip_validation_and_conversion):
    r"""Generate the :math:`S`-matrices for a given sample model and probe 
    model.

    See the documentation for the subpackage :mod:`prismatique.stem` for a
    discussion on :math:`S`-matrices.

    For each frozen phonon configuration subset, this Python function generates
    the corresponding :math:`S`-matrices and saves them to a set of files. See
    the documentation for the subpackage :mod:`prismatique.stem` for a
    discussion on :math:`S`-matrices, and the class
    :class:`prismatique.thermal.Params` for a discussion on frozen phonon
    configuration subsets.

    For each frozen phonon configuration subset, the corresponding
    :math:`S`-matrix data is written to a file with the basename
    ``"S_matrices_of_subset_"+str(i)+".h5"``, with ``i`` being the subset
    index.

    In addition to the :math:`S`-matrix files, three other files are generated
    as well: The first is a JSON file with the basename
    ``"sample_specification.json"`` which contains, in a serialized format, the
    simulation parameters related to the modelling of the sample; the second
    file is a JSON file with the basename ``"worker_params.json"`` which
    contains, in a serialized format, the simulation parameters related to GPU
    and CPU workers; and the third file is a JSON file with the basename
    ``"probe_model_params.json"`` which contains, in a serialized format, the
    simulation parameters related to the modelling of the probe.

    Parameters
    ----------
    sample_specification : :class:`prismatique.sample.ModelParams` | :class:`prismatique.sample.PotentialSliceSubsetIDs`
        The simulation parameters specifying the sample model. 

        If ``sample_specification`` is of the type
        :class:`prismatique.sample.ModelParams`, then ``sample_specifications``
        specifies sample model parameters that are used to construct the model
        from scratch, i.e. the potential slices for each frozen phonon
        configuration subset are calculated from said model parameters. See the
        documentation for the classes :class:`prismatique.discretization.Params`
        and :class:`prismatique.thermal.Params` for discussions on potential
        slices and frozen phonon configuration subsets respectively.

        Otherwise, if ``sample_specification`` is an instance of the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs` then
        ``sample_specification`` specifies a set of files, where each file
        stores the pre-calculated potential slices for a frozen phonon
        configuration subset. See the documentation for the aforementioned
        class for a further discussion on specifying pre-calculated
        potential slices. 
    probe_model_params : :class:`embeam.stem.probe.ModelParams` | `None`, optional
        The model parameters of the probe. See the documentation for the class
        :class:`embeam.stem.probe.ModelParams` for a discussion on said
        parameters. If ``probe_model_params`` is set to ``None`` [i.e. the
        default value], then the parameter will be reassigned to the value
        ``embeam.stem.probe.ModelParams()``.

        Note that of parameters stored in ``probe_model_params``, only the
        following are used:

        - probe_model_params

          * convergence_semiangle

          * gun_model_params

            + mean_beam_energy

    output_dirname : `str`, optional
        The relative or absolute path to the directory in which all output files
        are to be saved. If the directory doesn't exist upon saving the output 
        files, it will be created if possible.
    max_data_size : `int`, optional
        The data size limit, in bytes, of the :math:`S`-matrices to be
        generated. If the :math:`S`-matrices to be generated would require a
        data size larger than the aforementioned limit, then an exception is
        raised and the :math:`S`-matrices are not generated. Note that data size
        due to HDF5 file overhead and metadata are not taken into account.
    worker_params : :class:`prismatique.worker.Params` | `None`, optional
        The simulation parameters related to GPU and CPU workers. See the
        documentation for the class :class:`prismatique.worker.Params` for a
        discussion on said parameters. If ``worker_params`` is set to `None`
        [i.e. the default value], then the aforementioned simulation parameters 
        are set to default values.
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

    """
    params = locals()

    func_alias = _check_and_convert_skip_validation_and_conversion
    skip_validation_and_conversion = func_alias(params)

    if (skip_validation_and_conversion == False):
        params["accepted_types"] = (ModelParams,
                                    PotentialSliceSubsetIDs)
        global_symbol_table = globals()
        for param_name in params:
            if ((param_name in "skip_validation_and_conversion")
                or (param_name in "accepted_types")):
                continue
            func_name = "_check_and_convert_" + param_name
            func_alias = global_symbol_table[func_name]
            params[param_name] = func_alias(params)
        del params["accepted_types"]

    kwargs = params
    del kwargs["skip_validation_and_conversion"]
    _generate_S_matrices(**kwargs)
    
    return None



def _generate_S_matrices(sample_specification,
                         probe_model_params,
                         output_dirname,
                         worker_params,
                         max_data_size):
    kwargs = {"sample_specification": sample_specification}
    rng_seeds = _generate_rng_seeds_from_sample_specification(**kwargs)

    output_param_subset = {"output_dirname": output_dirname,
                           "save_potential_slices": False,
                           "save_S_matrices": True}

    kwargs = {"output_param_subset": output_param_subset,
              "subset_idx": 0,
              "first_or_last_call": True}
    _remove_temp_files(**kwargs)

    kwargs = \
        {"sample_specification": sample_specification,
         "probe_model_params": probe_model_params,
         "output_dirname": output_dirname,
         "worker_params": worker_params}
    pyprismatic_sim_obj = \
        _initialize_pyprismatic_sim_obj_for_S_matrix_generation(**kwargs)

    for subset_idx, rng_seed in enumerate(rng_seeds):
        kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
                  "sample_specification": sample_specification,
                  "output_param_subset": output_param_subset,
                  "subset_idx": subset_idx,
                  "rng_seeds": rng_seeds}
        _generate_S_matrix_subset(**kwargs)

    kwargs = {"filename": output_dirname + "/sample_specification.json",
              "overwrite": True}
    sample_specification.dump(**kwargs)
    
    kwargs["filename"] = output_dirname + "/probe_model_params.json"
    probe_model_params.dump(**kwargs)
    
    kwargs["filename"] = output_dirname + "/worker_params.json"
    worker_params.dump(**kwargs)
    
    return None



def _initialize_pyprismatic_sim_obj_for_S_matrix_generation(
        sample_specification,
        probe_model_params,
        output_dirname,
        worker_params):
    pyprismatic_sim_obj = pyprismatic.Metadata()

    _set_pyprismatic_sim_obj_attrs_to_default_values(pyprismatic_sim_obj,
                                                     output_dirname)

    _unpack_sample_specification_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                          output_dirname,
                                                          sample_specification)

    _unpack_probe_model_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                        probe_model_params,
                                                        output_dirname)

    pyprismatic_sim_obj.filenameOutput = output_dirname + "/prismatic_output.h5"
    
    _unpack_worker_params_into_pyprismatic_sim_obj(pyprismatic_sim_obj,
                                                   worker_params)

    pyprismatic_sim_obj.algorithm = "prism"
    
    return pyprismatic_sim_obj



def _generate_S_matrix_subset(pyprismatic_sim_obj,
                              sample_specification,
                              output_param_subset,
                              subset_idx,
                              rng_seeds):
    output_dirname = \
        output_param_subset["output_dirname"]
    new_output_filename = \
        output_dirname + "/S_matrices_of_subset_{}.h5".format(subset_idx)
        
    kwargs = {"pyprismatic_sim_obj": pyprismatic_sim_obj,
              "sample_specification": sample_specification,
              "output_param_subset": output_param_subset,
              "subset_idx": subset_idx,
              "defocus_idx": 0,
              "defocii": (0,),
              "rng_seeds": rng_seeds}
    _update_pyprismatic_sim_obj_for_next_prismatic_sim(**kwargs)

    unformatted_msg = ("Generating the S-matrices for the frozen phonon "
                       "configuration subset #{}...\n")
    msg = unformatted_msg.format(subset_idx)
    print(msg)

    # Run prismatic simulation.
    _call_pyprismatic_sim_obj_go(pyprismatic_sim_obj)
    gc.collect()
    
    func_alias = _move_sample_specification_output_to_separate_file
    kwargs = {"output_dirname": output_dirname,
              "sample_specification_type": "S_matrix_subset",
              "new_output_filename": new_output_filename}
    func_alias(**kwargs)
    gc.collect()

    num_subsets = len(rng_seeds)
    first_or_last_call = True if (subset_idx == num_subsets-1) else False
    
    kwargs = {"output_param_subset": output_param_subset,
              "subset_idx": subset_idx,
              "first_or_last_call": first_or_last_call}
    _remove_temp_files(**kwargs)
    
    unformatted_msg = ("Finished generating the S-matrices for the frozen "
                       "phonon configuration subset #{}.\n\n\n")
    msg = unformatted_msg.format(subset_idx)
    print()
    print(msg)

    return None



###########################
## Define error messages ##
###########################

_check_and_convert_unit_cell_tiling_err_msg_1 = \
    ("The object ``unit_cell_tiling`` must be a triplet of positive integers.")

_check_atomic_coords_file_format_err_msg_1 = \
    ("No file was found with the path ``'{}'``.")
_check_atomic_coords_file_format_err_msg_2 = \
    ("The last line in the atomic coordinates file ``'{}'`` should be '-1'.")
_check_atomic_coords_file_format_err_msg_3 = \
    ("The second line in the atomic coordinates file ``'{}'`` should be of the "
     "form 'a b c', where 'a', 'b', and 'c' are positive floating-point "
     "numbers.")
_check_atomic_coords_file_format_err_msg_4 = \
    ("Line #{} of the atomic coordinates file ``'{}'`` is not formatted "
     "correctly: it should be of the form 'Z x y z occ u_x_rms', where 'Z' is "
     "a positive integer, 'x', 'y', and 'z' are floating-point numbers, 'occ' "
     "is a floating-point number in the interval [0, 1]`, and 'u_x_rms' is a "
     "non-negative floating-point number.")

_check_and_convert_filenames_err_msg_1 = \
    ("The object ``filenames`` must satisfy ``len(filenames)>0``.")

_check_and_convert_max_num_frozen_phonon_configs_per_subset_err_msg_1 = \
    ("The object ``max_num_frozen_phonon_configs_per_subset`` must be an "
     "integer greater than or equal to ``2``.")

_check_and_convert_potential_slice_subset_filenames_err_msg_1 = \
    ("The object ``potential_slice_subset_filenames`` must satisfy "
     "``len(potential_slice_subset_filenames)>0``.")

_check_and_convert_S_matrix_subset_filenames_err_msg_1 = \
    ("The object ``S_matrix_subset_filenames`` must satisfy "
     "``len(S_matrix_subset_filenames)>0``.")

_check_that_hdf5_files_exist_err_msg_1 = \
    ("The object ``sample_specification.core_attrs['{}'][{}]``, which "
     "stores the string ``'{}'``, specifies an invalid HDF5 filename: see "
     "traceback for details.")

_load_hdf5_attr_err_msg_1 = \
    ("The HDF5 file ``'{}'``, specified by the object "
     "``sample_specification.core_attrs['{}'][{}]``, must store an HDF5 group "
     "at the HDF5 path "
     "'4DSTEM_simulation/metadata/metadata_0/original/simulation_parameters' "
     "with an HDF5 attribute ``'{}'``{}")

_check_hdf5_attr_err_msg_1 = \
    ("Each string in "
     "``sample_specification.core_attrs['potential_slice_subset_filenames']`` "
     "and "
     "``sample_specification.core_attrs['S_matrix_subset_filenames']`` must "
     "specify the filename of an HDF5 file, wherein there is an HDF5 group at "
     "the HDF5 path "
     "'4DSTEM_simulation/metadata/metadata_0/original/simulation_parameters' "
     "with an HDF5 attribute ``'{}'`` that is equal to that for every other "
     "string in the aforementioned string sequences.")
_check_hdf5_attr_err_msg_2 = \
    ("Each string in ``sample_specification.core_attrs['{}']`` must specify "
     "the filename of an HDF5 file, wherein there is an HDF5 group at the HDF5 "
     "path "
     "'4DSTEM_simulation/metadata/metadata_0/original/simulation_parameters' "
     "with an HDF5 attribute ``'{}'`` that is equal to that for every other "
     "string in ``sample_specification.core_attrs['{}']``.")

_check_datasets_in_hdf5_files_err_msg_1 = \
    ("The HDF5 file ``'{}'``, specified by the object "
     "``sample_specification.core_attrs['{}'][{}]``, must store a 3D HDF5 "
     "dataset of 32-bit floating-point numbers at the HDF5 path "
     "'4DSTEM_simulation/data/realslices/ppotential_fp0000/data' satisfying "
     "``np.abs(dataset.shape[0]*sim_params.attrs['px']"
     " - sim_params.attrs['c'][0]*sim_params.attrs['t'][0]) < epsilon``, "
     "``np.abs(dataset.shape[1]*sim_params.attrs['py']"
     " - sim_params.attrs['c'][1]*sim_params.attrs['t'][1]) < epsilon``, and "
     "``np.abs(dataset.shape[2]*sim_params.attrs['s']"
     " - sim_params.attrs['c'][2]*sim_params.attrs['t'][2]) < epsilon``, where "
     "``dataset`` is the HDF5 dataset, ``epsilon`` is the single-precision "
     "machine epsilon, and ``sim_params`` is an HDF5 group located at the HDF5 "
     "path "
     "'4DSTEM_simulation/metadata/metadata_0/original/simulation_parameters' "
     "in the same HDF5 file.")
_check_datasets_in_hdf5_files_err_msg_2 = \
    ("The HDF5 file ``'{}'``, specified by the object "
     "``sample_specification.core_attrs['{}'][{}]``, must store a 3D HDF5 "
     "dataset of 64-bit complex numbers at the HDF5 path "
     "'4DSTEM_simulation/data/realslices/smatrix_fp0000/data' satisfying "
     "``np.abs(2*dataset.shape[0]*sim_params.attrs['fx']*sim_params.attrs['px']"
     " - sim_params.attrs['c'][0]*sim_params.attrs['t'][0]) < epsilon``, and "
     "``np.abs(dataset.shape[1]*sim_params.attrs['fy']*sim_params.attrs['py']"
     " - sim_params.attrs['c'][1]*sim_params.attrs['t'][1]) < epsilon``, where "
     "``dataset`` is the HDF5 dataset, ``epsilon`` is the single-precision "
     "machine epsilon, and ``sim_params`` is an HDF5 group located at the HDF5 "
     "path "
     "'4DSTEM_simulation/metadata/metadata_0/original/simulation_parameters' "
     "in the same HDF5 file.")

_ppotential_fp_indices_err_msg_1 = \
    ("The HDF5 file ``'{}'``, specified by the object "
     "``sample_specification.core_attrs['{}'][{}]``, must store an HDF5 group "
     "at the HDF5 path '4DSTEM_simulation/data/realslices'.")

_smatrix_fp_indices_err_msg_1 = \
    _ppotential_fp_indices_err_msg_1

_check_that_probe_model_is_compatible_with_S_matrices_err_msg_1 = \
    ("The HDF5 file ``'{}'``, specified by the object "
     "``sample_specification.core_attrs['{}'][{}]``, stores a 3D HDF5 dataset "
     "at the HDF5 path "
     "'4DSTEM_simulation/data/realslices/smatrix_fp0000/data' with a shape "
     "that is incompatible with the probe model parameters specified the "
     "object ``probe_model_params``, namely the convergence semiangle and the "
     "mean beam energy: make sure that the pre-calculated S-matrices were "
     "calculated using the same convergence semiangle and mean beam energy "
     "as that specified by ``probe_model_params``.")

_check_and_convert_subset_idx_err_msg_1 = \
    ("The object ``subset_idx`` must be an integer satisfying "
     "``0<=subset_idx<prismatique.sample.num_frozen_phonon_config_subsets("
     "sample_specification)``.")

_check_data_size_err_msg_1 = \
    ("The data size of the set of {} to be generated is {} bytes, exceeding "
     "the maximum allowed size of {} bytes specified by the object "
     "``max_data_size``.")
