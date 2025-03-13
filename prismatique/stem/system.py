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
r"""For specifying simulation parameters related to the modelling of STEM
systems.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For validating and converting objects.
import czekitout.check
import czekitout.convert

# For defining classes that support enforced validation, updatability,
# pre-serialization, and de-serialization.
import fancytypes



# For validating, pre-serializing, and de-pre-serializing instances of the
# classes :class:`prismatique.sample.ModelParams`,
# :class:`prismatique.sample.PotentialSliceSubsetIDs`,
# :class:`prismatique.sample.SMatrixSubsetIDs`, and
# :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`.
import prismatique.sample

# For validating instances of the :class:`prismatique.scan.rectangular.Params`,
# and for generating probe positions.
import prismatique.scan



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params"]



def _check_and_convert_sample_specification(params):
    params["accepted_types"] = (prismatique.sample.ModelParams,
                                prismatique.sample.PotentialSliceSubsetIDs,
                                prismatique.sample.SMatrixSubsetIDs)
    
    module_alias = prismatique.sample
    func_alias = module_alias._check_and_convert_sample_specification
    sample_specification = func_alias(params)

    del params["accepted_types"]

    return sample_specification



def _pre_serialize_sample_specification(sample_specification):
    obj_to_pre_serialize = sample_specification
    module_alias = prismatique.sample
    func_alias = module_alias._pre_serialize_sample_specification
    serializable_rep = func_alias(obj_to_pre_serialize)
    
    return serializable_rep



def _de_pre_serialize_sample_specification(serializable_rep):
    module_alias = prismatique.sample
    func_alias = module_alias._de_pre_serialize_sample_specification
    sample_specification = func_alias(serializable_rep)

    return sample_specification



def _check_and_convert_probe_model_params(params):
    module_alias = prismatique.sample
    func_alias = module_alias._check_and_convert_probe_model_params
    probe_model_params = func_alias(params)

    return probe_model_params



def _pre_serialize_probe_model_params(probe_model_params):
    obj_to_pre_serialize = probe_model_params
    module_alias = prismatique.sample
    func_alias = module_alias._pre_serialize_probe_model_params
    serializable_rep = func_alias(obj_to_pre_serialize)
    
    return serializable_rep



def _de_pre_serialize_probe_model_params(serializable_rep):
    module_alias = prismatique.sample
    func_alias = module_alias._de_pre_serialize_probe_model_params
    probe_model_params = func_alias(serializable_rep)

    return probe_model_params



def _check_and_convert_specimen_tilt(params):
    obj_name = "specimen_tilt"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    specimen_tilt = czekitout.convert.to_pair_of_floats(**kwargs)
    
    return specimen_tilt



def _pre_serialize_specimen_tilt(specimen_tilt):
    obj_to_pre_serialize = specimen_tilt
    serializable_rep = specimen_tilt

    return serializable_rep



def _de_pre_serialize_specimen_tilt(serializable_rep):
    specimen_tilt = serializable_rep

    return specimen_tilt



def _check_and_convert_scan_config(params):
    module_alias = prismatique.scan
    func_alias = module_alias._check_and_convert_scan_config
    scan_config = func_alias(params)
    
    return scan_config



def _pre_serialize_scan_config(scan_config):
    obj_to_pre_serialize = scan_config
    if isinstance(obj_to_pre_serialize, str):
        serializable_rep = obj_to_pre_serialize
    elif isinstance(obj_to_pre_serialize, prismatique.scan.rectangular.Params):
        module_alias = prismatique.scan.rectangular
        func_alias = module_alias._pre_serialize_rectangular_scan_params
        kwargs = {"rectangular_scan_params": obj_to_pre_serialize}
        serializable_rep = func_alias(**kwargs)
    else:
        serializable_rep = scan_config

    return serializable_rep



def _de_pre_serialize_scan_config(serializable_rep):
    if isinstance(serializable_rep, str):
        scan_config = serializable_rep
    elif isinstance(serializable_rep, dict):
        module_alias = prismatique.scan.rectangular
        func_alias = module_alias._de_pre_serialize_rectangular_scan_params
        scan_config = func_alias(serializable_rep)
    else:
        func_alias = czekitout.convert.to_real_two_column_numpy_matrix
        kwargs = {"obj": serializable_rep, "obj_name": "serializable_rep"}
        scan_config = func_alias(**kwargs)

    return scan_config



_module_alias_1 = \
    prismatique.sample
_module_alias_2 = \
    prismatique.scan
_default_probe_model_params = \
    _module_alias_1._default_probe_model_params
_default_specimen_tilt = \
    (0, 0)
_default_scan_config = \
    _module_alias_2._default_scan_config
_default_skip_validation_and_conversion = \
    _module_alias_2._default_skip_validation_and_conversion



class ModelParams(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related to the modelling of STEM systems.

    Parameters
    ----------
    sample_specification : :class:`prismatique.sample.ModelParams` | :class:`prismatique.sample.PotentialSliceSubsetIDs` | :class:`prismatique.sample.SMatrixSubsetIDs`
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

        Moreover, note that ``sample_specification`` must be an instance of the
        class :class:`prismatique.sample.ModelParams`, or the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs` if using the
        multislice algorithm for the STEM simulation, or if the potential slices
        to be used in the simulation are to be saved as output, per the output
        parameters.
    probe_model_params : :class:`embeam.stem.probe.ModelParams` | `None`, optional
        The model parameters of the probe. See the documentation for the class
        :class:`embeam.stem.probe.ModelParams` for a discussion on said
        parameters. If ``probe_model_params`` is set to `None` [i.e. the default
        value], then the aforementioned model parameters are set to default
        values.
    specimen_tilt : `array_like` (`float`, shape=(``2``,)), optional
        ``specimen_tilt[0]`` specifies the specimen tilt along the
        :math:`x`-axis in mrads; ``specimen_tilt[1]`` specifies the specimen
        tilt along the :math:`y`-axis in mrads. Note that according to
        Ref. [Cowley1]_, the method used to simulate specimen tilt in
        ``prismatic`` is only valid for small tilts up to about 1 degree.
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
    ctor_param_names = ("sample_specification",
                        "probe_model_params",
                        "specimen_tilt",
                        "scan_config")
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
                 sample_specification,
                 probe_model_params=\
                 _default_probe_model_params,
                 specimen_tilt=\
                 _default_specimen_tilt,
                 scan_config=\
                 _default_scan_config,
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



def _check_and_convert_stem_system_model_params(params):
    obj_name = "stem_system_model_params"
    obj = params[obj_name]

    accepted_types = (ModelParams,)

    kwargs = {"obj": obj,
              "obj_name": obj_name,
              "accepted_types": accepted_types}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    kwargs = obj.get_core_attrs(deep_copy=False)
    stem_system_model_params = accepted_types[0](**kwargs)

    if "worker_params" not in params:
        stem_system_model_params_core_attrs = \
            stem_system_model_params.get_core_attrs(deep_copy=False)

        sample_specification = \
            stem_system_model_params_core_attrs["sample_specification"]
        probe_model_params = \
            stem_system_model_params_core_attrs["probe_model_params"]
        scan_config = \
            stem_system_model_params_core_attrs["scan_config"]

        kwargs = {"params": dict()}
        kwargs["params"]["sample_specification"] = sample_specification
        kwargs["params"]["probe_model_params"] = probe_model_params
        _check_and_convert_probe_model_params(**kwargs)

        kwargs = {"sample_specification": sample_specification,
                  "scan_config": scan_config,
                  "filename": None}
        prismatique.scan._generate_probe_positions(**kwargs)

    return stem_system_model_params



def _pre_serialize_stem_system_model_params(stem_system_model_params):
    obj_to_pre_serialize = stem_system_model_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_stem_system_model_params(serializable_rep):
    stem_system_model_params = ModelParams.de_pre_serialize(serializable_rep)
    
    return stem_system_model_params



###########################
## Define error messages ##
###########################
