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
r"""For specifying simulation parameters related to the modelling of HRTEM
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
# classes :class:`embeam.gun.ModelParams`, and :class:`embeam.lens.ModelParams`.
import embeam.gun
import embeam.lens



# For validating, pre-serializing, and de-pre-serializing instances of the
# classes :class:`prismatique.sample.ModelParams`,
# :class:`prismatique.sample.PotentialSliceSubsetIDs`,
# :class:`prismatique.tilt.Params`, and
# :class:`prismatique.aperture.Params`.
import prismatique.sample
import prismatique.tilt
import prismatique.aperture



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params"]



def _check_and_convert_sample_specification(params):
    params["accepted_types"] = (prismatique.sample.ModelParams,
                                prismatique.sample.PotentialSliceSubsetIDs)

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



def _check_and_convert_gun_model_params(params):
    obj_name = "gun_model_params"

    module_alias = embeam.stem.probe
    cls_alias = module_alias.ModelParams
    func_alias = cls_alias.get_validation_and_conversion_funcs()[obj_name]
    gun_model_params = func_alias(params)
    
    return gun_model_params



def _pre_serialize_gun_model_params(gun_model_params):
    obj_to_pre_serialize = gun_model_params

    obj_name = "gun_model_params"

    module_alias = embeam.stem.probe
    cls_alias = module_alias.ModelParams
    func_alias = cls_alias.get_pre_serialization_funcs()[obj_name]
    serializable_rep = func_alias(obj_to_pre_serialize)
    
    return serializable_rep



def _de_pre_serialize_gun_model_params(serializable_rep):
    obj_name = "gun_model_params"

    module_alias = embeam.stem.probe
    cls_alias = module_alias.ModelParams
    func_alias = cls_alias.get_de_pre_serialization_funcs()[obj_name]
    gun_model_params = func_alias(serializable_rep)

    return gun_model_params



def _check_and_convert_lens_model_params(params):
    obj_name = "lens_model_params"

    module_alias = embeam.stem.probe
    cls_alias = module_alias.ModelParams
    func_alias = cls_alias.get_validation_and_conversion_funcs()[obj_name]
    lens_model_params = func_alias(params)
    
    return lens_model_params



def _pre_serialize_lens_model_params(lens_model_params):
    obj_to_pre_serialize = lens_model_params

    obj_name = "lens_model_params"

    module_alias = embeam.stem.probe
    cls_alias = module_alias.ModelParams
    func_alias = cls_alias.get_pre_serialization_funcs()[obj_name]
    serializable_rep = func_alias(obj_to_pre_serialize)
    
    return serializable_rep



def _de_pre_serialize_lens_model_params(serializable_rep):
    obj_name = "lens_model_params"

    module_alias = embeam.stem.probe
    cls_alias = module_alias.ModelParams
    func_alias = cls_alias.get_de_pre_serialization_funcs()[obj_name]
    lens_model_params = func_alias(serializable_rep)

    return lens_model_params



def _check_and_convert_tilt_params(params):
    module_alias = prismatique.tilt
    func_alias = module_alias._check_and_convert_tilt_params
    tilt_params = func_alias(params)
    
    return tilt_params



def _pre_serialize_tilt_params(tilt_params):
    obj_to_pre_serialize = tilt_params
    module_alias = prismatique.tilt
    func_alias = module_alias._pre_serialize_tilt_params
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_tilt_params(serializable_rep):
    module_alias = prismatique.tilt
    func_alias = module_alias._de_pre_serialize_tilt_params
    tilt_params = func_alias(serializable_rep)

    return tilt_params



def _check_and_convert_objective_aperture_params(params):
    module_alias = prismatique.aperture
    func_alias = module_alias._check_and_convert_objective_aperture_params
    objective_aperture_params = func_alias(params)
    
    return objective_aperture_params



def _pre_serialize_objective_aperture_params(objective_aperture_params):
    obj_to_pre_serialize = objective_aperture_params
    module_alias = prismatique.aperture
    func_alias = module_alias._pre_serialize_objective_aperture_params
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_objective_aperture_params(serializable_rep):
    module_alias = prismatique.aperture
    func_alias = module_alias._de_pre_serialize_objective_aperture_params
    objective_aperture_params = func_alias(serializable_rep)

    return objective_aperture_params



def _check_and_convert_defocal_offset_supersampling(params):
    obj_name = "defocal_offset_supersampling"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    defocal_offset_supersampling = czekitout.convert.to_positive_int(**kwargs)

    return defocal_offset_supersampling



def _pre_serialize_defocal_offset_supersampling(defocal_offset_supersampling):
    obj_to_pre_serialize = defocal_offset_supersampling
    serializable_rep = obj_to_pre_serialize
    
    return serializable_rep



def _de_pre_serialize_defocal_offset_supersampling(serializable_rep):
    defocal_offset_supersampling = serializable_rep

    return defocal_offset_supersampling



_module_alias_1 = \
    prismatique.tilt
_module_alias_2 = \
    prismatique.aperture
_default_gun_model_params = \
    None
_default_lens_model_params = \
    None
_default_tilt_params = \
    _module_alias_1._default_tilt_params
_default_objective_aperture_params = \
    _module_alias_2._default_objective_aperture_params
_default_defocal_offset_supersampling = \
    9
_default_skip_validation_and_conversion = \
    _module_alias_2._default_skip_validation_and_conversion



class ModelParams(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related to the modelling of HRTEM systems.

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
    gun_model_params : :class:`embeam.gun.ModelParams` | `None`, optional
        The electron gun model parameters. See the documentation for the class
        :class:`embeam.gun.ModelParams` for a discussion on said parameters. If
        ``gun_model_params`` is set to ``None`` [i.e. the default value], then
        the aforementioned model parameters are set to default values.
    lens_model_params : :class:`embeam.lens.ModelParams` | `None`, optional
        The model parameters of the objective lens. See the documentation for
        the class :class:`embeam.lens.ModelParams` for a discussion on said
        parameters. If ``lens_model_params`` is set to ``None`` [i.e. the
        default value], then the aforementioned model parameters are set to
        default values.
    tilt_params : :class:`prismatique.tilt.Params` | `None`, optional
        The simulation parameters related to the beam tilt series in the HRTEM
        simulation to model a set of spatially coherent HRTEM experiments at
        different beam tilts, or to model a single spatially incoherent HRTEM
        beam. See the documentation for the class
        :class:`prismatique.tilt.Params` for a discussion on said parameters.
        If ``tilt_params`` is set to `None` [i.e. the default value], then the
        aforementioned simulation parameters are set to default values.
    objective_aperture_params : :class:`prismatique.aperture.Params` | `None`, optional
        The simulation parameters related to the objective aperture. See the
        documentation for the class :class:`prismatique.aperture.Params` for a
        discussion on said parameters.  If ``objective_aperture_params`` is set
        to `None` [i.e. the default value], then the aforementioned simulation
        parameters are set to default values.
    defocal_offset_supersampling : `int`, optional
        The number of points :math:`N_f` to use in the Gauss-Hermite quadrature
        scheme used to approximate the integration over the defocal offset
        :math:`\delta_{f}` in
        Eq. :eq:`mixed_state_operator_for_transmitted_electron`. Must be
        positive.
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
                        "gun_model_params",
                        "lens_model_params",
                        "tilt_params",
                        "objective_aperture_params",
                        "defocal_offset_supersampling")
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
                 gun_model_params=\
                 _default_gun_model_params,
                 lens_model_params=\
                 _default_lens_model_params,
                 tilt_params=\
                 _default_tilt_params,
                 objective_aperture_params=\
                 _default_objective_aperture_params,
                 defocal_offset_supersampling=\
                 _default_defocal_offset_supersampling,
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



def _check_and_convert_hrtem_system_model_params(params):
    obj_name = "hrtem_system_model_params"
    obj = params[obj_name]

    accepted_types = (ModelParams,)

    kwargs = {"obj": obj,
              "obj_name": obj_name,
              "accepted_types": accepted_types}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    kwargs = obj.get_core_attrs(deep_copy=False)
    hrtem_system_model_params = accepted_types[0](**kwargs)

    if "worker_params" not in params:
        hrtem_system_model_params_core_attrs = \
            hrtem_system_model_params.get_core_attrs(deep_copy=False)

        sample_specification = \
            hrtem_system_model_params_core_attrs["sample_specification"]

        kwargs = {"params": dict()}
        kwargs["params"]["sample_specification"] = sample_specification
        _check_and_convert_sample_specification(**kwargs)

    return hrtem_system_model_params



def _pre_serialize_hrtem_system_model_params(hrtem_system_model_params):
    obj_to_pre_serialize = hrtem_system_model_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_hrtem_system_model_params(serializable_rep):
    hrtem_system_model_params = ModelParams.de_pre_serialize(serializable_rep)
    
    return hrtem_system_model_params



###########################
## Define error messages ##
###########################
