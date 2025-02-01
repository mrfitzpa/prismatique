r"""For specifying simulation parameters related to the modelling of HRTEM
systems.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For performing deep copies of objects.
import copy



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
__all__ = ["Params"]



def _check_and_convert_sample_specification(ctor_params):
    sample_specification = copy.deepcopy(ctor_params["sample_specification"])
    
    accepted_types = (prismatique.sample.ModelParams,
                      prismatique.sample.PotentialSliceSubsetIDs)
    kwargs = {"obj": sample_specification,
              "obj_name": "sample_specification",
              "accepted_types": accepted_types}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)
    
    return sample_specification



def _pre_serialize_sample_specification(sample_specification):
    if "thermal_params" in sample_specification.core_attrs:
        pre_serialize_sample_specification = \
            prismatique.sample._pre_serialize_sample_model_params
    else:
        pre_serialize_sample_specification = \
            prismatique.sample._pre_serialize_potential_slice_subset_ids
        
    serializable_rep = pre_serialize_sample_specification(sample_specification)

    return serializable_rep



def _de_pre_serialize_sample_specification(serializable_rep):
    if "thermal_params" in serializable_rep:
        de_pre_serialize_sample_specification = \
            prismatique.sample._de_pre_serialize_sample_model_params
    else:
        de_pre_serialize_sample_specification = \
            prismatique.sample._de_pre_serialize_potential_slice_subset_ids
        
    sample_specification = \
        de_pre_serialize_sample_specification(serializable_rep)

    return sample_specification



def _check_and_convert_gun_model_params(ctor_params):
    check_and_convert_gun_model_params = \
        embeam.gun._check_and_convert_gun_model_params
    gun_model_params = \
        check_and_convert_gun_model_params(ctor_params)
    
    return gun_model_params



def _pre_serialize_gun_model_params(gun_model_params):
    pre_serialize_gun_model_params = \
        embeam.gun._pre_serialize_gun_model_params
    serializable_rep = \
        pre_serialize_gun_model_params(gun_model_params)

    return serializable_rep



def _de_pre_serialize_gun_model_params(serializable_rep):
    de_pre_serialize_gun_model_params = \
        embeam.gun._de_pre_serialize_gun_model_params
    gun_model_params = \
        de_pre_serialize_gun_model_params(serializable_rep)

    return gun_model_params



def _check_and_convert_lens_model_params(ctor_params):
    check_and_convert_lens_model_params = \
        embeam.lens._check_and_convert_lens_model_params
    lens_model_params = \
        check_and_convert_lens_model_params(ctor_params)
    
    return lens_model_params



def _pre_serialize_lens_model_params(lens_model_params):
    pre_serialize_lens_model_params = \
        embeam.lens._pre_serialize_lens_model_params
    serializable_rep = \
        pre_serialize_lens_model_params(lens_model_params)

    return serializable_rep



def _de_pre_serialize_lens_model_params(serializable_rep):
    de_pre_serialize_lens_model_params = \
        embeam.lens._de_pre_serialize_lens_model_params
    lens_model_params = \
        de_pre_serialize_lens_model_params(serializable_rep)

    return lens_model_params



def _check_and_convert_tilt_params(ctor_params):
    check_and_convert_tilt_params = \
        prismatique.tilt._check_and_convert_tilt_params
    tilt_params = \
        check_and_convert_tilt_params(ctor_params)
    
    return tilt_params



def _pre_serialize_tilt_params(tilt_params):
    pre_serialize_tilt_params = \
        prismatique.tilt._pre_serialize_tilt_params
    serializable_rep = \
        pre_serialize_tilt_params(tilt_params)

    return serializable_rep



def _de_pre_serialize_tilt_params(serializable_rep):
    de_pre_serialize_tilt_params = \
        prismatique.tilt._de_pre_serialize_tilt_params
    tilt_params = \
        de_pre_serialize_tilt_params(serializable_rep)

    return tilt_params



def _check_and_convert_objective_aperture_params(ctor_params):
    check_and_convert_objective_aperture_params = \
        prismatique.aperture._check_and_convert_objective_aperture_params
    objective_aperture_params = \
        check_and_convert_objective_aperture_params(ctor_params)
    
    return objective_aperture_params



def _pre_serialize_objective_aperture_params(objective_aperture_params):
    pre_serialize_objective_aperture_params = \
        prismatique.aperture._pre_serialize_objective_aperture_params
    serializable_rep = \
        pre_serialize_objective_aperture_params(objective_aperture_params)

    return serializable_rep



def _de_pre_serialize_objective_aperture_params(serializable_rep):
    de_pre_serialize_objective_aperture_params = \
        prismatique.aperture._de_pre_serialize_objective_aperture_params
    objective_aperture_params = \
        de_pre_serialize_objective_aperture_params(serializable_rep)

    return objective_aperture_params



def _check_and_convert_defocal_offset_supersampling(ctor_params):
    kwargs = {"obj": ctor_params["defocal_offset_supersampling"],
              "obj_name": "defocal_offset_supersampling"}
    defocal_offset_supersampling = czekitout.convert.to_positive_int(**kwargs)

    return defocal_offset_supersampling



def _pre_serialize_defocal_offset_supersampling(defocal_offset_supersampling):
    serializable_rep = defocal_offset_supersampling
    
    return serializable_rep



def _de_pre_serialize_defocal_offset_supersampling(serializable_rep):
    defocal_offset_supersampling = serializable_rep

    return defocal_offset_supersampling



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

    Attributes
    ----------
    core_attrs : `dict`, read-only
        A `dict` representation of the core attributes: each `dict` key is a
        `str` representing the name of a core attribute, and the corresponding
        `dict` value is the object to which said core attribute is set. The core
        attributes are the same as the construction parameters, except that 
        their values might have been updated since construction.

    """
    _validation_and_conversion_funcs = \
        {"sample_specification": _check_and_convert_sample_specification,
         "gun_model_params": _check_and_convert_gun_model_params,
         "lens_model_params": _check_and_convert_lens_model_params,
         "tilt_params": _check_and_convert_tilt_params,
         "objective_aperture_params": \
         _check_and_convert_objective_aperture_params,
         "defocal_offset_supersampling": \
         _check_and_convert_defocal_offset_supersampling}

    _pre_serialization_funcs = \
        {"sample_specification": _pre_serialize_sample_specification,
         "gun_model_params": _pre_serialize_gun_model_params,
         "lens_model_params": _pre_serialize_lens_model_params,
         "tilt_params": _pre_serialize_tilt_params,
         "objective_aperture_params": _pre_serialize_objective_aperture_params,
         "defocal_offset_supersampling": \
         _pre_serialize_defocal_offset_supersampling}

    _de_pre_serialization_funcs = \
        {"sample_specification": _de_pre_serialize_sample_specification,
         "gun_model_params": _de_pre_serialize_gun_model_params,
         "lens_model_params": _de_pre_serialize_lens_model_params,
         "tilt_params": _de_pre_serialize_tilt_params,
         "objective_aperture_params": \
         _de_pre_serialize_objective_aperture_params,
         "defocal_offset_supersampling": \
         _de_pre_serialize_defocal_offset_supersampling}
    
    def __init__(self,
                 sample_specification,
                 gun_model_params=None,
                 lens_model_params=None,
                 tilt_params=None,
                 objective_aperture_params=None,
                 defocal_offset_supersampling=1):
        ctor_params = \
            {"sample_specification": sample_specification,
             "gun_model_params": gun_model_params,
             "lens_model_params": lens_model_params,
             "tilt_params": tilt_params,
             "objective_aperture_params": objective_aperture_params,
             "defocal_offset_supersampling": defocal_offset_supersampling}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_hrtem_system_model_params(ctor_params):
    hrtem_system_model_params = \
        copy.deepcopy(ctor_params["hrtem_system_model_params"])    
    kwargs = {"obj": hrtem_system_model_params,
              "obj_name": "hrtem_system_model_params",
              "accepted_types": (ModelParams,)}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return hrtem_system_model_params



def _pre_serialize_hrtem_system_model_params(hrtem_system_model_params):
    serializable_rep = hrtem_system_model_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_hrtem_system_model_params(serializable_rep):
    hrtem_system_model_params = ModelParams.de_pre_serialize(serializable_rep)

    return hrtem_system_model_params
    


###########################
## Define error messages ##
###########################