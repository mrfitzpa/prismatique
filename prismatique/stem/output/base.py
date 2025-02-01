r"""For specifying the base output parameters for STEM simulations.

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
# class :class:`prismatique.cbed.Params`.
import prismatique.cbed



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



def _check_and_convert_output_dirname(ctor_params):
    kwargs = {"obj": ctor_params["output_dirname"],
              "obj_name": "output_dirname"}
    output_dirname = czekitout.convert.to_str_from_path_like(**kwargs)
    
    return output_dirname



def _pre_serialize_output_dirname(output_dirname):
    serializable_rep = output_dirname

    return serializable_rep



def _de_pre_serialize_output_dirname(serializable_rep):
    output_dirname = serializable_rep

    return output_dirname



def _check_and_convert_max_data_size(ctor_params):
    kwargs = {"obj": ctor_params["max_data_size"],
              "obj_name": "max_data_size"}
    max_data_size = czekitout.convert.to_positive_int(**kwargs)
    
    return max_data_size



def _pre_serialize_max_data_size(max_data_size):
    serializable_rep = max_data_size

    return serializable_rep



def _de_pre_serialize_max_data_size(serializable_rep):
    max_data_size = serializable_rep

    return max_data_size



def _check_and_convert_cbed_params(ctor_params):
    check_and_convert_cbed_params = \
        prismatique.cbed._check_and_convert_cbed_params
    cbed_params = \
        check_and_convert_cbed_params(ctor_params)
    
    return cbed_params



def _pre_serialize_cbed_params(cbed_params):
    pre_serialize_cbed_params = \
        prismatique.cbed._pre_serialize_cbed_params
    serializable_rep = \
        pre_serialize_cbed_params(cbed_params)

    return serializable_rep



def _de_pre_serialize_cbed_params(serializable_rep):
    de_pre_serialize_cbed_params = \
        prismatique.cbed._de_pre_serialize_cbed_params
    cbed_params = \
        de_pre_serialize_cbed_params(serializable_rep)

    return cbed_params



def _check_and_convert_radial_step_size_for_3d_stem(ctor_params):
    kwargs = {"obj": ctor_params["radial_step_size_for_3d_stem"],
              "obj_name": "radial_step_size_for_3d_stem"}
    radial_step_size_for_3d_stem = \
        czekitout.convert.to_nonnegative_float(**kwargs)
    
    return radial_step_size_for_3d_stem



def _pre_serialize_radial_step_size_for_3d_stem(radial_step_size_for_3d_stem):
    serializable_rep = radial_step_size_for_3d_stem

    return serializable_rep



def _de_pre_serialize_radial_step_size_for_3d_stem(serializable_rep):
    radial_step_size_for_3d_stem = serializable_rep

    return radial_step_size_for_3d_stem



def _check_and_convert_radial_range_for_2d_stem(ctor_params):
    kwargs = {"obj": ctor_params["radial_range_for_2d_stem"],
              "obj_name": "radial_range_for_2d_stem"}
    radial_range_for_2d_stem = czekitout.convert.to_pair_of_floats(**kwargs)

    if not (0 <= radial_range_for_2d_stem[0] <= radial_range_for_2d_stem[1]):
        raise TypeError(_check_and_convert_radial_range_for_2d_stem_err_msg_1)
    
    return radial_range_for_2d_stem



def _pre_serialize_radial_range_for_2d_stem(radial_range_for_2d_stem):
    serializable_rep = radial_range_for_2d_stem
    
    return serializable_rep



def _de_pre_serialize_radial_range_for_2d_stem(serializable_rep):
    radial_range_for_2d_stem = serializable_rep
    
    return radial_range_for_2d_stem



def _check_and_convert_save_com(ctor_params):
    kwargs = {"obj": ctor_params["save_com"],
              "obj_name": "save_com"}
    save_com = czekitout.convert.to_bool(**kwargs)
    
    return save_com



def _pre_serialize_save_com(save_com):
    serializable_rep = save_com

    return serializable_rep



def _de_pre_serialize_save_com(serializable_rep):
    save_com = serializable_rep

    return save_com



def _check_and_convert_save_potential_slices(ctor_params):
    kwargs = {"obj": ctor_params["save_potential_slices"],
              "obj_name": "save_potential_slices"}
    save_potential_slices = czekitout.convert.to_bool(**kwargs)
    
    return save_potential_slices



def _pre_serialize_save_potential_slices(save_potential_slices):
    serializable_rep = save_potential_slices

    return serializable_rep



def _de_pre_serialize_save_potential_slices(serializable_rep):
    save_potential_slices = serializable_rep

    return save_potential_slices



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The base output parameters for STEM simulations.

    For a general discussion on the possible output that can be generated from
    STEM simulations, see the documentation for the module
    :mod:`prismatique.stem.output`. That discussion provides important context
    to the description of the parameter set below.

    Parameters
    ----------
    output_dirname : `str`, optional
        The relative or absolute path to the directory in which all output files
        are to be saved. If the directory doesn't exist upon saving the output 
        files, it will be created if possible.
    max_data_size : `int`, optional
        The data size limit, in bytes, of the STEM simulation output to be
        generated. If the output to be generated would require a data size
        larger than the aforementioned limit, then an exception is raised and
        the STEM simulation is not performed. Note that data size due to HDF5
        file overhead and metadata are not taken into account.
    cbed_params : :class:`prismatique.cbed.Params` | `None`, optional
        The simulation parameters related to convengent beam electron
        diffraction patterns, which includes parameters specifying what kind of
        4D-STEM output should be saved, if any at all. If ``cbed_params`` is set
        to `None` [i.e. the default value], then the aforementioned parameters
        are set to default values.
    radial_step_size_for_3d_stem : `float`, optional
        The bin width in mrads of the annular detectors used in 3D-STEM data
        collection. If set to zero, then no 3D-STEM data is saved. Note that 
        ``radial_step_size_for_3d_stem`` must be nonnegative.
    radial_range_for_2d_stem : `array_like` (`float`, shape=(``2``,)), optional
        ``radial_range_for_2d_stem[0]`` and ``radial_range_for_2d_stem[1]`` are
        the lower and upper radial integration limits respectively in
        angular/Fourier space in units of mrads used to obtain the 2D-STEM
        intensity data. If
        ``0<=radial_range_for_2d_stem[0]<radial_range_for_2d_stem[1]``, then the
        intensity 2D-STEM data is written to a file with the basename
        ``"stem_sim_intensity_output.h5"``. If
        ``0<=radial_range_for_2d_stem[0]==radial_range_for_2d_stem[1]``, then no
        2D-STEM data is saved. In all other scenarios an error is raised.
    save_com : `bool`, optional
        If ``save_com`` is set to ``True``, then the center-of-mass (COM)
        momentum averaged over atomic configurations is written to a file with
        the basename ``"stem_sim_intensity_output.h5"``. If ``save_com`` is set 
        to ``False``, then no COM data is saved.
    save_potential_slices : `bool`, optional
        If ``save_potential_slices`` is set to ``True``, then for each frozen
        phonon configuration subset, the corresponding potential slice data is
        written to a file with the basename
        ``"potential_slices_of_subset_"+str(i)+".h5"``, with ``i`` being the 
        subset index. If ``save_potential_slices`` is set to ``False``, then no
        potential slice data is saved.

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
        {"output_dirname": _check_and_convert_output_dirname,
         "max_data_size": _check_and_convert_max_data_size,
         "cbed_params": _check_and_convert_cbed_params,
         "radial_step_size_for_3d_stem": \
         _check_and_convert_radial_step_size_for_3d_stem,
         "radial_range_for_2d_stem": \
         _check_and_convert_radial_range_for_2d_stem,
         "save_com": _check_and_convert_save_com,
         "save_potential_slices": _check_and_convert_save_potential_slices}

    _pre_serialization_funcs = \
        {"output_dirname": _pre_serialize_output_dirname,
         "max_data_size": _pre_serialize_max_data_size,
         "cbed_params": _pre_serialize_cbed_params,
         "radial_step_size_for_3d_stem": \
         _pre_serialize_radial_step_size_for_3d_stem,
         "radial_range_for_2d_stem": _pre_serialize_radial_range_for_2d_stem,
         "save_com": _pre_serialize_save_com,
         "save_potential_slices": _pre_serialize_save_potential_slices}

    _de_pre_serialization_funcs = \
        {"output_dirname": _de_pre_serialize_output_dirname,
         "max_data_size": _de_pre_serialize_max_data_size,
         "cbed_params": _de_pre_serialize_cbed_params,
         "radial_step_size_for_3d_stem": \
         _de_pre_serialize_radial_step_size_for_3d_stem,
         "radial_range_for_2d_stem": _de_pre_serialize_radial_range_for_2d_stem,
         "save_com": _de_pre_serialize_save_com,
         "save_potential_slices": _de_pre_serialize_save_potential_slices}
    
    def __init__(self,
                 output_dirname="sim_output_files",
                 max_data_size=2*10**9,
                 cbed_params=None,
                 radial_step_size_for_3d_stem=1,
                 radial_range_for_2d_stem=(0, 0),
                 save_com=False,
                 save_potential_slices=False):
        ctor_params = {"output_dirname": output_dirname,
                       "max_data_size": max_data_size,
                       "cbed_params": cbed_params,
                       "radial_step_size_for_3d_stem": \
                       radial_step_size_for_3d_stem,
                       "radial_range_for_2d_stem": radial_range_for_2d_stem,
                       "save_com": save_com,
                       "save_potential_slices": save_potential_slices}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_base_params(ctor_params):
    base_params = copy.deepcopy(ctor_params["base_params"])
    if base_params is None:
        base_params = Params()
        
    kwargs = {"obj": base_params,
              "obj_name": "base_params",
              "accepted_types": (Params, type(None))}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return base_params



def _pre_serialize_base_params(base_params):
    serializable_rep = base_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_base_params(serializable_rep):
    base_params = Params.de_pre_serialize(serializable_rep)

    return base_params
    


###########################
## Define error messages ##
###########################

_check_and_convert_radial_range_for_2d_stem_err_msg_1 = \
    ("The object ``radial_range_for_2d_stem`` must either be a pair of "
     "real numbers satisfying: "
     "``0<=radial_range_for_2d_stem[0]<=radial_range_for_2d_stem[1]``.")