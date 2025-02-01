r"""For specifying output parameters that are applicable only to the PRISM
implementation of STEM simulations.

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



def _check_and_convert_enable_S_matrix_refocus(ctor_params):
    kwargs = {"obj": ctor_params["enable_S_matrix_refocus"],
              "obj_name": "enable_S_matrix_refocus"}
    enable_S_matrix_refocus = czekitout.convert.to_bool(**kwargs)
    
    return enable_S_matrix_refocus



def _pre_serialize_enable_S_matrix_refocus(enable_S_matrix_refocus):
    serializable_rep = enable_S_matrix_refocus

    return serializable_rep



def _de_pre_serialize_enable_S_matrix_refocus(serializable_rep):
    enable_S_matrix_refocus = serializable_rep

    return enable_S_matrix_refocus



def _check_and_convert_save_S_matrices(ctor_params):
    kwargs = {"obj": ctor_params["save_S_matrices"],
              "obj_name": "save_S_matrices"}
    save_S_matrices = czekitout.convert.to_bool(**kwargs)
    
    return save_S_matrices



def _pre_serialize_save_S_matrices(save_S_matrices):
    serializable_rep = save_S_matrices

    return serializable_rep



def _de_pre_serialize_save_S_matrices(serializable_rep):
    save_S_matrices = serializable_rep

    return save_S_matrices



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The output parameters that are applicable only to the PRISM
    implementation of STEM simulations.

    For a general discussion on the possible output that can be generated from
    STEM simulations, see the documentation for the class
    :class:`prismatique.stem.output.Params`. For a general discussion on the
    PRISM algorithm and :math:`S`-matrices, see the documentation for the module
    :mod:`prismatique.stem`. These discussions provide important context to the
    description of the parameter set below.

    Parameters
    ----------
    enable_S_matrix_refocus : `bool`, optional
        If ``enable_S_matrix_refocus`` is set to ``True``, then :math:`S`-matrix
        refocusing is enabled. If ``enable_S_matrix_refocus`` is set to
        ``False``, then no :math:`S`-matrix refocusing is performed. 
    save_S_matrices : `bool`, optional
        If ``save_S_matrices`` is set to ``True``, then for each frozen phonon
        configuration subset, the corresponding :math:`S`-matrix data is written
        to the file with the basename ``"S_matrices_of_subset_"+str(i)+".h5"``,
        with ``i`` being the subset index. If ``save_S_matrices`` is set to
        ``False``, then no :math:`S`-matrix data is saved.

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
        {"enable_S_matrix_refocus": _check_and_convert_enable_S_matrix_refocus,
         "save_S_matrices": _check_and_convert_save_S_matrices}

    _pre_serialization_funcs = \
        {"enable_S_matrix_refocus": _pre_serialize_enable_S_matrix_refocus,
         "save_S_matrices": _pre_serialize_save_S_matrices}

    _de_pre_serialization_funcs = \
        {"enable_S_matrix_refocus": _de_pre_serialize_enable_S_matrix_refocus,
         "save_S_matrices": _de_pre_serialize_save_S_matrices}
    
    def __init__(self,
                 enable_S_matrix_refocus=False,
                 save_S_matrices=False):
        ctor_params = {"enable_S_matrix_refocus": enable_S_matrix_refocus,
                       "save_S_matrices": save_S_matrices}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_prism_output_params(ctor_params):
    prism_output_params = \
        copy.deepcopy(ctor_params["prism_output_params"])
    if prism_output_params is None:
        prism_output_params = Params()
        
    kwargs = {"obj": prism_output_params,
              "obj_name": "prism_output_params",
              "accepted_types": (Params, type(None))}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return prism_output_params



def _pre_serialize_prism_output_params(prism_output_params):
    serializable_rep = prism_output_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_prism_output_params(serializable_rep):
    prism_output_params = Params.de_pre_serialize(serializable_rep)

    return prism_output_params
    


###########################
## Define error messages ##
###########################