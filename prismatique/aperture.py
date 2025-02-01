r"""For specifying simulation parameters related to the objective aperture.

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



def _check_and_convert_offset(ctor_params):
    kwargs = {"obj": ctor_params["offset"],
              "obj_name": "offset"}
    offset = czekitout.convert.to_pair_of_floats(**kwargs)
    
    return offset



def _pre_serialize_offset(offset):
    serializable_rep = offset

    return serializable_rep



def _de_pre_serialize_offset(serializable_rep):
    offset = serializable_rep

    return offset



def _check_and_convert_window(ctor_params):
    try:
        kwargs = {"obj": ctor_params["window"],
                  "obj_name": "window"}
        window = \
            czekitout.convert.to_tuple_of_nonnegative_floats(**kwargs)

        if len(window) != 2:
            raise
        
        if not (window[0] <= window[1]):
            raise
    except:
        raise TypeError(_check_and_convert_window_err_msg_1)
    
    return window



def _pre_serialize_window(window):
    serializable_rep = window

    return serializable_rep



def _de_pre_serialize_window(serializable_rep):
    window = serializable_rep

    return window



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related objective aperture.

    The objective aperture is assumed to be annular, and is used in HRTEM
    simulations.

    The parameters below define an angular window within which all scattering
    angles that are not blocked by the objective aperture.

    Parameters
    ----------
    offset : array_like` (`float`, shape=(``2``,)), optional
        ``offset`` specifies the offset of the angular window of the objective
        aperture discussed above: ``offset[0]`` specifies the
        :math:`x`-coordinate of the offset in mrads; ``offset[1]`` specifies the
        :math:`y`-coordinate of the offset in mrads. 
    window : array_like` (`float`, shape=(``2``,)), optional
        If ``window`` is an array of length 2, then ``window`` specifies a
        radial window: ``window[0]`` and ``window[1]`` specify the minimum and
        maximum radial angles with respect to the offset angle in mrads.

    Attributes
    ----------
    core_attrs : `dict`, read-only
        A `dict` representation of the core attributes: each `dict` key is a
        `str` representing the name of a core attribute, and the corresponding
        `dict` value is the object to which said core attribute is set. The core
        attributes are the same as the construction parameters, except that 
        their values might have been updated since construction.

    """
    _validation_and_conversion_funcs = {"offset": _check_and_convert_offset,
                                        "window": _check_and_convert_window}

    _pre_serialization_funcs = {"offset": _pre_serialize_offset,
                                "window": _pre_serialize_window}

    _de_pre_serialization_funcs = {"offset": _de_pre_serialize_offset,
                                   "window": _de_pre_serialize_window}
    
    def __init__(self, offset=(0, 0), window=(0, float("inf"))):
        ctor_params = {"offset": offset, "window": window}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_objective_aperture_params(ctor_params):
    objective_aperture_params = \
        copy.deepcopy(ctor_params["objective_aperture_params"])
    if objective_aperture_params is None:
        objective_aperture_params = Params()
    
    kwargs = {"obj": objective_aperture_params,
              "obj_name": "objective_aperture_params",
              "accepted_types": (Params, type(None))}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return objective_aperture_params



def _pre_serialize_objective_aperture_params(objective_aperture_params):
    serializable_rep = objective_aperture_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_objective_aperture_params(serializable_rep):
    objective_aperture_params = Params.de_pre_serialize(serializable_rep)

    return objective_aperture_params
    


###########################
## Define error messages ##
###########################

_check_and_convert_window_err_msg_1 = \
    ("The object ``window`` must be a pair of real numbers satisfying: "
     "``0<=window[0]<=window[1]``.")