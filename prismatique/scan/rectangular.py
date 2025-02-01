r"""For specifying simulation parameters related to rectangular grid-like probe
scan patterns.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For performing deep copies of objects.
import copy



# For general array handling.
import numpy as np

# For validating and converting objects.
import czekitout.check
import czekitout.convert

# For defining classes that support enforced validation, updatability,
# pre-serialization, and de-serialization.
import fancytypes

# For checking whether two floating-point numbers are approximately equal.
import emconstants



# For validating, pre-serializing, and de-pre-serializing integers used as seeds
# to random number generators.
import prismatique.thermal

# For determining the dimensions of sample supercells.
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
__all__ = ["Params"]



def _check_and_convert_step_size(ctor_params):
    kwargs = {"obj": ctor_params["step_size"],
              "obj_name": "step_size"}
    step_size = czekitout.convert.to_pair_of_positive_floats(**kwargs)
    
    return step_size



def _pre_serialize_step_size(step_size):
    serializable_rep = step_size

    return serializable_rep



def _de_pre_serialize_step_size(serializable_rep):
    step_size = serializable_rep

    return step_size



def _check_and_convert_window(ctor_params):
    try:
        kwargs = {"obj": ctor_params["window"],
                  "obj_name": "window"}
        window = czekitout.convert.to_tuple_of_floats(**kwargs)

        if len(window) != 4:
            raise
        
        if ((not (window[0] < window[1])) or (not (window[2] < window[3]))):
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



def _check_and_convert_jitter(ctor_params):
    kwargs = {"obj": ctor_params["jitter"],
              "obj_name": "jitter"}
    jitter = czekitout.convert.to_nonnegative_float(**kwargs)
    
    return jitter



def _pre_serialize_jitter(jitter):
    serializable_rep = jitter

    return serializable_rep



def _de_pre_serialize_jitter(serializable_rep):
    jitter = serializable_rep

    return jitter



def _check_and_convert_rng_seed(ctor_params):
    rng_seed = prismatique.thermal._check_and_convert_rng_seed(ctor_params)
    
    return rng_seed



def _pre_serialize_rng_seed(rng_seed):
    serializable_rep = prismatique.thermal._pre_serialize_rng_seed(rng_seed)

    return serializable_rep



def _de_pre_serialize_rng_seed(serializable_rep):
    rng_seed = prismatique.thermal._de_pre_serialize_rng_seed(serializable_rep)

    return rng_seed



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related to rectangular grid-like probe scan
    patterns.

    A common scanning pattern in STEM is a rectangular grid. Sometimes it is
    useful to introduce small random deviations to such a scanning pattern. This
    is done in e.g. ptychography.

    The class :class:`prismatique.scan.rectangular.Params` specifies a scanning
    pattern that is obtained by generating an underlying rectangular grid to
    which to apply a random positional deviation to each point of the original
    grid, thus yielding an irregular scanning pattern that is rectangular
    grid-like, assuming the deviations are small. The radial components of the
    random positional deviations are sampled from a uniform distribution on an
    interval determined by the parameters defined below, and the polar angular
    components of the random position deviations are sampled from a uniform
    distribution on the interval :math:`\left[0, 2\pi\right)`.

    Both the
    :math:`x`- and :math:`y`-components of the random positional deviations are
    sampled from a normal distribution with zero-mean and a standard deviation
    which we refer to below as the "jitter" [the :math:`xy`-plane is assumed to
    coincide with the object plane].

    Parameters
    ----------
    step_size : `array_like` (`float`, shape=(``2``,)), optional
        ``step_size[0]`` specifies the probe scanning step size in the
        :math:`x`-direction in units of angstroms of the underlying rectangular
        grid; ``step_size[1]`` specifies the probe scanning step size in the
        :math:`y`-direction in units of angstroms of the underlying rectangular
        grid.
    window : `array_like` (`float`, shape=(``4``,)), optional
        ``window`` specifies a set of fractional coordinates that define the
        area of the underlying rectangular grid. The fractional coordinates are
        defined with respect to the spatial dimensions of the sample's supercell
        [see the documentation for the class
        :class:`prismatique.discretization.Params` for a discussion on sample
        supercells]. E.g. a fractional :math:`x`-coordinate of 0.5 corresponds
        to the midpoint of some supercell, along the :math:`x`-axis; a
        fractional :math:`x`-coordinate of 1.5 corresponds to the midpoint of
        the supercell to the right of the previous supercell, along the
        :math:`x`-axis; and a fractional :math:`x`-coordinate of -0.5
        corresponds to the midpoint of the supercell to the left of the first
        supercell mentioned, along the :math:`x`-axis.

        ``window[0]`` specifies the minimum fractional :math:`x`-coordinate of
        the scanning window of the underlying rectangular grid; ``window[1]``
        specifies the maximum fractional :math:`x`-coordinate of the scanning
        window of the underlying rectangular grid; ``window[2]`` specifies the
        minimum fractional :math:`y`-coordinate of the scanning window of the
        underlying rectangular grid; ``window[3]`` specifies the maximum
        fractional :math:`y`-coordinate of the scanning window of the underlying
        rectangular grid. Note that the tiling starts at the fractional
        coordinates ``(window[0], window[2])``.
    jitter : `int`, optional
        ``jitter``, along with the parameter ``step_size`` above, determine the
        interval from which the radial components of the random positional 
        deviations are sampled uniformly. Let :math:`\mathcal{J}`, 
        :math:`\Delta x`, and :math:`\Delta y` denote ``jitter``, 
        ``step_size[0]``, and ``step_size[1]`` respectively, then the interval 
        from which the radial components of the random positional deviations are
        sampled uniformly is 
        :math:`\left[0, \mathcal{J}\max\left(\Delta x, \Delta y\right)\right)`.
    rng_seed : `int` | `None`, optional
        A seed to pass to the random number generator used to sample the frozen
        phonon configurations. If set, ``rng_seed`` must be a non-negative 
        integer or `None`.

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
        {"step_size": _check_and_convert_step_size,
         "window": _check_and_convert_window,
         "jitter": _check_and_convert_jitter,
         "rng_seed": _check_and_convert_rng_seed}

    _pre_serialization_funcs = \
        {"step_size": _pre_serialize_step_size,
         "window": _pre_serialize_window,
         "jitter": _pre_serialize_jitter,
         "rng_seed": _pre_serialize_rng_seed}

    _de_pre_serialization_funcs = \
        {"step_size": _de_pre_serialize_step_size,
         "window": _de_pre_serialize_window,
         "jitter": _de_pre_serialize_jitter,
         "rng_seed": _de_pre_serialize_rng_seed}
    
    def __init__(self,
                 step_size=(0.25, 0.25),
                 window=(0, 1, 0, 1),
                 jitter=0,
                 rng_seed=None):
        ctor_params = {"step_size": step_size,
                       "window": window,
                       "jitter": jitter,
                       "rng_seed": rng_seed}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_rectangular_scan_params(ctor_params):
    rectangular_scan_params = \
        copy.deepcopy(ctor_params["rectangular_scan_params"])
    if rectangular_scan_params is None:
        rectangular_scan_params = Params()
    
    kwargs = {"obj": rectangular_scan_params,
              "obj_name": "rectangular_scan_params",
              "accepted_types": (Params, type(None))}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return rectangular_scan_params



def _pre_serialize_rectangular_scan_params(rectangular_scan_params):
    serializable_rep = rectangular_scan_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_rectangular_scan_params(serializable_rep):
    rectangular_scan_params = Params.de_pre_serialize(serializable_rep)

    return rectangular_scan_params



def _generate_probe_positions(rectangular_scan_params, sample_specification):
    step_size = rectangular_scan_params.core_attrs["step_size"]
    window = rectangular_scan_params.core_attrs["window"]
    jitter = rectangular_scan_params.core_attrs["jitter"]
    rng_seed = rectangular_scan_params.core_attrs["rng_seed"]

    Delta_X, Delta_Y, _ = \
        prismatique.sample._supercell_dims(sample_specification)

    min_x_probe_coord = Delta_X * window[0]
    max_x_probe_coord = Delta_X * window[1]
    min_y_probe_coord = Delta_Y * window[2]
    max_y_probe_coord = Delta_Y * window[3]

    x_probe_coord_step = step_size[0]
    y_probe_coord_step = step_size[1]
    tol = emconstants.tol

    x_coords = np.arange(min_x_probe_coord,
                         max_x_probe_coord+tol,
                         x_probe_coord_step)
    y_coords = np.arange(min_y_probe_coord,
                         max_y_probe_coord+tol,
                         y_probe_coord_step)

    np.random.seed(rng_seed)

    probe_positions = []
    for x_coord in x_coords:
        for y_coord in y_coords:
            dr = np.random.uniform(0, jitter*max(step_size))
            theta = np.random.uniform(0, 2*np.pi)
            dx = dr * np.cos(theta)
            dy = dr * np.sin(theta)
            probe_positions.append((x_coord+dx, y_coord+dy))
    probe_positions = tuple(probe_positions)

    np.random.seed()

    return probe_positions
    


###########################
## Define error messages ##
###########################

_check_and_convert_window_err_msg_1 = \
    ("The object ``window`` must be an array-like object of length 4, with "
     "each element being of type `float`, and the array as a whole satisfying: "
     "``window[0]<window[1]`` and ``window[2]<window[3]``.")