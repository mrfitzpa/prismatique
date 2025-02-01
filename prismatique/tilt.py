r"""For specifying simulation parameters and calculating quantities related to
beam tilt series in HRTEM.

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

# For calculating the electron beam wavelength and validating attributes of the
# class :class:`embeam.gun.ModelParams`.
import embeam



# For validating instances of the classes
# :class:`prismatique.sample.ModelParams`, and
# :class:`prismatique.sample.PotentialSliceSubsetIDs``; and for calculating
# quantities related to the modelling of the sample.
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
__all__ = ["Params",
           "step_size",
           "series"]



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
        window = czekitout.convert.to_tuple_of_nonnegative_floats(**kwargs)

        if (len(window) != 2) and (len(window) != 4):
            raise
        
        if ((not (0 <= window[0] <= window[1]))
            or (not (0 <= window[-2] <= window[-1]))):
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



def _check_and_convert_spread(ctor_params):
    kwargs = {"obj": ctor_params["spread"],
              "obj_name": "spread"}
    spread = czekitout.convert.to_nonnegative_float(**kwargs)
    
    return spread



def _pre_serialize_spread(spread):
    serializable_rep = spread

    return serializable_rep



def _de_pre_serialize_spread(serializable_rep):
    spread = serializable_rep

    return spread



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related to the beam tilt series in the HRTEM
    simulation.

    As discussed in the documentation for the class
    :class:`prismatique.thermal.Params`, due to finite source effects from the
    electron gun, the electron beam in HRTEM experiments will generally
    illuminate the sample with a small distribution of angles. The resulting
    decoherence effects are encoded in the expression we use for the state
    operator for a transmitted beam electron,
    Eq. :eq:`mixed_state_operator_for_transmitted_electron`, wherein an
    incoherent averaging over a set of beam tilt angles is performed.

    The incident beam's tilt angles with the :math:`x`- and :math:`y`-axes,
    :math:`\theta_x` and :math:`\theta_y`, are related to the incident beam
    electron momentum by:

    .. math ::
        \theta_x = \lambda k_x,\quad\theta_y = \lambda k_y,
        :label: k_to_theta_in_tilt_params
    
    where :math:`\lambda` is the incident beam electron wavelength, and
    :math:`k_x` and :math:`k_y` are the :math:`x`- and :math:`y`-components of
    the incident beam electron momentum.

    As discussed in the documentation for the class
    :class:`prismatique.discretization.Params`, real-space [and thus
    momentum/Fourier/angular-space] need to be discretized in order to handle
    wavefunctions and probabilities numerically. The smallest possible
    Fourier-space pixel sizes in the :math:`k_{x}`- and
    :math:`k_{y}`-directions, which we denote here as :math:`\Delta k_{x}` and
    :math:`\Delta k_{y}` respectively, are determined by the :math:`x`- and
    :math:`y`-dimensions of the sample's supercell:

    .. math ::
        \Delta k_{x}=\frac{1}{\Delta X},\quad\Delta k_{y}=\frac{1}{\Delta Y},
        :label: Delta_ks_in_tilt_params

    where :math:`\Delta X` and :math:`\Delta Y` are the :math:`x`- and
    :math:`y`-dimensions of the sample's supercell in units of length [see the
    documentation for the class :class:`prismatique.discretization.Params` for a
    discussion on supercells]. An angular space
    :math:`\boldsymbol{\Theta}_{f_{x},f_{y}}` that is useful to our current
    discussion is the set of angles:

    .. math ::
        \boldsymbol{\Theta}_{f_{x},f_{y}}=\left\{ \left.\left(l_{x}f_{x}\lambda
        \Delta k_{x},l_{y}f_{y}\lambda
        \Delta k_{y}\right)\right|l_{x},l_{y}\in\mathbb{Z}\right\},
        :label: discretized_angular_space

    where :math:`f_{x}` and :math:`f_{y}` are positive integers called the
    :math:`x` and :math:`y` interpolation factors, which are introduced in the
    documentation for the class :class:`prismatique.discretization.Params`.

    In ``prismatique``, users can select all the angles from either a
    rectangular or radial window/region of the discretized angular-space
    :math:`\boldsymbol{\Theta}_{f_{x},f_{y}}` as a set of beam tilt angles for
    which to perform HRTEM simulations. Users can generate beam tilt series to
    conveniently model spatially coherent HRTEM experiments at different beam
    tilts, or to model a single spatially incoherent HRTEM beam, for which to
    use in simulations.

    Parameters
    ----------
    offset : array_like` (`float`, shape=(``2``,)), optional
        ``offset`` specifies the offset of the window for the beam tilt series
        discussed above: ``offset[0]`` specifies the :math:`x`-coordinate of the
        offset in mrads; ``offset[1]`` specifies the :math:`y`-coordinate of the
        offset in mrads.
    window : array_like` (`float`, shape=(``2``,)) | array_like` (`float`, shape=(``4``,)), optional
        If ``window`` is an array of length 2, then ``window`` specifies a
        radial window for the beam tilt series: ``window[0]`` and ``window[1]``
        specify the minimum and maximum possible radial tilts with respect to 
        the tilt offset in mrads.

        If ``window`` is an array of length 4, then ``window`` specifies a
        rectangular window for the beam tilt series: ``window[0]`` and
        ``window[1]`` specify the minimum and maximum possible absolute values
        of the :math:`x`-tilt with respect to the tilt offset in mrads;
        ``window[2]`` and ``window[3]`` specify the minimum and maximum possible
        absolute values of the :math:`y`-tilt with respect to the tilt offset in
        mrads.
    spread : `float`, optional
        The beam tilt spread :math:`\sigma_{\beta}` in units of mrads,
        introduced in Eq. :eq:`p_sigma_beta`. Must be nonnegative.

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
                                        "window": _check_and_convert_window,
                                        "spread": _check_and_convert_spread}

    _pre_serialization_funcs = {"offset": _pre_serialize_offset,
                                "window": _pre_serialize_window,
                                "spread": _pre_serialize_spread}

    _de_pre_serialization_funcs = {"offset": _de_pre_serialize_offset,
                                   "window": _de_pre_serialize_window,
                                   "spread": _de_pre_serialize_spread}
    
    def __init__(self, offset=(0, 0), window=(0, 0), spread=0):
        ctor_params = {"offset": offset, "window": window, "spread": spread}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_tilt_params(ctor_params):
    tilt_params = copy.deepcopy(ctor_params["tilt_params"])
    if tilt_params is None:
        tilt_params = Params()
    
    kwargs = {"obj": tilt_params,
              "obj_name": "tilt_params",
              "accepted_types": (Params, type(None))}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return tilt_params



def _pre_serialize_tilt_params(tilt_params):
    serializable_rep = tilt_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_tilt_params(serializable_rep):
    tilt_params = Params.de_pre_serialize(serializable_rep)

    return tilt_params



def step_size(sample_specification, mean_beam_energy):
    r"""Determine beam tilt series step size in the HRTEM simulation from a 
    subset of the simulation parameters.

    For additional context on the tilt step size, see the documentation for the
    class :class:`prismatique.tilt.Params`.

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

          * atomic_coords_filename

          * unit_cell_tiling

          * discretization_params

            + interpolation_factors

        Otherwise, if ``sample_specification`` is an instance of the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs` then
        ``sample_specification`` specifies a set of files, where each file
        stores the pre-calculated potential slices for a frozen phonon
        configuration subset. See the documentation for the aforementioned
        class for a further discussion on specifying pre-calculated
        potential slices. 
    mean_beam_energy : `float`
        The mean electron beam energy :math:`E` in units of keV. Must be
        positive.

    Returns
    -------
    tilt_step_size : `array_like` (`float`, shape=(``2``,))
        ``tilt_step_size[0]`` and ``tilt_step_size[1]`` are the tilt step sizes 
        in the :math:`x`- and :math:`y`-directions in units of mrads 
        respectively.

    """
    sample_specification, mean_beam_energy = \
        _check_sample_specification_and_mean_beam_energy(sample_specification,
                                                         mean_beam_energy)

    tilt_step_size = _step_size(sample_specification, mean_beam_energy)

    return tilt_step_size



def _check_sample_specification_and_mean_beam_energy(sample_specification,
                                                     mean_beam_energy):
    accepted_types = (prismatique.sample.ModelParams,
                      prismatique.sample.PotentialSliceSubsetIDs)
    prismatique.sample._check_sample_specification(sample_specification,
                                                   accepted_types)
        
    temp_ctor_params = \
        {"mean_beam_energy": mean_beam_energy}
    mean_beam_energy = \
        embeam.gun._check_and_convert_mean_beam_energy(temp_ctor_params)

    return sample_specification, mean_beam_energy



def _step_size(sample_specification, mean_beam_energy):
    if "thermal_params" in sample_specification.core_attrs:
        discretization_params = \
            sample_specification.core_attrs["discretization_params"]
        interpolation_factors = \
            np.array(discretization_params.core_attrs["interpolation_factors"])
    else:
        interpolation_factors = \
            np.array(sample_specification.core_attrs["interpolation_factors"])

    smallest_possible_angular_space_pixel_size = \
        _smallest_possible_angular_space_pixel_size(sample_specification,
                                                    mean_beam_energy)
    tilt_step_size = \
        smallest_possible_angular_space_pixel_size * interpolation_factors
    tilt_step_size = \
        tuple(float(elem) for elem in tilt_step_size)

    return tilt_step_size



def _smallest_possible_angular_space_pixel_size(sample_specification,
                                                mean_beam_energy):
    sample_supercell_dims = \
        prismatique.sample._supercell_dims(sample_specification)
    electron_beam_wavelength = \
        embeam.wavelength(mean_beam_energy)
    result = electron_beam_wavelength / sample_supercell_dims[:2] * 1000
    result = tuple(float(elem) for elem in result)

    return result



def series(sample_specification, mean_beam_energy, tilt_params=None):
    r"""Determine the beam tilt series in the HRTEM simulation from a subset of 
    the simulation parameters.

    For additional context on beam tilts, see the documentation for the class
    :class:`prismatique.tilt.Params`.

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

          * atomic_coords_filename

          * unit_cell_tiling

          * discretization_params

            + sample_supercell_reduced_xy_dims_in_pixels

            + interpolation_factors

        Otherwise, if ``sample_specification`` is an instance of the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs` then
        ``sample_specification`` specifies a set of files, where each file
        stores the pre-calculated potential slices for a frozen phonon
        configuration subset. See the documentation for the aforementioned
        class for a further discussion on specifying pre-calculated
        potential slices. 
    mean_beam_energy : `float`
        The mean electron beam energy :math:`E` in units of keV. Must be
        positive.
    tilt_params : :class:`prismatique.tilt.Params` | `None`, optional
        The simulation parameters related to the beam tilt series in the HRTEM
        simulation to model a set of spatially coherent HRTEM experiments at
        different beam tilts, or to model a single spatially incoherent HRTEM
        beam. See the documentation for the class
        :class:`prismatique.tilt.Params` for a discussion on said parameters.
        If ``tilt_params`` is set to `None` [i.e. the default value], then the
        aforementioned simulation parameters are set to default values.

    Returns
    -------
    tilt_series : `array_like` (`float`, shape=(``num_tilts``, ``2``))
        If we let ``num_tilts`` be the number of beam tilts, then
        ``tilt_series[i][0]`` and ``tilt_series[i][1]`` are the beam tilt angles
        along the :math:`x`- and :math:`y`-axes respectively for the ``i`` th
        beam tilt in units of mrad, where ``0<=i<num_tilts``.

    """
    sample_specification, mean_beam_energy = \
        _check_sample_specification_and_mean_beam_energy(sample_specification,
                                                         mean_beam_energy)
    tilt_params = \
        _check_and_convert_tilt_params({"tilt_params": tilt_params})

    tilt_series = _series(sample_specification, mean_beam_energy, tilt_params)

    return tilt_series



def _series(sample_specification, mean_beam_energy, tilt_params):
    tilt_offset = tilt_params.core_attrs["offset"]
    tilt_window = tilt_params.core_attrs["window"]
    
    for_a_hrtem_calculation = True

    angular_mesh, beam_mask = \
        prismatique.sample._angular_mesh_and_beam_mask(sample_specification,
                                                       mean_beam_energy,
                                                       tilt_offset,
                                                       tilt_window,
                                                       for_a_hrtem_calculation)

    tilt_series = []
    for k_x_idx in range(beam_mask.shape[0]):
        for k_y_idx in range(beam_mask.shape[1]):
            if beam_mask[k_x_idx][k_y_idx]:
                x_tilt = angular_mesh[0][k_x_idx][k_y_idx]  # In rads.
                y_tilt = angular_mesh[1][k_x_idx][k_y_idx]  # In rads.
                tilt_series.append((x_tilt*1000, y_tilt*1000))  # In mrads.

    tilt_series = np.array(tilt_series)
    reordered_indices = np.lexsort((tilt_series[:, 1], tilt_series[:, 0]))
    tilt_series = tilt_series[reordered_indices]

    return tilt_series



###########################
## Define error messages ##
###########################

_check_and_convert_window_err_msg_1 = \
    ("The object ``window`` must be an array-like object of length 2 or 4, "
     "with each element being of type `float`, and the array as a whole "
     "satisfying: ``0<=window[0]<=window[1]`` and "
     "``0<=window[-2]<=window[-1]``.")