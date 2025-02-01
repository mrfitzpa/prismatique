r"""For specifying simulation parameters related to HRTEM image wavefunctions 
and intensities.

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



# For postprocessing HRTEM intensity images.
import prismatique._signal



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
           "blank_unprocessed_image_signal"]



def _check_and_convert_postprocessing_seq(ctor_params):
    check_and_convert_postprocessing_seq = \
        prismatique._signal._check_and_convert_postprocessing_seq
    postprocessing_seq = \
        check_and_convert_postprocessing_seq(ctor_params)
    
    return postprocessing_seq



def _pre_serialize_postprocessing_seq(postprocessing_seq):
    pre_serialize_postprocessing_seq = \
        prismatique._signal._pre_serialize_postprocessing_seq
    serializable_rep = \
        pre_serialize_postprocessing_seq(postprocessing_seq)

    return serializable_rep



def _de_pre_serialize_postprocessing_seq(serializable_rep):
    de_pre_serialize_postprocessing_seq = \
        prismatique._signal._de_pre_serialize_postprocessing_seq
    postprocessing_seq = \
        de_pre_serialize_postprocessing_seq(serializable_rep)

    return postprocessing_seq



def _check_and_convert_avg_num_electrons_per_postprocessed_image(ctor_params):
    kwargs = \
        {"obj": ctor_params["avg_num_electrons_per_postprocessed_image"],
         "obj_name": "avg_num_electrons_per_postprocessed_image"}
    avg_num_electrons_per_postprocessed_image = \
        czekitout.convert.to_positive_float(**kwargs)
    
    return avg_num_electrons_per_postprocessed_image



def _pre_serialize_avg_num_electrons_per_postprocessed_image(
        avg_num_electrons_per_postprocessed_image):
    serializable_rep = avg_num_electrons_per_postprocessed_image

    return serializable_rep



def _de_pre_serialize_avg_num_electrons_per_postprocessed_image(
        serializable_rep):
    avg_num_electrons_per_postprocessed_image = serializable_rep

    return avg_num_electrons_per_postprocessed_image



def _check_and_convert_apply_shot_noise(ctor_params):
    kwargs = {"obj": ctor_params["apply_shot_noise"],
              "obj_name": "apply_shot_noise"}
    apply_shot_noise = czekitout.convert.to_bool(**kwargs)
    
    return apply_shot_noise



def _pre_serialize_apply_shot_noise(apply_shot_noise):
    serializable_rep = apply_shot_noise

    return serializable_rep



def _de_pre_serialize_apply_shot_noise(serializable_rep):
    apply_shot_noise = serializable_rep

    return apply_shot_noise



def _check_and_convert_save_wavefunctions(ctor_params):
    kwargs = {"obj": ctor_params["save_wavefunctions"],
              "obj_name": "save_wavefunctions"}
    save_wavefunctions = czekitout.convert.to_bool(**kwargs)
    
    return save_wavefunctions



def _pre_serialize_save_wavefunctions(save_wavefunctions):
    serializable_rep = save_wavefunctions

    return serializable_rep



def _de_pre_serialize_save_wavefunctions(serializable_rep):
    save_wavefunctions = serializable_rep

    return save_wavefunctions



def _check_and_convert_save_final_intensity(ctor_params):
    kwargs = {"obj": ctor_params["save_final_intensity"],
              "obj_name": "save_final_intensity"}
    save_final_intensity = czekitout.convert.to_bool(**kwargs)
    
    return save_final_intensity



def _pre_serialize_save_final_intensity(
        save_final_intensity):
    serializable_rep = save_final_intensity

    return serializable_rep



def _de_pre_serialize_save_final_intensity(serializable_rep):
    save_final_intensity = serializable_rep

    return save_final_intensity



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related to HRTEM image wavefunctions and
    intensities.

    As discussed in the documentation for the class
    :class:`prismatique.thermal.Params`, the intensity image for a given probe
    collected by the TEM detector measures the diagonal elements of the state
    operator :math:`\hat{\rho}_{t}` of a transmitted beam electron, in the
    electron's transverse position basis. See Eq. :eq:`I_HRTEM` for a
    mathematical expression of the above. We model :math:`\hat{\rho}_{t}` by
    Eq. :eq:`mixed_state_operator_for_transmitted_electron`, wherein
    :math:`\hat{\rho}_{t}` is expressed as mixed state. Specifically, it is
    expressed as an incoherent average of pure (i.e. coherent) states
    :math:`\left|\psi_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
    \boldsymbol{\delta}_{\beta}\right)\right\rangle` over a range of defocii,
    beam tils, and a set of frozen phonon configurations. See the documentation
    for the class :class:`prismatique.thermal.Params` for further discussion on
    :math:`\hat{\rho}_{t}`. To approximate the average, ``prismatic`` calculates
    the
    :math:`\left|\psi_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
    \boldsymbol{\delta}_{\beta}\right)\right\rangle` for a discrete set of
    defocii, beam tilts, and randomly sampled frozen phonon configurations.

    Prior to any postprocessing, the pixel size of any HRTEM image wavefunction
    or image intensity is given by:
    
    .. math ::
        \Delta \tilde{x}=2 \Delta x,\quad\Delta \tilde{y}=2 \Delta y,
        :label: HRTEM_image_pixel_sizes

    where :math:`\Delta x` and :math:`\Delta y` are the potential slice or
    sample supercell pixel sizes along the :math:`x`- and :math:`y`-directions
    respectively [see the documentation for the class
    :class:`prismatique.discretization.Params` for a discussion on potential
    slices and sample supercells]. The :math:`x`- and :math:`y`-dimensions of
    any HRTEM image wavefunction or image intensity in units of pixels is given
    by respectively:

    .. math ::
        n_x=\frac{N_x}{2},\quad n_y=\frac{N_y}{2},
        :label: HRTEM_image_dims_in_pixels

    where :math:`N_x` and :math:`N_y` are the :math:`x`- and
    :math:`y`-dimensions of the sample's supercell in units of sample supercell
    pixels respectively. The factors of 2 in Eqs. :eq:`HRTEM_image_pixel_sizes`
    and :eq:`HRTEM_image_dims_in_pixels` are the result of an anti-aliasing
    operation performed in ``prismatic``.

    Parameters
    ----------
    postprocessing_seq : `array_like` (:class:`empix.OptionalCroppingParams` | :class:`empix.OptionalDownsamplingParams` | :class:`empix.OptionalResamplingParams`, ndim=1), optional
        Each item in ``postprocessing_seq`` specifies a postprocessing step to
        be applied to the HRTEM intensity image. Each item must be an instance
        of one of the following classes: :class:`empix.OptionalCroppingParams`;
        :class:`empix.OptionalDownsamplingParams`; or
        :class:`empix.OptionalResamplingParams`. If for example the
        :math:`i^{\text{th}}` item is an instance of the class
        :class:`empix.OptionalCroppingParams`, then said item specifies that at
        the :math:`i^{\text{th}}` postprocessing step, the output from the
        previous step is converted to a ``hyperspy`` signal, and passed as the
        first parameter to the function :func:`empix.crop`, with the item being
        passed as the second parameter, i.e. the optional parameters to the
        function. If the item is an instance of the class
        :class:`empix.OptionalDownsamplingParams`, then the function used is
        :func:`empix.downsample`. If the item is an instance of the class
        :class:`empix.OptionalResamplingParams`, then the function used is
        :func:`empix.resample`. Of course, for ``postprocessing_seq[0]``, the
        unprocessed HRTEM intensity image generated by the simulation is used as
        the first parameter to the implied postprocessing function, after being
        converted to a ``hyperspy`` signal. The convention used in prismatique
        is that, when converted to a ``hyperspy`` signal, the HRTEM intensity
        image or image wavefunction is visualized with the :math:`x`-axis being
        the horizontal axis, increasing from left to right, and the
        :math:`y`-axis being the vertical axis, increasing from bottom to top,
        both expressed in units of :math:`Å`.

        Blank [i.e. zeroed] unprocessed HRTEM intensity images can be generated
        as ``hyperspy`` signals using the function
        :func:`prismatique.hrtem.image.blank_unprocessed_image_signal`. This
        function may help users determine what postprocessing sequence they
        require to obtain postprocessed HRTEM intensity images with the desired
        pixel sizes, number of pixels, etc.

        Note that the parameter ``postprocessing_seq[idx].core_attrs["title"]``
        is effectively ignored for all integers ``idx`` satisfying
        ``0<=idx<len(postprocessing_seq)``.  Moreover, if
        ``postprocessing_seq[idx]`` is an instance of the class
        :class:`empix.OptionalDownsamplingParams`, then said object must satisfy
        ``postprocessing_seq[idx].core_attrs["downsample_mode"]=="mean``,
        otherwise an exception will be raised.
    avg_num_electrons_per_postprocessed_image : `float`, optional
        The average number of electrons per postprocessed HRTEM intensity image.
    apply_shot_noise : `bool`, optional
        If ``apply_shot_noise`` is set to ``True`` and ``save`` is set to
        ``"intensity"``, then simulated shot noise is applied to each HRTEM
        intensity image as a final postprocessing step, i.e. after all the
        postprocessing steps, specified by ``postprocessing_seq``, have been
        applied. Otherwise, no simulated shot noise is applied.

        Shot noise is simulated as follows: for each pixel in each HRTEM
        intensity image, the numerical value stored therein is used as the
        variance of a Poisson distribution, from which to sample a new value of
        said pixel.
    save_wavefunctions : `bool`, optional
        If ``save_wavefunctions`` is set to ``True``, then the unprocessed
        :math:`\left|\psi_{t}\left(\delta_{f};\mathbf{u}_{1},
        \ldots,\mathbf{u}_{N}; \boldsymbol{\delta}_{\beta}\right)\right\rangle`,
        represented in the electron's transverse position basis, are saved,
        where the wavefunction data corresponding to the ``i`` th frozen phonon
        configuration subset is saved to a file with the basename
        ``"hrtem_sim_wavefunction_output_of_subset_"+str(i)+".h5"``. Otherwise,
        no wavefunction data is saved.
    save_final_intensity : `bool`, optional
        If ``save_final_intensity`` is set to ``True``, then the postprocessed
        HRTEM intensity image, obtained by performing the incoherent average
        described further above, is saved to a file with basename
        ``"hrtem_sim_intensity_output.h5"``. Otherwise, it is not saved.

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
        {"postprocessing_seq": _check_and_convert_postprocessing_seq,
         "avg_num_electrons_per_postprocessed_image": \
         _check_and_convert_avg_num_electrons_per_postprocessed_image,
         "apply_shot_noise": _check_and_convert_apply_shot_noise,
         "save_wavefunctions": _check_and_convert_save_wavefunctions,
         "save_final_intensity": _check_and_convert_save_final_intensity}

    _pre_serialization_funcs = \
        {"postprocessing_seq": _pre_serialize_postprocessing_seq,
         "avg_num_electrons_per_postprocessed_image": \
         _pre_serialize_avg_num_electrons_per_postprocessed_image,
         "apply_shot_noise": _pre_serialize_apply_shot_noise,
         "save_wavefunctions": _pre_serialize_save_wavefunctions,
         "save_final_intensity": _pre_serialize_save_final_intensity}

    _de_pre_serialization_funcs = \
        {"postprocessing_seq": _de_pre_serialize_postprocessing_seq,
         "avg_num_electrons_per_postprocessed_image": \
         _de_pre_serialize_avg_num_electrons_per_postprocessed_image,
         "apply_shot_noise": _de_pre_serialize_apply_shot_noise,
         "save_wavefunctions": _de_pre_serialize_save_wavefunctions,
         "save_final_intensity": _de_pre_serialize_save_final_intensity}
    
    def __init__(self,
                 postprocessing_seq=tuple(),
                 avg_num_electrons_per_postprocessed_image=1,
                 apply_shot_noise=False,
                 save_wavefunctions=False,
                 save_final_intensity=False):
        ctor_params = {"postprocessing_seq": postprocessing_seq,
                       "avg_num_electrons_per_postprocessed_image": \
                       avg_num_electrons_per_postprocessed_image,
                       "apply_shot_noise": apply_shot_noise,
                       "save_wavefunctions": save_wavefunctions,
                       "save_final_intensity": save_final_intensity}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_image_params(ctor_params):
    image_params = copy.deepcopy(ctor_params["image_params"])
    if image_params is None:
        image_params = Params()
    
    kwargs = {"obj": image_params,
              "obj_name": "image_params",
              "accepted_types": (Params, type(None))}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return image_params



def _pre_serialize_image_params(image_params):
    serializable_rep = image_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_image_params(serializable_rep):
    image_params = Params.de_pre_serialize(serializable_rep)

    return image_params



def blank_unprocessed_image_signal(sample_specification):
    r"""Generate a blank unprocessed HRTEM intensity image as a ``hyperspy`` 
    signal.

    This Python function may help users determine what postprocessing sequence
    they require to obtain postprocessed HRTEM intensity images with the desired
    pixel sizes, number of pixels, etc. For a discussion on postprocessing HRTEM
    intensity images, see the documentation for the class
    :class:`prismatique.hrtem.image.Params`.

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

    Returns
    -------
    blank_unprocessed_image_signal : :class:`hyperspy._signals.signal2d.Signal2D`
        The blank unprocessed image, represented as a ``hyperspy`` signal.  The
        convention used in prismatique is that, when converted to a ``hyperspy``
        signal, the HRTEM intensity image is visualized with the :math:`x`-axis
        being the horizontal axis, increasing from left to right, and the
        :math:`y`-axis being the vertical axis, increasing from bottom to top,
        both expressed in units of :math:`Å`.

    """
    kwargs = \
        {"sample_specification": sample_specification,
         "navigation_dims": tuple(),
         "signal_is_cbed_pattern_set": False,
         "signal_dtype": "float"}
    blank_unprocessed_image_signal = \
        prismatique._signal._blank_unprocessed_2d_signal(**kwargs)

    return blank_unprocessed_image_signal



###########################
## Define error messages ##
###########################