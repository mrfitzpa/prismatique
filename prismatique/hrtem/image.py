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
r"""For specifying simulation parameters related to HRTEM image wavefunctions 
and intensities.

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



# For postprocessing HRTEM intensity images.
import prismatique._signal

# For recycling helper functions and/or constants.
import prismatique.cbed
import prismatique.sample



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params",
           "blank_unprocessed_image_signal"]



def _check_and_convert_postprocessing_seq(params):
    module_alias = prismatique.cbed
    func_alias = module_alias._check_and_convert_postprocessing_seq
    postprocessing_seq = func_alias(params)

    return postprocessing_seq



def _pre_serialize_postprocessing_seq(postprocessing_seq):
    obj_to_pre_serialize = postprocessing_seq
    module_alias = prismatique.cbed
    func_alias = module_alias._pre_serialize_postprocessing_seq
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_postprocessing_seq(serializable_rep):
    module_alias = prismatique.cbed
    func_alias = module_alias._de_pre_serialize_postprocessing_seq
    postprocessing_seq = func_alias(serializable_rep)

    return postprocessing_seq



def _check_and_convert_avg_num_electrons_per_postprocessed_image(params):
    obj_name = \
        "avg_num_electrons_per_postprocessed_image"
    kwargs = \
        {"obj": params[obj_name], "obj_name": obj_name}
    avg_num_electrons_per_postprocessed_image = \
        czekitout.convert.to_positive_float(**kwargs)
    
    return avg_num_electrons_per_postprocessed_image



def _pre_serialize_avg_num_electrons_per_postprocessed_image(
        avg_num_electrons_per_postprocessed_image):
    obj_to_pre_serialize = avg_num_electrons_per_postprocessed_image
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_avg_num_electrons_per_postprocessed_image(
        serializable_rep):
    avg_num_electrons_per_postprocessed_image = serializable_rep

    return avg_num_electrons_per_postprocessed_image



def _check_and_convert_apply_shot_noise(params):
    module_alias = prismatique.cbed
    func_alias = module_alias._check_and_convert_apply_shot_noise
    apply_shot_noise = func_alias(params)

    return apply_shot_noise



def _pre_serialize_apply_shot_noise(apply_shot_noise):
    obj_to_pre_serialize = apply_shot_noise
    module_alias = prismatique.cbed
    func_alias = module_alias._pre_serialize_apply_shot_noise
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_apply_shot_noise(serializable_rep):
    module_alias = prismatique.cbed
    func_alias = module_alias._de_pre_serialize_apply_shot_noise
    apply_shot_noise = func_alias(serializable_rep)

    return apply_shot_noise



def _check_and_convert_save_wavefunctions(params):
    module_alias = prismatique.cbed
    func_alias = module_alias._check_and_convert_save_wavefunctions
    save_wavefunctions = func_alias(params)

    return save_wavefunctions



def _pre_serialize_save_wavefunctions(save_wavefunctions):
    obj_to_pre_serialize = save_wavefunctions
    module_alias = prismatique.cbed
    func_alias = module_alias._pre_serialize_save_wavefunctions
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_save_wavefunctions(serializable_rep):
    module_alias = prismatique.cbed
    func_alias = module_alias._de_pre_serialize_save_wavefunctions
    save_wavefunctions = func_alias(serializable_rep)

    return save_wavefunctions



def _check_and_convert_save_final_intensity(params):
    module_alias = prismatique.cbed
    func_alias = module_alias._check_and_convert_save_final_intensity
    save_final_intensity = func_alias(params)

    return save_final_intensity



def _pre_serialize_save_final_intensity(save_final_intensity):
    obj_to_pre_serialize = save_final_intensity
    module_alias = prismatique.cbed
    func_alias = module_alias._pre_serialize_save_final_intensity
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_save_final_intensity(serializable_rep):
    module_alias = prismatique.cbed
    func_alias = module_alias._de_pre_serialize_save_final_intensity
    save_final_intensity = func_alias(serializable_rep)

    return save_final_intensity



_module_alias = \
    prismatique.cbed
_default_postprocessing_seq = \
    _module_alias._default_postprocessing_seq
_default_avg_num_electrons_per_postprocessed_image = \
    1
_default_apply_shot_noise = \
    _module_alias._default_apply_shot_noise
_default_save_wavefunctions = \
    _module_alias._default_save_wavefunctions
_default_save_final_intensity = \
    _module_alias._default_save_final_intensity
_default_skip_validation_and_conversion = \
    _module_alias._default_skip_validation_and_conversion



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
    ctor_param_names = ("postprocessing_seq",
                        "avg_num_electrons_per_postprocessed_image",
                        "apply_shot_noise",
                        "save_wavefunctions",
                        "save_final_intensity")
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
                 postprocessing_seq=\
                 _default_postprocessing_seq,
                 avg_num_electrons_per_postprocessed_image=\
                 _default_avg_num_electrons_per_postprocessed_image,
                 apply_shot_noise=\
                 _default_apply_shot_noise,
                 save_wavefunctions=\
                 _default_save_wavefunctions,
                 save_final_intensity=\
                 _default_save_final_intensity,
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



def _check_and_convert_image_params(params):
    obj_name = "image_params"
    obj = params[obj_name]

    accepted_types = (Params, type(None))
    
    if isinstance(obj, accepted_types[-1]):
        image_params = accepted_types[0]()
    else:
        kwargs = {"obj": obj,
                  "obj_name": obj_name,
                  "accepted_types": accepted_types}
        czekitout.check.if_instance_of_any_accepted_types(**kwargs)

        kwargs = obj.get_core_attrs(deep_copy=False)
        image_params = accepted_types[0](**kwargs)

    return image_params



def _pre_serialize_image_params(image_params):
    obj_to_pre_serialize = image_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_image_params(serializable_rep):
    image_params = Params.de_pre_serialize(serializable_rep)
    
    return image_params



_default_image_params = None



def _check_and_convert_sample_specification(params):
    params["accepted_types"] = (prismatique.sample.ModelParams,
                                prismatique.sample.PotentialSliceSubsetIDs)

    module_alias = prismatique.sample
    func_alias = module_alias._check_and_convert_sample_specification
    sample_specification = func_alias(params)

    del params["accepted_types"]

    return sample_specification



def _check_and_convert_skip_validation_and_conversion(params):
    module_alias = prismatique.sample
    func_alias = module_alias._check_and_convert_skip_validation_and_conversion
    skip_validation_and_conversion = func_alias(params)

    return skip_validation_and_conversion



def blank_unprocessed_image_signal(sample_specification,
                                   skip_validation_and_conversion=\
                                   _default_skip_validation_and_conversion):
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
    image_signal : :class:`hyperspy._signals.signal2d.Signal2D`
        The blank unprocessed image, represented as a ``hyperspy`` signal.  The
        convention used in prismatique is that, when converted to a ``hyperspy``
        signal, the HRTEM intensity image is visualized with the :math:`x`-axis
        being the horizontal axis, increasing from left to right, and the
        :math:`y`-axis being the vertical axis, increasing from bottom to top,
        both expressed in units of :math:`Å`.

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
    image_signal = _blank_unprocessed_image_signal(**kwargs)

    return image_signal



def _blank_unprocessed_image_signal(sample_specification):
    kwargs = {"sample_specification": sample_specification,
              "navigation_dims": tuple(),
              "signal_is_cbed_pattern_set": False,
              "signal_dtype": "float"}
    image_signal = prismatique._signal._blank_unprocessed_2d_signal(**kwargs)

    return image_signal



###########################
## Define error messages ##
###########################
