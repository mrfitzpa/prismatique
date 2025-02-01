r"""For specifying simulation parameters related to convergent beam electron
diffraction patterns.

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



# For postprocessing CBED intensity patterns.
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
           "blank_unprocessed_pattern_signal"]



def _check_and_convert_postprocessing_seq(ctor_params):
    check_and_convert_postprocessing_seq = \
        prismatique._signal._check_and_convert_postprocessing_seq
    postprocessing_seq = check_and_convert_postprocessing_seq(ctor_params)
    
    return postprocessing_seq



def _pre_serialize_postprocessing_seq(postprocessing_seq):
    pre_serialize_postprocessing_seq = \
        prismatique._signal._pre_serialize_postprocessing_seq
    serializable_rep = pre_serialize_postprocessing_seq(postprocessing_seq)

    return serializable_rep



def _de_pre_serialize_postprocessing_seq(serializable_rep):
    de_pre_serialize_postprocessing_seq = \
        prismatique._signal._de_pre_serialize_postprocessing_seq
    postprocessing_seq = \
        de_pre_serialize_postprocessing_seq(serializable_rep)

    return postprocessing_seq



def _check_and_convert_avg_num_electrons_per_postprocessed_dp(ctor_params):
    kwargs = \
        {"obj": ctor_params["avg_num_electrons_per_postprocessed_dp"],
         "obj_name": "avg_num_electrons_per_postprocessed_dp"}
    avg_num_electrons_per_postprocessed_dp = \
        czekitout.convert.to_positive_float(**kwargs)
    
    return avg_num_electrons_per_postprocessed_dp



def _pre_serialize_avg_num_electrons_per_postprocessed_dp(
        avg_num_electrons_per_postprocessed_dp):
    serializable_rep = avg_num_electrons_per_postprocessed_dp

    return serializable_rep



def _de_pre_serialize_avg_num_electrons_per_postprocessed_dp(serializable_rep):
    avg_num_electrons_per_postprocessed_dp = serializable_rep

    return avg_num_electrons_per_postprocessed_dp



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
    r"""The simulation parameters related to convergent beam electron
    diffraction patterns.

    In performing STEM, a probe is scanned across a plane (i.e. 2D geometry)
    where for each probe position a 2D convergent beam electron diffraction
    (CBED) pattern is collected. The coordinates :math:`\left(\varphi_x,
    \varphi_y\right)` in the diffraction plane specify the scattering angle,
    which is related to the beam electron's transverse momentum
    :math:`\left(k_x, k_y\right)` by

    .. math ::
        \left(\varphi_x, \varphi_y\right) = \lambda \left(k_x, k_y\right)
        :label: varphi_to_k

    where :math:`\lambda` is the beam electron's wavelength. This set of CBED
    patterns is what is sometimes referred to as the 4D-STEM data/output, the
    "4D" referring to the two spatial dimensions associated with each probe
    position and the two angular dimensions associated with each CBED
    pattern. Experimentally CBED intensity patterns are collected.

    As discussed in the documentation for the class
    :class:`prismatique.thermal.Params`, the intensity pattern for a given probe
    collected by the STEM detector measures the diagonal elements of the state
    operator :math:`\hat{\rho}_{t}` of a transmitted beam electron, in the
    electron's transverse momentum basis. See Eq. :eq:`I_STEM` for a
    mathematical expression of the above. We model :math:`\hat{\rho}_{t}` by
    Eq. :eq:`mixed_state_operator_for_transmitted_electron`, wherein
    :math:`\hat{\rho}_{t}` is expressed as mixed state. Specifically, it is
    expressed as an incoherent average of pure (i.e. coherent) states
    :math:`\left|\psi_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
    \boldsymbol{\delta}_{\beta}\right)\right\rangle` over a range of defocii,
    and a set of frozen phonon configurations. See the documentation for the
    class :class:`primsatique.thermal.Params` for further discussion on
    :math:`\hat{\rho}_{t}`. To approximate the average, ``prismatic`` calculates
    the
    :math:`\left|\psi_{t}\left(\delta_{f};\mathbf{u}_{1},\ldots,\mathbf{u}_{N};
    \boldsymbol{\delta}_{\beta}\right)\right\rangle` for a discrete set of
    defocii and randomly sampled frozen phonon configurations. 

    Prior to any postprocessing, the pixel size of the CBED patterns and the
    dimensions of the CBED patterns in units of said pixels are given by
    Eq. :eq:`Delta_tilde_ks_in_prism_params` and
    :eq:`stem_detector_dims_in_pixels` respectively. See the documentation for
    the subpackage :mod:`prismatique.stem` for relevant context to the above
    equations.

    Parameters
    ----------
    postprocessing_seq : `array_like` (:class:`empix.OptionalCroppingParams` | :class:`empix.OptionalDownsamplingParams` | :class:`empix.OptionalResamplingParams`, ndim=1), optional
        Each item in ``postprocessing_seq`` specifies a postprocessing step to
        be applied to each CBED intensity pattern . Each item must be an
        instance of one of the following classes:
        :class:`empix.OptionalCroppingParams`;
        :class:`empix.OptionalDownsamplingParams`;
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
        unprocessed CBED intensity pattern set generated by the simulation is
        used as the first parameter to the implied postprocessing function,
        after being converted to a ``hyperspy`` signal. The convention used in
        prismatique is that, when converted to a ``hyperspy`` signal, the CBED
        pattern is visualized with the :math:`k_x`-axis being the horizontal
        axis, increasing from left to right, and the :math:`k_y`-axis being the
        vertical axis, increasing from bottom to top, both expressed in units of
        :math:`1/Å`.

        Blank [i.e. zeroed] unprocessed CBED patterns can be generated as
        ``hyperspy`` signals using the function
        :func:`prismatique.cbed.blank_unprocessed_pattern_signal`. This function
        may help users determine what postprocessing sequence they require to
        obtain postprocessed CBED intensity patterns with the desired pixel
        sizes, number of pixels, etc.

        Note that the parameter ``postprocessing_seq[idx].core_attrs["title"]``
        is effectively ignored for all integers ``idx`` satisfying
        ``0<=idx<len(postprocessing_seq)``.  Moreover, if
        ``postprocessing_seq[idx]`` is an instance of the class
        :class:`empix.OptionalDownsamplingParams`, then said object must satisfy
        ``postprocessing_seq[idx].core_attrs["downsample_mode"]=="mean``,
        otherwise an exception will be raised.
    avg_num_electrons_per_postprocessed_dp : `float`, optional
        The average number of electrons per postprocessed CBED intensity 
        pattern.
    apply_shot_noise : `bool`, optional
        If ``apply_shot_noise`` is set to ``True`` and ``save`` is set to
        ``"intensity"``, then simulated shot noise is applied to each CBED
        intensity pattern as a final postprocessing step, i.e. after all the
        postprocessing steps, specified by ``postprocessing_seq``, have been
        applied. Otherwise, no simulated shot noise is applied.

        Shot noise is simulated as follows: for each pixel in each CBED
        intensity pattern, the numerical value stored therein is used as the
        variance of a Poisson distribution, from which to sample a new value of
        said pixel.
    save_wavefunctions : `bool`, optional
        If ``save_wavefunctions`` is set to ``True``, then the unprocessed
        :math:`\left|\psi_{t}\left(\delta_{f};\mathbf{u}_{1},
        \ldots,\mathbf{u}_{N}; \boldsymbol{\delta}_{\beta}\right)\right\rangle`,
        represented in the electron's transverse momentum basis, are saved,
        where the wavefunction data corresponding to the ``i`` th frozen phonon
        configuration subset is saved to a file with the basename
        ``"stem_sim_wavefunction_output_of_subset_"+str(i)+".h5"``. Otherwise,
        no wavefunction data is saved.
    save_final_intensity : `bool`, optional
        If ``save_final_intensity`` is set to ``True``, then the postprocessed
        CBED intensity patterns, obtained by performing the incoherent average
        described further above, are saved to a file with basename
        ``"stem_sim_intensity_output.h5"``. Otherwise, it is not saved.

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
         "avg_num_electrons_per_postprocessed_dp": \
         _check_and_convert_avg_num_electrons_per_postprocessed_dp,
         "apply_shot_noise": _check_and_convert_apply_shot_noise,
         "save_wavefunctions": _check_and_convert_save_wavefunctions,
         "save_final_intensity": _check_and_convert_save_final_intensity}

    _pre_serialization_funcs = \
        {"postprocessing_seq": _pre_serialize_postprocessing_seq,
         "avg_num_electrons_per_postprocessed_dp": \
         _pre_serialize_avg_num_electrons_per_postprocessed_dp,
         "apply_shot_noise": _pre_serialize_apply_shot_noise,
         "save_wavefunctions": _pre_serialize_save_wavefunctions,
         "save_final_intensity": _pre_serialize_save_final_intensity}

    _de_pre_serialization_funcs = \
        {"postprocessing_seq": _de_pre_serialize_postprocessing_seq,
         "avg_num_electrons_per_postprocessed_dp": \
         _de_pre_serialize_avg_num_electrons_per_postprocessed_dp,
         "apply_shot_noise": _de_pre_serialize_apply_shot_noise,
         "save_wavefunctions": _de_pre_serialize_save_wavefunctions,
         "save_final_intensity": _de_pre_serialize_save_final_intensity}
    
    def __init__(self,
                 postprocessing_seq=tuple(),
                 avg_num_electrons_per_postprocessed_dp=1,
                 apply_shot_noise=False,
                 save_wavefunctions=False,
                 save_final_intensity=False):
        ctor_params = {"postprocessing_seq": postprocessing_seq,
                       "avg_num_electrons_per_postprocessed_dp": \
                       avg_num_electrons_per_postprocessed_dp,
                       "apply_shot_noise": apply_shot_noise,
                       "save_wavefunctions": save_wavefunctions,
                       "save_final_intensity": save_final_intensity}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_cbed_params(ctor_params):
    cbed_params = copy.deepcopy(ctor_params["cbed_params"])
    if cbed_params is None:
        cbed_params = Params()
    
    kwargs = {"obj": cbed_params,
              "obj_name": "cbed_params",
              "accepted_types": (Params, type(None))}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return cbed_params



def _pre_serialize_cbed_params(cbed_params):
    serializable_rep = cbed_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_cbed_params(serializable_rep):
    cbed_params = Params.de_pre_serialize(serializable_rep)

    return cbed_params



def blank_unprocessed_pattern_signal(sample_specification):
    r"""Generate a blank unprocessed CBED pattern as a ``hyperspy`` signal.

    This Python function may help users determine what postprocessing sequence
    they require to obtain postprocessed CBED intensity patterns with the
    desired pixel sizes, number of pixels, etc. For a discussion on
    postprocessing CBED patterns, see the documentation for the class
    :class:`prismatique.cbed.Params`.

    Parameters
    ----------
    sample_specification : :class:`prismatique.sample.ModelParams` | :class:`prismatique.sample.PotentialSliceSubsetIDs` | :class:`prismatique.sample.SMatrixSubsetIDs` | :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`
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
        :class:`prismatique.sample.PotentialSliceSubsetIDs`, the class
        :class:`prismatique.sample.SMatrixSubsetIDs`, or the class
        :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`, then
        ``sample_specification`` specifies a set of files, where each file
        stores either the pre-calculated potential slices or :math:`S`-matrices
        for a frozen phonon configuration subset. See the documentation for the
        aforementioned classes for further discussions on specifying
        pre-calculated objects. See the documentation for the subpackage
        :mod:`prismatique.stem` for a discussion on :math:`S`-matrices. 

    Returns
    -------
    blank_unprocessed_pattern_signal : :class:`hyperspy._signals.signal2d.Signal2D`
        The blank unprocessed pattern, represented as a ``hyperspy`` signal.
        The convention used in prismatique is that, when converted to a
        ``hyperspy`` signal, the CBED pattern is visualized with the
        :math:`k_x`-axis being the horizontal axis, increasing from left to
        right, and the :math:`k_y`-axis being the vertical axis, increasing from
        bottom to top, both expressed in units of :math:`1/Å`.

    """
    kwargs = \
        {"sample_specification": sample_specification,
         "navigation_dims": tuple(),
         "signal_is_cbed_pattern_set": True,
         "signal_dtype": "float"}
    blank_unprocessed_pattern_signal = \
        prismatique._signal._blank_unprocessed_2d_signal(**kwargs)

    return blank_unprocessed_pattern_signal



###########################
## Define error messages ##
###########################