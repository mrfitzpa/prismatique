r"""Contains helper functions for postprocessing CBED intensity patterns and
HRTEM intensity images.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For general array handling.
import numpy as np

# For checking whether floats are approximately zero.
import emconstants

# For postprocessing CBED intensity patterns and HRTEM intensity images.
import empix

# For creating ``hyperspy`` axes and signals.
import hyperspy.axes
import hyperspy.signals

# For validating objects.
import czekitout.check



# For validating instances of the classes
# :class:`prismatique.sample.ModelParams`,
# :class:`prismatique.sample.PotentialSliceSubsetIDs`,
# :class:`prismatique.sample.SMatrixSubsetIDs`, and
# :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`; and for
# calculating quantities related to the modelling of the sample.
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
__all__ = []



def _check_and_convert_postprocessing_seq(ctor_params):
    try:
        postprocessing_seq = ctor_params["postprocessing_seq"]
        accepted_types = (empix.OptionalCroppingParams,
                          empix.OptionalDownsamplingParams,
                          empix.OptionalResamplingParams)

        for postprocessing_step in postprocessing_seq:
            kwargs = {"obj": postprocessing_step,
                      "obj_name": "postprocessing_step",
                      "accepted_types": accepted_types}
            czekitout.check.if_instance_of_any_accepted_types(**kwargs)
    except:
        raise TypeError(_check_and_convert_postprocessing_seq_err_msg_1)

    for idx, postprocessing_step in enumerate(postprocessing_seq):
        if isinstance(postprocessing_step, empix.OptionalDownsamplingParams):
            if postprocessing_step.core_attrs["downsample_mode"] != "mean":
                unformatted_err_msg = \
                    _check_and_convert_postprocessing_seq_err_msg_2
                err_msg = \
                    unformatted_err_msg.format(idx, idx)
                raise ValueError(err_msg)
    
    return postprocessing_seq



def _pre_serialize_postprocessing_seq(postprocessing_seq):
    serializable_rep = []
    
    for postprocessing_step in postprocessing_seq:
        serializable_rep.append(postprocessing_step.pre_serialize())
        
    serializable_rep = tuple(serializable_rep)

    return serializable_rep



def _de_pre_serialize_postprocessing_seq(serializable_rep):
    postprocessing_seq = []
    
    for serialized_postprocessing_step in serializable_rep:
        if "center" in serialized_postprocessing_step:
            cls = empix.OptionalCroppingParams
        elif "block_dims" in serialized_postprocessing_step:
            cls = empix.OptionalDownsamplingParams
        else:
            cls = empix.OptionalResamplingParams
            
        postprocessing_step = \
            cls.de_pre_serialize(serialized_postprocessing_step)
        postprocessing_seq.append(postprocessing_step)
        
    postprocessing_seq = tuple(postprocessing_seq)

    return postprocessing_seq



def _blank_unprocessed_2d_signal(sample_specification,
                                 navigation_dims,
                                 signal_is_cbed_pattern_set,
                                 signal_dtype):
    kwargs = {"sample_specification": sample_specification,
              "signal_is_cbed_pattern_set": signal_is_cbed_pattern_set}
    _check_sample_specification(**kwargs)
    signal_space_axes = _signal_space_axes_of_unprocessed_2d_signal(**kwargs)

    signal_shape = (list(navigation_dims)
                    + [signal_space_axes[1].size, signal_space_axes[0].size])
    zeros = np.zeros(signal_shape)

    if signal_dtype == "float":
        blank_unprocessed_2d_signal = \
            hyperspy.signals.Signal2D(data=zeros)
    else:
        blank_unprocessed_2d_signal = \
            hyperspy.signals.ComplexSignal2D(data=zeros)

    for idx, axis in enumerate(signal_space_axes):
        blank_unprocessed_2d_signal.axes_manager[-2+idx].update_from(axis)
        blank_unprocessed_2d_signal.axes_manager[-2+idx].name = axis.name

    return blank_unprocessed_2d_signal



def _check_sample_specification(sample_specification,
                                signal_is_cbed_pattern_set):
    mod_alias = prismatique.sample
    accepted_types = [mod_alias.ModelParams,
                      mod_alias.PotentialSliceSubsetIDs]
    if signal_is_cbed_pattern_set:
        accepted_types.append(mod_alias.SMatrixSubsetIDs)
        accepted_types.append(mod_alias.PotentialSliceAndSMatrixSubsetIDs)
    accepted_types = tuple(accepted_types)
    mod_alias._check_sample_specification(sample_specification, accepted_types)

    return None



def _signal_space_axes_of_unprocessed_2d_signal(sample_specification,
                                                signal_is_cbed_pattern_set):
    if signal_is_cbed_pattern_set:
        signal_space_axes = \
            _k_x_and_k_y_axes_of_unprocessed_dp_signal(sample_specification)
    else:
        signal_space_axes = \
            _r_x_and_r_y_axes_of_unprocessed_image_signal(sample_specification)

    return signal_space_axes



def _k_x_and_k_y_axes_of_unprocessed_dp_signal(sample_specification):
    k_x, k_y = _k_x_and_k_y_of_unprocessed_dp_signal(sample_specification)

    axes_labels = (r"$k_x$", r"$k_y$")
    sizes = (len(k_x), len(k_y))
    scales = (k_x[1]-k_x[0], -(k_y[-2]-k_y[-1]))
    offsets = (k_x[0], k_y[0])
    units = ("1/Å", "1/Å")

    if abs(scales[0]+scales[1]) < emconstants.tol:
        scales = (scales[0], -scales[0])

    k_x_and_k_y_axes = []
    for axis_idx in range(len(units)):
        axis = hyperspy.axes.UniformDataAxis(size=sizes[axis_idx],
                                             scale=scales[axis_idx],
                                             offset=offsets[axis_idx],
                                             units=units[axis_idx])
        axis.name = axes_labels[axis_idx]
        k_x_and_k_y_axes.append(axis)

    k_x_axis, k_y_axis = k_x_and_k_y_axes

    return k_x_axis, k_y_axis



def _k_x_and_k_y_of_unprocessed_dp_signal(sample_specification):
    mod_alias = prismatique.sample

    kwargs = \
        {"sample_specification": sample_specification}
    sample_supercell_xy_dims_in_pixels = \
        mod_alias._supercell_xy_dims_in_pixels(**kwargs)
    sample_supercell_lateral_pixel_size = \
        mod_alias._supercell_lateral_pixel_size(**kwargs)
    f_x, f_y = \
        mod_alias._interpolation_factors_from_sample_specification(**kwargs)

    # Fourier coordinates before anti-aliasing crop.
    k_x = mod_alias._FFT_1D_freqs(sample_supercell_xy_dims_in_pixels[0],
                                  sample_supercell_lateral_pixel_size[0])
    k_y = mod_alias._FFT_1D_freqs(sample_supercell_xy_dims_in_pixels[1],
                                  sample_supercell_lateral_pixel_size[1])

    # Fourier coordinates after anti-aliasing crop.
    N_x = len(k_x)
    k_x = np.concatenate((k_x[:(N_x//4)], k_x[-(N_x//4):]))
    N_y = len(k_y)
    k_y = np.concatenate((k_y[:(N_y//4)], k_y[-(N_y//4):]))

    # Fourier coordinates after anti-aliasing crop and interpolation.
    k_x = k_x[::f_x]
    k_y = k_y[::f_y]

    # Sort ``k_x`` and ``k_y`` in ascending and descending order respectively.
    k_x = np.sort(k_x)
    k_y = np.sort(k_y)[::-1]

    return k_x, k_y



def _r_x_and_r_y_axes_of_unprocessed_image_signal(sample_specification):
    r_x, r_y = _r_x_and_r_y_of_unprocessed_image_signal(sample_specification)

    axes_labels = (r"$x$", r"$y$")
    sizes = (len(r_x), len(r_y))
    scales = (r_x[1]-r_x[0], -(r_y[-2]-r_y[-1]))
    offsets = (r_x[0], r_y[0])
    units = ("Å", "Å")

    if abs(scales[0]+scales[1]) < emconstants.tol:
        scales = (scales[0], -scales[0])

    r_x_and_r_y_axes = []
    for axis_idx in range(len(units)):
        axis = hyperspy.axes.UniformDataAxis(size=sizes[axis_idx],
                                             scale=scales[axis_idx],
                                             offset=offsets[axis_idx],
                                             units=units[axis_idx])
        axis.name = axes_labels[axis_idx]
        r_x_and_r_y_axes.append(axis)

    r_x_axis, r_y_axis = r_x_and_r_y_axes

    return r_x_axis, r_y_axis



def _r_x_and_r_y_of_unprocessed_image_signal(sample_specification):
    kwargs = \
        {"sample_specification": sample_specification}
    N_x, N_y = \
        prismatique.sample._supercell_xy_dims_in_pixels(**kwargs)
    Delta_x, Delta_y = \
        prismatique.sample._supercell_lateral_pixel_size(**kwargs)

    n_x = N_x / 2
    n_y = N_y / 2

    Delta_tilde_x = 2 * Delta_x
    Delta_tilde_y = 2 * Delta_y

    r_x = Delta_tilde_x * np.arange(n_x) - (((n_x-1) * Delta_tilde_x) / 2)
    r_y = -(Delta_tilde_y * np.arange(n_y) - (((n_y-1) * Delta_tilde_y) / 2))

    return r_x, r_y



def _num_pixels_in_postprocessed_2d_signal_space(sample_specification,
                                                 signal_is_cbed_pattern_set,
                                                 postprocessing_seq):
    kwargs = {"sample_specification": sample_specification,
              "postprocessing_seq": postprocessing_seq,
              "navigation_dims": tuple(),
              "signal_is_cbed_pattern_set": signal_is_cbed_pattern_set,
              "signal_dtype": "float"}
    blank_postprocessed_2d_signal = _blank_postprocessed_2d_signal(**kwargs)

    num_pixels_in_postprocessed_2d_signal_space = \
        np.prod(blank_postprocessed_2d_signal.data.shape[-2:])

    return num_pixels_in_postprocessed_2d_signal_space



def _blank_postprocessed_2d_signal(sample_specification,
                                   postprocessing_seq,
                                   navigation_dims,
                                   signal_is_cbed_pattern_set,
                                   signal_dtype):
    kwargs = {"sample_specification": sample_specification,
              "navigation_dims": navigation_dims,
              "signal_is_cbed_pattern_set": signal_is_cbed_pattern_set,
              "signal_dtype": signal_dtype}
    blank_unprocessed_2d_signal = _blank_unprocessed_2d_signal(**kwargs)

    kwargs = {"input_signal": blank_unprocessed_2d_signal,
              "postprocessing_seq": postprocessing_seq}
    blank_postprocessed_2d_signal = _postprocess_2d_signal(**kwargs)

    return blank_postprocessed_2d_signal



def _postprocess_2d_signal(input_signal, postprocessing_seq):
    pixel_area = np.abs(input_signal.axes_manager[-2].scale
                        * input_signal.axes_manager[-1].scale)
    input_signal.data /= pixel_area
    
    for idx, postprocessing_step in enumerate(postprocessing_seq):
        optional_params = postprocessing_step
        if isinstance(optional_params, empix.OptionalCroppingParams):
            output_signal = empix.crop(input_signal, optional_params)
        elif isinstance(optional_params, empix.OptionalDownsamplingParams):
            output_signal = empix.downsample(input_signal, optional_params)
        else:
            output_signal = empix.resample(input_signal, optional_params)

        if idx == len(postprocessing_seq)-1:
            postprocess_2d_signal = output_signal
        else:
            input_signal = output_signal

    if len(postprocessing_seq) == 0:
        postprocess_2d_signal = input_signal

    pixel_area = np.abs(postprocess_2d_signal.axes_manager[-2].scale
                        * postprocess_2d_signal.axes_manager[-1].scale)
    postprocess_2d_signal.data *= pixel_area

    return postprocess_2d_signal



def _num_pixels_in_unprocessed_2d_signal_space(sample_specification,
                                               signal_is_cbed_pattern_set):
    kwargs = {"sample_specification": sample_specification,
              "navigation_dims": tuple(),
              "signal_is_cbed_pattern_set": signal_is_cbed_pattern_set,
              "signal_dtype": "float"}
    blank_unprocessed_2d_signal = _blank_unprocessed_2d_signal(**kwargs)

    num_pixels_in_unprocessed_2d_signal_space = \
        np.prod(blank_unprocessed_2d_signal.data.shape[-2:])

    return num_pixels_in_unprocessed_2d_signal_space



###########################
## Define error messages ##
###########################

_check_and_convert_postprocessing_seq_err_msg_1 = \
    ("The object ``postprocessing_seq`` must be a sequence of objects of any "
     "of the following types: (`empix.OptionalCroppingParams`, "
     "`empix.OptionalDownsamplingParams`, `empix.OptionalResamplingParams`).")
_check_and_convert_postprocessing_seq_err_msg_2 = \
    ("The object ``postprocessing_seq[{}]``, being an instance of the class "
     " `empix.OptionalDownsamplingParams`, must satisfy "
     "``postprocessing_seq[{}].core_attrs['downsample_mode']=='mean'`` in this "
     "context.")
