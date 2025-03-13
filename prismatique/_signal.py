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
r"""Contains helper functions for postprocessing CBED intensity patterns and
HRTEM intensity images.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For general array handling.
import numpy as np

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



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = []



def _check_and_convert_postprocessing_seq(params):
    obj_name = "postprocessing_seq"

    current_func_name = "_check_and_convert_postprocessing_seq"

    try:
        obj = params[obj_name]
        accepted_types = (empix.OptionalCroppingParams,
                          empix.OptionalDownsamplingParams,
                          empix.OptionalResamplingParams)

        for obj_elem in obj:
            kwargs = {"obj": obj_elem,
                      "obj_name": "postprocessing_step",
                      "accepted_types": accepted_types}
            czekitout.check.if_instance_of_any_accepted_types(**kwargs)
    except:
        err_msg = globals()[current_func_name+"_err_msg_1"]
        raise TypeError(err_msg)

    postprocessing_seq = tuple()
    for obj_elem_idx, obj_elem in enumerate(obj):
        obj_elem_core_attrs = obj_elem.get_core_attrs(deep_copy=False)
        if isinstance(obj_elem, empix.OptionalDownsamplingParams):
            if obj_elem_core_attrs["downsample_mode"] != "mean":
                unformatted_err_msg = globals()[current_func_name+"_err_msg_2"]
                err_msg = unformatted_err_msg.format(obj_elem_idx, obj_elem_idx)
                raise ValueError(err_msg)

        kwargs = obj_elem_core_attrs
        postprocessing_step = type(obj_elem)(**kwargs)

        postprocessing_seq += (postprocessing_step,)
    
    return postprocessing_seq



def _pre_serialize_postprocessing_seq(postprocessing_seq):
    obj_to_pre_serialize = postprocessing_seq
    serializable_rep = tuple()
    for elem_of_obj_to_pre_serialize in obj_to_pre_serialize:
        serializable_rep += (elem_of_obj_to_pre_serialize.pre_serialize(),)

    return serializable_rep



def _de_pre_serialize_postprocessing_seq(serializable_rep):
    postprocessing_seq = tuple()
    
    for elem_of_serializable_rep in serializable_rep:
        if "center" in elem_of_serializable_rep:
            cls_alias = empix.OptionalCroppingParams
        elif "block_dims" in elem_of_serializable_rep:
            cls_alias = empix.OptionalDownsamplingParams
        else:
            cls_alias = empix.OptionalResamplingParams
            
        postprocessing_step = \
            cls_alias.de_pre_serialize(elem_of_serializable_rep)
        postprocessing_seq += \
            (postprocessing_step,)

    return postprocessing_seq



def _blank_unprocessed_2d_signal(sample_specification,
                                 navigation_dims,
                                 signal_is_cbed_pattern_set,
                                 signal_dtype):
    kwargs = {"sample_specification": sample_specification,
              "signal_is_cbed_pattern_set": signal_is_cbed_pattern_set}
    signal_space_axes = _signal_space_axes_of_unprocessed_2d_signal(**kwargs)
    
    signal_shape = (navigation_dims
                    + (signal_space_axes[1].size, signal_space_axes[0].size))
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



def _signal_space_axes_of_unprocessed_2d_signal(
        sample_specification,
        signal_is_cbed_pattern_set):
    func_alias = (_k_x_and_k_y_axes_of_unprocessed_dp_signal
                  if signal_is_cbed_pattern_set
                  else _r_x_and_r_y_axes_of_unprocessed_image_signal)
    signal_space_axes = func_alias(sample_specification)

    return signal_space_axes



def _k_x_and_k_y_axes_of_unprocessed_dp_signal(sample_specification):
    k_x, k_y = \
        _k_x_and_k_y_of_unprocessed_dp_signal(sample_specification)

    axes_labels = (r"$k_x$", r"$k_y$")
    sizes = (len(k_x), len(k_y))
    scales = (k_x[1]-k_x[0], -(k_y[-2]-k_y[-1]))
    offsets = (k_x[0], k_y[0])
    units = ("1/Å", "1/Å")

    tol = 10*np.finfo(np.float32).eps

    scales = ((scales[0], -scales[0])
              if (abs(scales[0]+scales[1]) < tol)
              else scales)
    
    k_x_and_k_y_axes = tuple()
    for axis_idx in range(len(units)):
        axis = hyperspy.axes.UniformDataAxis(size=sizes[axis_idx],
                                             scale=scales[axis_idx],
                                             offset=offsets[axis_idx])
        axis.name = axes_labels[axis_idx]
        axis.units = units[axis_idx]
        k_x_and_k_y_axes += (axis,)
    k_x_axis, k_y_axis = k_x_and_k_y_axes

    return k_x_axis, k_y_axis



def _k_x_and_k_y_of_unprocessed_dp_signal(sample_specification):
    module_alias = \
        prismatique.sample
    kwargs = \
        {"sample_specification": sample_specification}
    sample_supercell_xy_dims_in_pixels = \
        module_alias._supercell_xy_dims_in_pixels(**kwargs)
    sample_supercell_lateral_pixel_size = \
        module_alias._supercell_lateral_pixel_size(**kwargs)
    f_x, f_y = \
        module_alias._interpolation_factors_from_sample_specification(**kwargs)

    # Fourier coordinates before anti-aliasing crop.
    k_x = module_alias._FFT_1D_freqs(sample_supercell_xy_dims_in_pixels[0],
                                     sample_supercell_lateral_pixel_size[0])
    k_y = module_alias._FFT_1D_freqs(sample_supercell_xy_dims_in_pixels[1],
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

    tol = 10*np.finfo(np.float32).eps

    scales = ((scales[0], -scales[0])
              if (abs(scales[0]+scales[1]) < tol)
              else scales)

    r_x_and_r_y_axes = tuple()
    for axis_idx in range(len(units)):
        axis = hyperspy.axes.UniformDataAxis(size=sizes[axis_idx],
                                             scale=scales[axis_idx],
                                             offset=offsets[axis_idx])
        axis.name = axes_labels[axis_idx]
        axis.units = units[axis_idx]
        r_x_and_r_y_axes += (axis,)
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



def _num_pixels_in_postprocessed_2d_signal_space(
        sample_specification,
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
