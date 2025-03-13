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
r"""For specifying simulation parameters related to GPU workers.

Note that the documentation in this module draws from Refs. [Pryor1]_ and
[Hinitt1]_, directly copying certain passages when convenient.
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



# For recycling helper functions and/or constants.
import prismatique.worker.cpu



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params"]



def _check_and_convert_num_gpus(params):
    obj_name = "num_gpus"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    num_gpus = czekitout.convert.to_nonnegative_int(**kwargs)

    return num_gpus



def _pre_serialize_num_gpus(num_gpus):
    obj_to_pre_serialize = num_gpus
    serializable_rep = obj_to_pre_serialize
    
    return serializable_rep



def _de_pre_serialize_num_gpus(serializable_rep):
    num_gpus = serializable_rep

    return num_gpus



def _check_and_convert_batch_size(params):
    module_alias = prismatique.worker.cpu
    func_alias = module_alias._check_and_convert_batch_size
    batch_size = func_alias(params)

    return batch_size



def _pre_serialize_batch_size(batch_size):
    obj_to_pre_serialize = batch_size
    module_alias = prismatique.worker.cpu
    func_alias = module_alias._pre_serialize_batch_size
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_batch_size(serializable_rep):
    module_alias = prismatique.worker.cpu
    func_alias = module_alias._de_pre_serialize_batch_size
    batch_size = func_alias(serializable_rep)

    return batch_size



def _check_and_convert_data_transfer_mode(params):
    obj_name = "data_transfer_mode"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    obj = czekitout.convert.to_str_from_str_like(**kwargs)

    kwargs = {"obj": obj,
              "obj_name": obj_name,
              "accepted_strings": ("single-transfer", "streaming", "auto")}
    czekitout.check.if_one_of_any_accepted_strings(**kwargs)

    data_transfer_mode = obj
    
    return data_transfer_mode



def _pre_serialize_data_transfer_mode(data_transfer_mode):
    obj_to_pre_serialize = data_transfer_mode
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_data_transfer_mode(serializable_rep):
    data_transfer_mode = serializable_rep

    return data_transfer_mode



def _check_and_convert_num_streams_per_gpu(params):
    obj_name = "num_streams_per_gpu"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    num_streams_per_gpu = czekitout.convert.to_positive_int(**kwargs)

    return num_streams_per_gpu



def _pre_serialize_num_streams_per_gpu(num_streams_per_gpu):
    obj_to_pre_serialize = num_streams_per_gpu
    serializable_rep = obj_to_pre_serialize
    
    return serializable_rep



def _de_pre_serialize_num_streams_per_gpu(serializable_rep):
    num_streams_per_gpu = serializable_rep

    return num_streams_per_gpu



_default_num_gpus = \
    4
_default_batch_size = \
    1
_default_data_transfer_mode = \
    "auto"
_default_num_streams_per_gpu = \
    3
_default_skip_validation_and_conversion = \
    False



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related to GPU workers.

    Parameters
    ----------
    num_gpus : `int`, optional
        Let ``num_available_gpus`` be the number of GPU devices available for
        the simulation. ``max(num_available_gpus, num_gpus)`` determines the
        number of GPUs that are to be used in the simulation. See the
        documentation for the class :class:`prismatique.worker.Params` for a 
        discussion on how the number of GPU devices affects performance.
    batch_size : `int`, optional
        The calculation of the transmission of a single probe or plane wave 
        through the entire sample (i.e. from the incident surface to the exit 
        surface) involves a series of fast-Fourier transform (FFT) operations. 
        FFTs are calculated using a divide-and-conquer algorithm that 
        recursively breaks down a discrete Fourier transform (DFT) into smaller 
        DFTs and performs multiplications involving complex roots of unity 
        called twiddle factors. Thus, a given FFT in this scheme is calculated 
        in multiple steps. The libraries used in ``prismatic`` that implement 
        FFTs support batch FFTs, whereby multiple Fourier transforms of the same
        size can be computed simultaneously. By simultaneously, we mean that 
        step :math:`i+1` of a given FFT in a given batch cannot be executed 
        until step :math:`i` has been executed for all FFTs in said batch. This 
        order of operations allows for reuse of intermediate twiddle factors, 
        resulting in a faster overall computation than performing individual
        transforms one-by-one at the expense of requiring a larger block of
        memory to store the multiple arrays. We can therefore use this batch
        FFT method to calculate the transmission of a batch of probes or plane
        waves simultaneously in the same sense as that articulated above.

        If ``num_gpus`` has been set to a positive integer`, then ``batch_size``
        specifies the number of probes or plane waves to transmit simultaneously
        per GPU device. If ``num_gpus`` has been set to ``0``, then the
        parameter ``batch_size`` is ignored upon configuring the simulation. 
    data_transfer_mode : ``"single-transfer"`` | ``"streaming"`` | ``"auto"``, optional
        The preferred way to perform simulations is to transfer large data
        structures such as the projected potential array or the compact 
        scattering matrices to each GPU only once, where they can then be read
        from repeatedly over the course of the calculation. However, this
        requires that the arrays fit into the limited GPU memory. For 
        simulations that are too large, ``prismatic`` has implemented an
        asynchronous streaming version for simulations. A stream is a sequence
        of operations which are processed in order; however, different streams
        can execute out of order with respect to one another. These operations
        include kernel executions and memory transfers. Each GPU device
        can manage multiple streams, where each stream may use some subset of
        the threads in said GPU device. Since only one kernel is able to run on
        a given GPU device at any one time, a queue of streams can be formed
        such that the memory copies of one stream can overlap with the kernel
        execution of another stream as depicted in 
        :numref:`worker_gpu_params_illustrating_streaming`.

        .. _worker_gpu_params_illustrating_streaming:
        .. figure:: ../_images/illustrating_streaming.png

           Depiction of streaming execution. Figure taken from Ref. [Hinitt1]_.

        In streaming mode, rather than allocate and transfer a single read-only 
        copy of large arrays, buffers are allocated to each stream large enough 
        to hold only the relevant subset of the data for the current step in the
        calculation, and the job itself triggers asynchronous streaming of the 
        data it requires for the next step. The use of asynchronous memory 
        copies and CUDA streams permits the partial hiding of memory transfer 
        latencies behind kernel execution.

        By default, ``data_transfer_mode`` is set to ``"auto"``, which signals
        ``prismatic`` to use an automatic procedure to determine whether to use
        the single-transfer or streaming mode, whereby the input parameters are
        used to estimate how much memory will be consumed on the device, and if
        this estimate is too large compared with the available device memory
        then the streaming mode is used. Users can manually select streaming
        mode by setting ``data_transfer_mode`` to ``"streaming"``, or if memory
        permits so, users can also manually select single-transfer mode by
        setting ``data_transfer_mode`` to ``"single-transfer"``. If ``num_gpus``
        has been set to ``0``, then the parameter ``data_transfer_mode`` is 
        ignored upon configuring the simulation.
    num_streams_per_gpu : `int`, optional
        If ``num_gpus`` has been set to a positive integer` and streaming
        mode has been enabled, then ``num_streams_per_gpu`` specifies the number
        of CUDA streams per GPU device. If ``num_gpus`` has been set to ``0`` or
        streaming mode has not been enabled, then the parameter 
        ``num_streams_per_gpu`` is ignored upon configuring the simulation. 
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
    ctor_param_names = ("num_gpus",
                        "batch_size",
                        "data_transfer_mode",
                        "num_streams_per_gpu")
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
                 num_gpus=\
                 _default_num_gpus,
                 batch_size=\
                 _default_batch_size,
                 data_transfer_mode=\
                 _default_data_transfer_mode,
                 num_streams_per_gpu=\
                 _default_num_streams_per_gpu,
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



def _check_and_convert_gpu_params(params):
    obj_name = "gpu_params"
    obj = params[obj_name]

    accepted_types = (Params, type(None))
    
    if isinstance(obj, accepted_types[-1]):
        gpu_params = accepted_types[0]()
    else:
        kwargs = {"obj": obj,
                  "obj_name": obj_name,
                  "accepted_types": accepted_types}
        czekitout.check.if_instance_of_any_accepted_types(**kwargs)

        kwargs = obj.get_core_attrs(deep_copy=False)
        gpu_params = accepted_types[0](**kwargs)

    return gpu_params



def _pre_serialize_gpu_params(gpu_params):
    obj_to_pre_serialize = gpu_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_gpu_params(serializable_rep):
    gpu_params = Params.de_pre_serialize(serializable_rep)
    
    return gpu_params



_default_gpu_params = None
    


###########################
## Define error messages ##
###########################
