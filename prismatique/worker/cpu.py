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
r"""For specifying simulation parameters related to CPU workers.

Note that the documentation in this module draws from reference [Pryor1]_.
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



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params"]



def _check_and_convert_enable_workers(params):
    obj_name = "enable_workers"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    enable_workers = czekitout.convert.to_bool(**kwargs)

    return enable_workers



def _pre_serialize_enable_workers(enable_workers):
    obj_to_pre_serialize = enable_workers
    serializable_rep = obj_to_pre_serialize
    
    return serializable_rep



def _de_pre_serialize_enable_workers(serializable_rep):
    enable_workers = serializable_rep

    return enable_workers



def _check_and_convert_num_worker_threads(params):
    obj_name = "num_worker_threads"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    num_worker_threads = czekitout.convert.to_nonnegative_int(**kwargs)

    return num_worker_threads



def _pre_serialize_num_worker_threads(num_worker_threads):
    obj_to_pre_serialize = num_worker_threads
    serializable_rep = obj_to_pre_serialize
    
    return serializable_rep



def _de_pre_serialize_num_worker_threads(serializable_rep):
    num_worker_threads = serializable_rep

    return num_worker_threads



def _check_and_convert_batch_size(params):
    obj_name = "batch_size"
    kwargs = {"obj": params[obj_name], "obj_name": obj_name}
    batch_size = czekitout.convert.to_positive_int(**kwargs)

    return batch_size



def _pre_serialize_batch_size(batch_size):
    obj_to_pre_serialize = batch_size
    serializable_rep = obj_to_pre_serialize
    
    return serializable_rep



def _de_pre_serialize_batch_size(serializable_rep):
    batch_size = serializable_rep

    return batch_size



def _check_and_convert_early_stop_count(params):
    obj_name = "early_stop_count"

    current_func_name = "_check_and_convert_early_stop_count"
    
    try:
        kwargs = {"obj": params[obj_name], "obj_name": obj_name}
        early_stop_count = czekitout.convert.to_positive_int(**kwargs)
    except:
        err_msg = globals()[current_func_name+"_err_msg_1"]
        raise TypeError(err_msg)
    
    if early_stop_count < 2:
        err_msg = globals()[current_func_name+"_err_msg_1"]
        raise ValueError(err_msg)
    
    return early_stop_count



def _pre_serialize_early_stop_count(early_stop_count):
    obj_to_pre_serialize = early_stop_count
    serializable_rep = obj_to_pre_serialize

    return serializable_rep



def _de_pre_serialize_early_stop_count(serializable_rep):
    early_stop_count = serializable_rep

    return early_stop_count



_default_enable_workers = True
_default_num_worker_threads = 12
_default_batch_size = 1
_default_early_stop_count = 100
_default_skip_validation_and_conversion = False



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related to CPU workers.

    Parameters
    ----------
    enable_workers : `bool`, optional
        If set to ``True``, then the simulation will also make use of CPU 
        workers, in addition to GPU workers (if available).
    num_worker_threads : `int`, optional
        Let ``num_available_worker_threads`` be the number of CPU worker threads
        available for the simulation. If ``enable_workers`` has been set to
        ``True``, then ``max(num_available_worker_threads, num_worker_threads)``
        determines the number of CPU worker threads that are to be used in the
        simulation. If ``enable_workers`` has been set to ``False``, then the
        parameter ``num_worker_threads`` is ignored upon configuring the
        simulation. See the documentation for the class
        :class:`prismatique.worker.Params` for a discussion on how the number of
        CPU worker threads affects performance.
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

        If ``enable_workers`` has been set to ``True``, then ``batch_size`` 
        specifies the number of probes or plane waves to transmit simultaneously
        per CPU worker. If ``enable_workers`` has been set to ``False``, then 
        the parameter ``batch_size`` is ignored upon configuring the simulation. 
    early_stop_count : `int`, optional
        Assuming that GPUs have been enabled to do work in the simulation, then
        the work dispatcher will cease providing work to the CPU workers
        ``early_stop_count`` jobs from the end. This is to prevent the program
        waiting for CPU workers to complete that are slower than the enabled
        GPU workers. What qualifies as a job in this context is unclear from
        the documentation of the ``prismatic`` library. Presumably, a typical
        job could be a batch FFT operation, or a series of element-wise
        multiplication operations. In any case, the user could try experimenting
        with this parameter to see if it yields appreciable performance gains.
        If no GPUs have been enabled to do  work in the simulation, then the
        parameter ``early_stop_count`` is ignored upon configuring the 
        simulation.
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
    ctor_param_names = ("enable_workers",
                        "num_worker_threads",
                        "batch_size",
                        "early_stop_count")
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
                 enable_workers=\
                 _default_enable_workers,
                 num_worker_threads=\
                 _default_num_worker_threads,
                 batch_size=\
                 _default_batch_size,
                 early_stop_count=\
                 _default_early_stop_count,
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



def _check_and_convert_cpu_params(params):
    obj_name = "cpu_params"
    obj = params[obj_name]

    accepted_types = (Params, type(None))
    
    if isinstance(obj, accepted_types[-1]):
        cpu_params = accepted_types[0]()
    else:
        kwargs = {"obj": obj,
                  "obj_name": obj_name,
                  "accepted_types": accepted_types}
        czekitout.check.if_instance_of_any_accepted_types(**kwargs)

        kwargs = obj.get_core_attrs(deep_copy=False)
        cpu_params = accepted_types[0](**kwargs)

    return cpu_params



def _pre_serialize_cpu_params(cpu_params):
    obj_to_pre_serialize = cpu_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_cpu_params(serializable_rep):
    cpu_params = Params.de_pre_serialize(serializable_rep)
    
    return cpu_params



_default_cpu_params = None
    


###########################
## Define error messages ##
###########################

_check_and_convert_early_stop_count_err_msg_1 = \
    ("The object ``early_stop_count`` must be an integer greater than 1.")
