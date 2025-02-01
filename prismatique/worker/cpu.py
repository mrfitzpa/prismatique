r"""For specifying simulation parameters related to CPU workers.

Note that the documentation in this module draws from reference [Pryor1]_.
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



def _check_and_convert_enable_workers(ctor_params):
    kwargs = {"obj": ctor_params["enable_workers"],
              "obj_name": "enable_workers"}
    enable_workers = czekitout.convert.to_bool(**kwargs)
    
    return enable_workers



def _pre_serialize_enable_workers(enable_workers):
    serializable_rep = enable_workers

    return serializable_rep



def _de_pre_serialize_enable_workers(serializable_rep):
    enable_workers = serializable_rep

    return enable_workers



def _check_and_convert_num_worker_threads(ctor_params):
    kwargs = {"obj": ctor_params["num_worker_threads"],
              "obj_name": "num_worker_threads"}
    num_worker_threads = czekitout.convert.to_nonnegative_int(**kwargs)
    
    return num_worker_threads



def _pre_serialize_num_worker_threads(num_worker_threads):
    serializable_rep = num_worker_threads

    return serializable_rep



def _de_pre_serialize_num_worker_threads(serializable_rep):
    num_worker_threads = serializable_rep

    return num_worker_threads



def _check_and_convert_batch_size(ctor_params):
    kwargs = {"obj": ctor_params["batch_size"],
              "obj_name": "batch_size"}
    batch_size = czekitout.convert.to_positive_int(**kwargs)
    
    return batch_size



def _pre_serialize_batch_size(batch_size):
    serializable_rep = batch_size

    return serializable_rep



def _de_pre_serialize_batch_size(serializable_rep):
    batch_size = serializable_rep

    return batch_size



def _check_and_convert_early_stop_count(ctor_params):
    try:
        kwargs = {"obj": ctor_params["early_stop_count"],
                  "obj_name": "early_stop_count"}
        early_stop_count = czekitout.convert.to_int(**kwargs)
    except:
        raise TypeError(_check_and_convert_early_stop_count_err_msg_1)
    
    if early_stop_count < 2:
        raise ValueError(_check_and_convert_early_stop_count_err_msg_1)
    
    return early_stop_count



def _pre_serialize_early_stop_count(early_stop_count):
    serializable_rep = early_stop_count

    return serializable_rep



def _de_pre_serialize_early_stop_count(serializable_rep):
    early_stop_count = serializable_rep

    return early_stop_count



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
        {"enable_workers": _check_and_convert_enable_workers,
         "num_worker_threads": _check_and_convert_num_worker_threads,
         "batch_size": _check_and_convert_batch_size,
         "early_stop_count": _check_and_convert_early_stop_count}

    _pre_serialization_funcs = \
        {"enable_workers": _pre_serialize_enable_workers,
         "num_worker_threads": _pre_serialize_num_worker_threads,
         "batch_size": _pre_serialize_batch_size,
         "early_stop_count": _pre_serialize_early_stop_count}

    _de_pre_serialization_funcs = \
        {"enable_workers": _de_pre_serialize_enable_workers,
         "num_worker_threads": _de_pre_serialize_num_worker_threads,
         "batch_size": _de_pre_serialize_batch_size,
         "early_stop_count": _de_pre_serialize_early_stop_count}
    
    def __init__(self,
                 enable_workers=True,
                 num_worker_threads=12,
                 batch_size=1,
                 early_stop_count=100):
        ctor_params = {"enable_workers": enable_workers,
                       "num_worker_threads": num_worker_threads,
                       "batch_size": batch_size,
                       "early_stop_count": early_stop_count}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_cpu_params(ctor_params):
    cpu_params = copy.deepcopy(ctor_params["cpu_params"])
    if cpu_params is None:
        cpu_params = Params()
    
    kwargs = {"obj": cpu_params,
              "obj_name": "cpu_params",
              "accepted_types": (Params, type(None))}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return cpu_params



def _pre_serialize_cpu_params(cpu_params):
    serializable_rep = cpu_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_cpu_params(serializable_rep):
    cpu_params = Params.de_pre_serialize(serializable_rep)

    return cpu_params
    


###########################
## Define error messages ##
###########################

_check_and_convert_early_stop_count_err_msg_1 = \
    ("The object ``early_stop_count`` must be an integer greater than 1.")