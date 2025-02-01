r"""For specifying simulation parameters related to GPU workers.

Note that the documentation in this module draws from Refs. [Pryor1]_ and
[Hinitt1]_, directly copying certain passages when convenient.
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



def _check_and_convert_num_gpus(ctor_params):
    kwargs = {"obj": ctor_params["num_gpus"],
              "obj_name": "num_gpus"}
    num_gpus = czekitout.convert.to_nonnegative_int(**kwargs)
    
    return num_gpus



def _pre_serialize_num_gpus(num_gpus):
    serializable_rep = num_gpus

    return serializable_rep



def _de_pre_serialize_num_gpus(serializable_rep):
    num_gpus = serializable_rep

    return num_gpus



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



def _check_and_convert_data_transfer_mode(ctor_params):
    try:
        kwargs = {"obj": ctor_params["data_transfer_mode"],
                  "obj_name": "data_transfer_mode"}
        data_transfer_mode = czekitout.convert.to_str_from_str_like(**kwargs)
    except:
        raise TypeError(_check_and_convert_data_transfer_mode_err_msg_1)

    accepted_values = ("single-transfer", "streaming", "auto")
    if data_transfer_mode not in accepted_values:
        raise ValueError(_check_and_convert_data_transfer_mode_err_msg_1)
    
    return data_transfer_mode



def _pre_serialize_data_transfer_mode(data_transfer_mode):
    serializable_rep = data_transfer_mode

    return serializable_rep



def _de_pre_serialize_data_transfer_mode(serializable_rep):
    data_transfer_mode = serializable_rep

    return data_transfer_mode



def _check_and_convert_num_streams_per_gpu(ctor_params):
    kwargs = {"obj": ctor_params["num_streams_per_gpu"],
              "obj_name": "num_streams_per_gpu"}
    num_streams_per_gpu = czekitout.convert.to_positive_int(**kwargs)
    
    return num_streams_per_gpu



def _pre_serialize_num_streams_per_gpu(num_streams_per_gpu):
    serializable_rep = num_streams_per_gpu

    return serializable_rep



def _de_pre_serialize_num_streams_per_gpu(serializable_rep):
    num_streams_per_gpu = serializable_rep

    return num_streams_per_gpu



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
        {"num_gpus": _check_and_convert_num_gpus,
         "batch_size": _check_and_convert_batch_size,
         "data_transfer_mode": _check_and_convert_data_transfer_mode,
         "num_streams_per_gpu": _check_and_convert_num_streams_per_gpu}

    _pre_serialization_funcs = \
        {"num_gpus": _pre_serialize_num_gpus,
         "batch_size": _pre_serialize_batch_size,
         "data_transfer_mode": _pre_serialize_data_transfer_mode,
         "num_streams_per_gpu": _pre_serialize_num_streams_per_gpu}

    _de_pre_serialization_funcs = \
        {"num_gpus": _de_pre_serialize_num_gpus,
         "batch_size": _de_pre_serialize_batch_size,
         "data_transfer_mode": _de_pre_serialize_data_transfer_mode,
         "num_streams_per_gpu": _de_pre_serialize_num_streams_per_gpu}
    
    def __init__(self,
                 num_gpus=4,
                 batch_size=1,
                 data_transfer_mode="auto",
                 num_streams_per_gpu=3):
        ctor_params = {"num_gpus": num_gpus,
                       "batch_size": batch_size,
                       "data_transfer_mode": data_transfer_mode,
                       "num_streams_per_gpu": num_streams_per_gpu}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_gpu_params(ctor_params):
    gpu_params = copy.deepcopy(ctor_params["gpu_params"])
    if gpu_params is None:
        gpu_params = Params()
    
    kwargs = {"obj": gpu_params,
              "obj_name": "gpu_params",
              "accepted_types": (Params, type(None))}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return gpu_params



def _pre_serialize_gpu_params(gpu_params):
    serializable_rep = gpu_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_gpu_params(serializable_rep):
    gpu_params = Params.de_pre_serialize(serializable_rep)

    return gpu_params
    


###########################
## Define error messages ##
###########################

_check_and_convert_data_transfer_mode_err_msg_1 = \
    ("The object ``data_transfer_mode`` must be one of three strings: "
     "``'single-transfer'``, ``'streaming'``, or ``'auto'``.")