r"""For specifying simulation parameters related to workers, i.e. GPU and CPU
workers.

Note that the documentation in this module draws from Ref. [Pryor1]_, directly 
copying certain passages when convenient.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For performing deep copies of objects.
import copy



# For validating objects.
import czekitout.check

# For defining classes that support enforced validation, updatability,
# pre-serialization, and de-serialization.
import fancytypes



# Import child modules and packages of current package.
import prismatique.worker.cpu
import prismatique.worker.gpu



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



def _check_and_convert_cpu_params(ctor_params):
    cpu_params = \
        prismatique.worker.cpu._check_and_convert_cpu_params(ctor_params)

    return cpu_params



def _pre_serialize_cpu_params(cpu_params):
    serializable_rep = \
        prismatique.worker.cpu._pre_serialize_cpu_params(cpu_params)

    return serializable_rep



def _de_pre_serialize_cpu_params(serializable_rep):
    cpu_params = \
        prismatique.worker.cpu._de_pre_serialize_cpu_params(serializable_rep)

    return cpu_params



def _check_and_convert_gpu_params(ctor_params):
    gpu_params = \
        prismatique.worker.gpu._check_and_convert_gpu_params(ctor_params)

    return gpu_params



def _pre_serialize_gpu_params(gpu_params):
    serializable_rep = \
        prismatique.worker.gpu._pre_serialize_gpu_params(gpu_params)

    return serializable_rep



def _de_pre_serialize_gpu_params(serializable_rep):
    gpu_params = \
        prismatique.worker.gpu._de_pre_serialize_gpu_params(serializable_rep)

    return gpu_params



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related to GPU and CPU workers.

    In Ref. [Pryor1]_, the authors compared the performance of the multislice
    and PRISM algorithms [as implemented in ``prismatic``]. They also studied
    how performance scales with the number of CPU worker threads and GPU devices
    used. :numref:`worker_params_gpu_and_cpu_scaling`, which was taken from Ref.
    Pryor1]_, presents benchmarking results for STEM simulations of amorphous
    carbon. From this figure, we see that when using only CPU worker threads,
    the wall time for both algorithms scales approximately like
    :math:`1/N_{\text{CPU}}`, where :math:`N_{\text{CPU}}` is the number of CPU
    worker threads used during the simulation. In other words, if one doubles
    the number of CPU worker threads used, then one should expect approximately
    half the wall time. We also see from
    :numref:`worker_params_gpu_and_cpu_scaling` that the addition of a single
    GPU device improves both algorithms by approximately a factor of 8 in this
    case, but in general, the relative improvement varies depending on the
    quality and number of the CPUs vs GPUs. The addition of a second GPU device
    roughly doubles the performance, and then doubling the number of GPU devices
    again to 4 roughly doubles the performance. Users should use these 
    benchmarking results as a guide when deciding how many CPU worker threads
    and GPU devices to use in their simulations. That being said, users are
    recommended to do their own benchmarking for different simulation cases.

    .. _worker_params_gpu_and_cpu_scaling:
    .. figure:: ../_images/gpu_and_cpu_performance_scaling.png

       Comparison of the implementations of multislice and PRISM for varying
       combinations of CPU threads and GPUs. The simulation was performed on a
       :math:`100 \times 100 \times 100\ \text{Å}^3` amorphous carbon cell with
       :math:`5\ \text{Å}` thick slices, a :math:`0.1\ \text{Å}` pixel size, and
       a 20 mrad probe convergence semi-angle. All simulations were performed on
       compute nodes with dual Intel Xeon E5-2650 processors, four Tesla K20
       GPUs, and 64 GB RAM. Calculation time of rightmost data point is labeled
       for all curves. Figure taken from Ref. [Pryor1]_.

    Parameters
    ----------
    cpu_params : :class:`prismatique.worker.cpu.Params` | `None`, optional
        The simulation parameters related to CPU workers. If ``cpu`` is set to
        `None` [i.e. the default value], then the simulation parameters related
        to CPU workers are set to default values.
    gpu_params : :class:`prismatique.worker.gpu.Params` | `None`, optional
        The simulation parameters related to GPU workers. If ``gpu`` is set to
        `None` [i.e. the default value], then the simulation parameters related
        to GPU workers are set to default values.

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
        {"cpu_params": _check_and_convert_cpu_params,
         "gpu_params": _check_and_convert_gpu_params}

    _pre_serialization_funcs = \
        {"cpu_params": _pre_serialize_cpu_params,
         "gpu_params": _pre_serialize_gpu_params}

    _de_pre_serialization_funcs = \
        {"cpu_params": _de_pre_serialize_cpu_params,
         "gpu_params": _de_pre_serialize_gpu_params}
    
    def __init__(self, cpu_params=None, gpu_params=None):
        ctor_params = {"cpu_params": cpu_params, "gpu_params": gpu_params}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_worker_params(ctor_params):
    worker_params = copy.deepcopy(ctor_params["worker_params"])
    if worker_params is None:
        worker_params = Params()
    
    kwargs = {"obj": worker_params,
              "obj_name": "worker_params",
              "accepted_types": (Params, type(None))}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return worker_params



def _pre_serialize_worker_params(worker_params):
    serializable_rep = worker_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_worker_params(serializable_rep):
    worker_params = Params.de_pre_serialize(serializable_rep)

    return worker_params



###########################
## Define error messages ##
###########################