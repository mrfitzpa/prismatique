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
r"""For specifying simulation parameters related to workers, i.e. GPU and CPU
workers.

Note that the documentation in this module draws from Ref. [Pryor1]_, directly 
copying certain passages when convenient.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For validating objects.
import czekitout.check

# For defining classes that support enforced validation, updatability,
# pre-serialization, and de-serialization.
import fancytypes



# Import child modules and packages of current package.
import prismatique.worker.cpu
import prismatique.worker.gpu



##################################
## Define classes and functions ##
##################################

# List of public objects in objects.
__all__ = ["Params"]



def _check_and_convert_cpu_params(params):
    module_alias = prismatique.worker.cpu
    func_alias = module_alias._check_and_convert_cpu_params
    cpu_params = func_alias(params)

    return cpu_params



def _pre_serialize_cpu_params(cpu_params):
    obj_to_pre_serialize = cpu_params
    module_alias = prismatique.worker.cpu
    func_alias = module_alias._pre_serialize_cpu_params
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_cpu_params(serializable_rep):
    module_alias = prismatique.worker.cpu
    func_alias = module_alias._de_pre_serialize_cpu_params
    cpu_params = func_alias(serializable_rep)

    return cpu_params



def _check_and_convert_gpu_params(params):
    module_alias = prismatique.worker.gpu
    func_alias = module_alias._check_and_convert_gpu_params
    gpu_params = func_alias(params)

    return gpu_params



def _pre_serialize_gpu_params(gpu_params):
    obj_to_pre_serialize = gpu_params
    module_alias = prismatique.worker.gpu
    func_alias = module_alias._pre_serialize_gpu_params
    serializable_rep = func_alias(obj_to_pre_serialize)

    return serializable_rep



def _de_pre_serialize_gpu_params(serializable_rep):
    module_alias = prismatique.worker.gpu
    func_alias = module_alias._de_pre_serialize_gpu_params
    gpu_params = func_alias(serializable_rep)

    return gpu_params



_default_cpu_params = \
    None
_default_gpu_params = \
    None
_default_skip_validation_and_conversion = \
    False



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
    ctor_param_names = ("cpu_params",
                        "gpu_params")
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
                 cpu_params=\
                 _default_cpu_params,
                 gpu_params=\
                 _default_gpu_params,
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



def _check_and_convert_worker_params(params):
    obj_name = "worker_params"
    obj = params[obj_name]

    accepted_types = (Params, type(None))
    
    if isinstance(obj, accepted_types[-1]):
        worker_params = accepted_types[0]()
    else:
        kwargs = {"obj": obj,
                  "obj_name": obj_name,
                  "accepted_types": accepted_types}
        czekitout.check.if_instance_of_any_accepted_types(**kwargs)

        kwargs = obj.get_core_attrs(deep_copy=False)
        worker_params = accepted_types[0](**kwargs)

    return worker_params



def _pre_serialize_worker_params(worker_params):
    obj_to_pre_serialize = worker_params
    serializable_rep = obj_to_pre_serialize.pre_serialize()
    
    return serializable_rep



def _de_pre_serialize_worker_params(serializable_rep):
    worker_params = Params.de_pre_serialize(serializable_rep)
    
    return worker_params



_default_worker_params = None



###########################
## Define error messages ##
###########################
