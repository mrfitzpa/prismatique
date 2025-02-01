r"""For specifying output parameters that are applicable only to the multislice
implementation of STEM simulations.

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



def _check_and_convert_num_slices_per_output(ctor_params):
    kwargs = {"obj": ctor_params["num_slices_per_output"],
              "obj_name": "num_slices_per_output"}
    num_slices_per_output = czekitout.convert.to_positive_int(**kwargs)
    
    return num_slices_per_output



def _pre_serialize_num_slices_per_output(num_slices_per_output):
    serializable_rep = num_slices_per_output

    return serializable_rep



def _de_pre_serialize_num_slices_per_output(serializable_rep):
    num_slices_per_output = serializable_rep

    return num_slices_per_output



def _check_and_convert_z_start_output(ctor_params):
    kwargs = {"obj": ctor_params["z_start_output"],
              "obj_name": "z_start_output"}
    z_start_output = czekitout.convert.to_nonnegative_float(**kwargs)
    
    return z_start_output



def _pre_serialize_z_start_output(z_start_output):
    serializable_rep = z_start_output

    return serializable_rep



def _de_pre_serialize_z_start_output(serializable_rep):
    z_start_output = serializable_rep

    return z_start_output



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The output parameters that are applicable only to the multislice
    implementation of STEM simulations.

    For a general discussion on the possible output that can be generated from
    STEM simulations, see the documentation for the class
    :class:`prismatique.stem.output.Params`. That discussion provides important
    context to the description of the parameter set below.

    Parameters
    ----------
    num_slices_per_output : `int`, optional
        In the multislice algorithm, the sample is partitioned into slices along
        the :math:`z`-axis [i.e. the optic axis]. One can optionally save output
        for exit waves that emerge immediately after select intermediate slices,
        in addition to the output for the final exit wave. The output generated
        after a given slice is saved into an "output layer", i.e. a given output
        layer will contain the output for an exit wave that emerged from an
        intermediate slice or the final slice. Note that the exit wave emerging
        from the :math:`n^{\text{th}}` slice, where :math:`n=0` is the first
        slice, emerges at a depth of :math:`\left(n+1\right)\delta z`, where
        :math:`\delta z` is the slice thickness, given by
        Eq. :eq:`slice_thickness_in_potential_params`. ``num_slices_per_output``
        specifies the number of slices between intermediate outputs. For
        example, if ``num_slices_per_output`` is set to ``1``, then the output
        from each slice is saved; if ``num_slices_per_output`` is set to ``2``,
        then the output from every second slice is saved; and so on. Note that
        the output for the exit wave that emerges from the final slice is always
        saved to the last output layer, irrespective of the value of
        ``num_slices_per_output``. Moreover, note that ``num_slices_per_output``
        must be a positive `int`.
    z_start_output : `float`, optional
        Continuing from above, ``z_start_output`` specifies the depth in
        angstroms along the :math:`z`-axis at which intermediate output
        collection begins. Let ``m = int(math.ceil(z_start_output /
        (num_slices_per_output*slice_thickness))`` and ``n = max(1,
        m)*slices_per_output - 1``, where ``slice_thickness`` is the slice
        thickness :math:`\delta z` in angstroms. If ``n<=N_slice-1``, where
        ``N_slice`` is the total number of slices used to partition the sample,
        then intermediate output is stored with the output from the exit wave
        that emerges from the ``n`` th slice being saved into the first output
        layer. Note that the output for the exit wave that emerges from the
        final slice is always saved to the last output layer, irrespective of
        the value of ``z_start_output``. Note that ``z_start_output`` must be a
        nonnegative `float`.

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
        {"num_slices_per_output": _check_and_convert_num_slices_per_output,
         "z_start_output": _check_and_convert_z_start_output}

    _pre_serialization_funcs = \
        {"num_slices_per_output": _pre_serialize_num_slices_per_output,
         "z_start_output": _pre_serialize_z_start_output}

    _de_pre_serialization_funcs = \
        {"num_slices_per_output": _de_pre_serialize_num_slices_per_output,
         "z_start_output": _de_pre_serialize_z_start_output}
    
    def __init__(self,
                 num_slices_per_output=1,
                 z_start_output=float("inf")):
        ctor_params = {"num_slices_per_output": num_slices_per_output,
                       "z_start_output": z_start_output}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_multislice_output_params(ctor_params):
    multislice_output_params = \
        copy.deepcopy(ctor_params["multislice_output_params"])
    if multislice_output_params is None:
        multislice_output_params = Params()
        
    kwargs = {"obj": multislice_output_params,
              "obj_name": "multislice_output_params",
              "accepted_types": (Params, type(None))}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return multislice_output_params



def _pre_serialize_multislice_output_params(multislice_output_params):
    serializable_rep = multislice_output_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_multislice_output_params(serializable_rep):
    multislice_output_params = Params.de_pre_serialize(serializable_rep)

    return multislice_output_params
    


###########################
## Define error messages ##
###########################
