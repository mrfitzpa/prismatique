r"""For specifying simulation parameters related to the modelling of STEM
systems.

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

# For validating, pre-serializing, and de-pre-serializing instances of the
# class :class:`embeam.stem.probe.ModelParams`.
import embeam.stem.probe



# For validating, pre-serializing, and de-pre-serializing instances of the
# classes :class:`prismatique.sample.ModelParams`,
# :class:`prismatique.sample.PotentialSliceSubsetIDs`,
# :class:`prismatique.sample.SMatrixSubsetIDs`, and
# :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`.
import prismatique.sample

# For validating instances of the :class:`prismatique.scan.rectangular.Params`,
# and for generating probe positions.
import prismatique.scan



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



def _check_and_convert_sample_specification(ctor_params):
    sample_specification = copy.deepcopy(ctor_params["sample_specification"])
    
    accepted_types = (prismatique.sample.ModelParams,
                      prismatique.sample.PotentialSliceSubsetIDs,
                      prismatique.sample.SMatrixSubsetIDs,
                      prismatique.sample.PotentialSliceAndSMatrixSubsetIDs)
    kwargs = {"obj": sample_specification,
              "obj_name": "sample_specification",
              "accepted_types": accepted_types}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)
    
    return sample_specification



def _pre_serialize_sample_specification(sample_specification):
    mod_alias = prismatique.sample
    if "thermal_params" in sample_specification.core_attrs:
        pre_serialize_sample_specification = \
            mod_alias._pre_serialize_sample_model_params
    elif "interpolation_factors" in sample_specification.core_attrs:
        pre_serialize_sample_specification = \
            mod_alias._pre_serialize_potential_slice_subset_ids
    elif "S_matrix_set_filenames" in sample_specification.core_attrs:
        pre_serialize_sample_specification = \
            mod_alias._pre_serialize_potential_slice_set_and_S_matrix_subset_ids
    else:
        pre_serialize_sample_specification = \
            mod_alias._pre_serialize_S_matrix_subset_ids
        
    serializable_rep = pre_serialize_sample_specification(sample_specification)

    return serializable_rep



def _de_pre_serialize_sample_specification(serializable_rep):
    if "thermal_params" in serializable_rep:
        de_pre_serialize_sample_specification = \
            prismatique.sample._de_pre_serialize_sample_model_params
    elif "interpolation_factors" in serializable_rep:
        de_pre_serialize_sample_specification = \
            prismatique.sample._de_pre_serialize_potential_slice_subset_ids
    elif "S_matrix_set_filenames" in serializable_rep:
        alias = prismatique.sample  # To maintain 80 character line width below.
        de_pre_serialize_sample_specification = \
            alias._de_pre_serialize_potential_slice_set_and_S_matrix_subset_ids
    else:
        de_pre_serialize_sample_specification = \
            prismatique.sample._de_pre_serialize_S_matrix_subset_ids
        
    sample_specification = \
        de_pre_serialize_sample_specification(serializable_rep)

    return sample_specification



def _check_and_convert_probe_model_params(ctor_params):
    probe_model_params = \
        embeam.stem.probe._check_and_convert_probe_model_params(ctor_params)
    
    return probe_model_params



def _pre_serialize_probe_model_params(probe_model_params):
    pre_serialize_probe_model_params = \
        embeam.stem.probe._pre_serialize_probe_model_params
    serializable_rep = \
        pre_serialize_probe_model_params(probe_model_params)

    return serializable_rep



def _de_pre_serialize_probe_model_params(serializable_rep):
    de_pre_serialize_probe_model_params = \
        embeam.stem.probe._de_pre_serialize_probe_model_params
    probe_model_params = \
        de_pre_serialize_probe_model_params(serializable_rep)

    return probe_model_params



def _check_and_convert_specimen_tilt(ctor_params):
    kwargs = {"obj": ctor_params["specimen_tilt"],
              "obj_name": "specimen_tilt"}
    specimen_tilt = czekitout.convert.to_pair_of_floats(**kwargs)
    
    return specimen_tilt



def _pre_serialize_specimen_tilt(specimen_tilt):
    serializable_rep = specimen_tilt

    return serializable_rep



def _de_pre_serialize_specimen_tilt(serializable_rep):
    specimen_tilt = serializable_rep

    return specimen_tilt



def _check_and_convert_scan_config(ctor_params):
    scan_config = prismatique.scan._check_and_convert_scan_config(ctor_params)
    
    return scan_config



def _pre_serialize_scan_config(scan_config):
    if isinstance(scan_config, str):
        serializable_rep = scan_config
    elif isinstance(scan_config, prismatique.scan.rectangular.Params):
        pre_serialize_rectangular_scan_params = \
            prismatique.scan.rectangular._pre_serialize_rectangular_scan_params
        serializable_rep = \
            pre_serialize_rectangular_scan_params(scan_config)
    else:
        serializable_rep = scan_config.tolist()
        for idx, row in enumerate(serializable_rep):
            serializable_rep[idx] = tuple(row)
        serializable_rep = tuple(serializable_rep)

    return serializable_rep



def _de_pre_serialize_scan_config(serializable_rep):
    if isinstance(serializable_rep, str):
        scan_config = serializable_rep
    elif isinstance(serializable_rep, dict):
        mod_alias = prismatique.scan.rectangular
        de_pre_serialize_rectangular_scan_params = \
            mod_alias._de_pre_serialize_rectangular_scan_params
        scan_config = de_pre_serialize_rectangular_scan_params(serializable_rep)
    else:
        kwargs = {"obj": serializable_rep,
                  "obj_name": "serializable_rep"}
        scan_config = \
            czekitout.convert.to_real_two_column_numpy_matrix(**kwargs)

    return scan_config



class ModelParams(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related to the modelling of STEM systems.

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
        slices and frozen phonon configuration subsets respectively.

        Otherwise, if ``sample_specification`` is an instance of the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs`, the class
        :class:`prismatique.sample.SMatrixSubsetIDs`, or the class
        :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`, then
        ``sample_specification`` specifies a set of files, where each file
        stores either the pre-calculated potential slices or :math:`S`-matrices
        for a frozen phonon configuration subset. See the documentation for the
        aforementioned classes for further discussions on specifying
        pre-calculated objects. See the documentation for the subpackage
        :mod:`prismatique.stem` for a discussion on :math:`S`-matrices. Note
        that any :math:`S`-matrices specified in ``sample_specification`` must
        be pre-calculated using the same convergence semiangle and mean beam
        energy as that specified by the model parameters of the probe
        ``probe_model_params`` below.

        Moreover, note that ``sample_specification`` must be an instance of the
        class :class:`prismatique.sample.ModelParams`, or the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs` if using the
        multislice algorithm for the STEM simulation, or if the potential slices
        to be used in the simulation are to be saved as output, per the output
        parameters.
    probe_model_params : :class:`embeam.stem.probe.ModelParams` | `None`, optional
        The model parameters of the probe. See the documentation for the class
        :class:`embeam.stem.probe.ModelParams` for a discussion on said
        parameters. If ``probe_model_params`` is set to `None` [i.e. the default
        value], then the aforementioned model parameters are set to default
        values.
    specimen_tilt : `array_like` (`float`, shape=(``2``,)), optional
        ``specimen_tilt[0]`` specifies the specimen tilt along the
        :math:`x`-axis in mrads; ``specimen_tilt[1]`` specifies the specimen
        tilt along the :math:`y`-axis in mrads. Note that according to
        Ref. [Cowley1]_, the method used to simulate specimen tilt in
        ``prismatic`` is only valid for small tilts up to about 1 degree.
    scan_config : `array_like` (`float`, shape=(``num_positions``, ``2``)) | :class:`prismatique.scan.rectangular.Params` | `str` | `None`, optional
        If ``scan_config`` is a real-valued two-column matrix, then it specifies
        a set of probe positions, where ``scan_config[i][0]`` and
        ``scan_config[i][1]`` specify respectively the :math:`x`- and
        :math:`y`-coordinates of the probe position indexed by ``i``, in units
        of angstroms, with ``0<=i<num_positions`` and ``num_positions`` being
        the number of probe positions. If ``scan_config`` is of the type
        :class:`prismatique.scan.rectangular.Params`, then it specifies a
        rectangular grid-like pattern of probe positions. See the documentation
        for this class for more details. If ``scan_config`` is a string, then
        ``scan_config`` is a path to a file that specifies a set of probe
        positions. The file must be encoded as ASCII text (UTF-8).  The file
        should be formatted as follows: the first line can be whatever header
        the user desires, i.e. the first line is treated as a comment; each
        subsequent line except for the last is of the form "x y", where "x" and
        "y" are the :math:`x`- and :math:`y`-coordinates of a probe position, in
        units of angstroms; and the last line in the file should be "-1".

        Since periodic boundaries conditions (PBCs) are assumed in the
        :math:`x`- and :math:`y`-dimensions, :math:`x`- and
        :math:`y`-coordinates can take on any real numbers. If ``scan_config``
        is set to `None`, [i.e.  the default value], then the probe is to be
        scanned across the entire unit cell of the simulation, in steps of 0.25
        angstroms in both the :math:`x`- and :math:`y`-directions.

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
        {"sample_specification": _check_and_convert_sample_specification,
         "probe_model_params": _check_and_convert_probe_model_params,
         "specimen_tilt": _check_and_convert_specimen_tilt,
         "scan_config": _check_and_convert_scan_config}

    _pre_serialization_funcs = \
        {"sample_specification": _pre_serialize_sample_specification,
         "probe_model_params": _pre_serialize_probe_model_params,
         "specimen_tilt": _pre_serialize_specimen_tilt,
         "scan_config": _pre_serialize_scan_config}

    _de_pre_serialization_funcs = \
        {"sample_specification": _de_pre_serialize_sample_specification,
         "probe_model_params": _de_pre_serialize_probe_model_params,
         "specimen_tilt": _de_pre_serialize_specimen_tilt,
         "scan_config": _de_pre_serialize_scan_config}
    
    def __init__(self,
                 sample_specification,
                 probe_model_params=None,
                 specimen_tilt=(0, 0),
                 scan_config=None):
        ctor_params = {"sample_specification": sample_specification,
                       "probe_model_params": probe_model_params,
                       "specimen_tilt": specimen_tilt,
                       "scan_config": scan_config}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_stem_system_model_params(ctor_params):
    stem_system_model_params = \
        copy.deepcopy(ctor_params["stem_system_model_params"])    
    kwargs = {"obj": stem_system_model_params,
              "obj_name": "stem_system_model_params",
              "accepted_types": (ModelParams,)}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return stem_system_model_params



def _pre_serialize_stem_system_model_params(stem_system_model_params):
    serializable_rep = stem_system_model_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_stem_system_model_params(serializable_rep):
    stem_system_model_params = ModelParams.de_pre_serialize(serializable_rep)

    return stem_system_model_params
    


###########################
## Define error messages ##
###########################