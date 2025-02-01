r"""For specifying the output parameters for STEM simulations.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For performing deep copies of objects.
import copy



# For general array handling and special math functions.
import numpy as np

# For validating and converting objects.
import czekitout.check
import czekitout.convert

# For defining classes that support enforced validation, updatability,
# pre-serialization, and de-serialization.
import fancytypes

# For checking whether certain numbers are approximately zero.
import emconstants

# For validating, pre-serializing, and de-pre-serializing instances of the
# class :class:`embeam.stem.probe.ModelParams`.
import embeam.stem.probe

# For calculating the radial range from zero to the largest radial distance
# within the signal space boundaries.
import empix



# Import child modules and packages of current package.
import prismatique.stem.output.base
import prismatique.stem.output.multislice
import prismatique.stem.output.prism

# For validating instances of the class
# :class:`prismatique.stem.system.ModelParams`.
import prismatique.stem.system

# For validating instances of the classes
# :class:`prismatique.sample.ModelParams`,
# :class:`prismatique.sample.PotentialSliceSubsetIDs`,
# :class:`prismatique.sample.SMatrixSubsetIDs`, and
# :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`; and for
# calculating quantities related to the modelling of the sample.
import prismatique.sample

# For generating probe positions.
import prismatique.scan

# For postprocessing CBED intensity patterns.
import prismatique._signal



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
__all__ = ["Params",
           "layer_depths",
           "data_size"]



def _check_and_convert_base_params(ctor_params):
    check_and_convert_base_params = \
        prismatique.stem.output.base._check_and_convert_base_params
    base_params = \
        check_and_convert_base_params(ctor_params)
    
    return base_params



def _pre_serialize_base_params(base_params):
    pre_serialize_base_params = \
        prismatique.stem.output.base._pre_serialize_base_params
    serializable_rep = \
        pre_serialize_base_params(base_params)

    return serializable_rep



def _de_pre_serialize_base_params(serializable_rep):
    de_pre_serialize_base_params = \
        prismatique.stem.output.base._de_pre_serialize_base_params
    base_params = \
        de_pre_serialize_base_params(serializable_rep)

    return base_params



def _check_and_convert_alg_specific_params(ctor_params):
    temp_ctor_params = \
        {"multislice_output_params": ctor_params["alg_specific_params"]}
    mod_alias = \
        prismatique.stem.output
    
    try:
        check_and_convert_multislice_output_params = \
            mod_alias.multislice._check_and_convert_multislice_output_params
        alg_specific_params = \
            check_and_convert_multislice_output_params(temp_ctor_params)
    except:
        try:
            temp_ctor_params = \
                {"prism_output_params": ctor_params["alg_specific_params"]}
            check_and_convert_prism_output_params = \
                mod_alias.prism._check_and_convert_prism_output_params
            alg_specific_params = \
                check_and_convert_prism_output_params(temp_ctor_params)
        except:
            raise TypeError(_check_and_convert_alg_specific_params_err_msg_1)
    
    return alg_specific_params



def _pre_serialize_alg_specific_params(alg_specific_params):
    mod_alias = prismatique.stem.output

    if "z_start_output" in alg_specific_params.core_attrs:
        pre_serialize_multislice_output_params = \
            mod_alias.multislice._pre_serialize_multislice_output_params
        serializable_rep = \
            pre_serialize_multislice_output_params(alg_specific_params)
    else:
        pre_serialize_prism_output_params = \
            mod_alias.prism._pre_serialize_prism_output_params
        serializable_rep = \
            pre_serialize_prism_output_params(alg_specific_params)

    return serializable_rep



def _de_pre_serialize_alg_specific_params(serializable_rep):
    mod_alias = prismatique.stem.output

    if "num_slices_per_output" in serializable_rep:
        de_pre_serialize_multislice_output_params = \
            mod_alias.multislice._de_pre_serialize_multislice_output_params
        alg_specific_params = \
            de_pre_serialize_multislice_output_params(serializable_rep)
    else:
        de_pre_serialize_prism_output_params = \
            mod_alias.prism._de_pre_serialize_prism_output_params
        alg_specific_params = \
            de_pre_serialize_prism_output_params(serializable_rep)

    return alg_specific_params



class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The output parameters for STEM simulations.

    See the documentation for the class :class:`prismatique.cbed.Params` for a
    discussion on CBED patterns that is relevant to the discussion on this page.
    We refer to CBED patterns generated from a STEM simulation as 4D-STEM data.

    ``prismatic`` simulates annular detectors shaped in concentric rings, where
    each ring has a fixed bin width that can be chosen by the user. This allows
    one to collect the integrated intensity over each ring in a given
    DP. Integrating over the rings this way for each CBED intensity pattern
    yields what we will refer to as the "3D_STEM intensity data/output": there
    are still two spatial dimensions associated with each probe position, but
    now there is only one dimension associated with the diffraction data for
    each probe position, i.e.  a radial dimension in momentum space.

    We can collect a further reduced set of data by integrating/summing over a
    consecutive set of rings for each CBED intensity pattern. For example, for
    each CBED intensity pattern we could integrate over the first so many rings
    that encompasses the direct beam, thus collecting bright field
    data. Therefore, for each CBED intensity pattern we extract a single number,
    which is the integrated intensity over a solid ring [e.g. in dark field
    imaging] or a disk [e.g. in bright field imaging]. We will call the
    resulting data the "2D_STEM intensity data/output".

    From the 4D-STEM intensity data, we can also calculate the "center-of-mass"
    (COM) momentum for each probe position:

    .. math ::
        \left\langle \mathbf{k}_{\text{COM}}\right\rangle =
        \frac{\int d\mathbf{k}\,
        I_{\text{STEM}}\left(k_{x},k_{y}\right)\mathbf{k}}{\int d\mathbf{k}\,
        I_{\text{STEM}}\left(k_{x},k_{y}\right)},
        :label: COM

    where :math:`I_{\text{STEM}}\left(k_{x},k_{y}\right)` is the CBED intensity.

    Upon the completion of a STEM simulation, ``prismatique`` can optionally
    save a variety of different STEM data to a set of HDF5 files based on the
    specifications of the user. This is done by taking the original output file
    generated by ``prismatic``, and then restructuring it into one or more
    output files. If the output parameters specify that intensity data be saved,
    then said data is extracted from the original ``prismatic`` output files,
    postprocessed, and written to a new file with basename
    ``"stem_sim_intensity_output.h5"`` in a more readable layout. Assuming all
    intensity data is to be saved, this new output file has the following
    structure:

    - metadata: <HDF5 group>

      * probe_positions: <HDF5 2D dataset>
        
        + dim 1: "probe idx"
        + dim 2: "vector component idx [0->x, 1->y]"
        + units: "Å"
        + pattern type: "rectangular grid" | "jittered rectangular grid" | "no underlying rectangular grid"
        + grid dimensions in units of probe shifts: <ordered pair> | "N/A"

      * output_layer_depths: <HDF5 1D dataset>
        
        + dim 1: "output layer idx"
        + units: "Å"

      * k_x: <HDF5 1D dataset>
        
        + dim 1: "k_x idx"
        + units: "1/Å"

      * k_y: <HDF5 1D dataset>
        
        + dim 1: "k_y idx"
        + units: "1/Å"

      * k_xy: <HDF5 1D dataset>
        
        + dim 1: "k_xy idx"
        + units: "1/Å"

    - data: <HDF5 group>

      - 2D_STEM: <HDF5 group>

        * integrated_intensities: <HDF5 2D dataset>

          + dim 1: "output layer idx"
          + dim 2: "probe idx"
          + units: "dimensionless"
          + lower integration limit in mrads: <scalar>
          + upper integration limit in mrads: <scalar>

      - 3D_STEM: <HDF5 group>

        * integrated_intensities: <HDF5 3D dataset>

          + dim 1: "output layer idx"
          + dim 2: "probe idx"
          + dim 3: "k_xy idx"
          + units: "dimensionless"

      - 4D_STEM: <HDF5 group>

        * intensity_DPs: <HDF5 4D dataset>

          + dim 1: "output layer idx"
          + dim 2: "probe idx"
          + dim 3: "k_y idx"
          + dim 4: "k_x idx"
          + units: "dimensionless"

      * center_of_mass_momentum: <HDF5 3D dataset>
        
          + dim 1: "output layer idx"
          + dim 2: "vector component idx [0->x, 1->y]"
          + dim 3: "probe idx"
          + units: "1/Å"

    Note that the sub-bullet points listed immediately below a given HDF5
    dataset display the HDF5 attributes associated with said HDF5 dataset. Each
    HDF5 scalar and dataset has a ``"units"`` attribute which, as the name
    suggests, indicates the units in which said data [i.e. the scalar or
    dataset] is expressed. Each HDF5 dataset will also have a set of attributes
    with names of the form ``"dim {}".format(i)`` with ``i`` being an integer
    ranging from 1 to the rank of said HDF5 dataset. Attribute ``"dim
    {}".format(i)`` of a given HDF5 dataset labels the :math:`i^{\text{th}}`
    dimension of the underlying array of the dataset. Most of these dimension
    labels should be self-explanatory but for clarification: "idx" is short for
    "index"; "k_x" and "k_y" refer to the :math:`x`- and :math:`y`-components of
    the electron scattering momentum; and the notion of an output layer is discussed
    in the documentation for the class
    :class:`prismatique.stem.output.multislice.Params`.

    If the output parameters specify that complex-valued wavefunction data be
    saved, then said data is extracted from the original ``prismatic`` output
    files, and written to a new set of files: one file per frozen phonon/atomic
    configuration subset. Note that, unlike the intensity data, the
    complex-valued wavefunction data is not postprocessed. See the documentation
    for the class :class:`prismatique.thermal.Params` for a discussion on frozen
    phonon configurations and their grouping into subsets. For the ``i`` th
    frozen phonon configuration subset, the corresponding wavefunction data is
    saved to the file with basename
    ``"stem_sim_wavefunction_output_of_subset_"+str(i)+".h5"``. Each one of
    these new output files has the following structure:

    - metadata: <HDF5 group>

      * probe_positions: <HDF5 2D dataset>
        
        + dim 1: "probe idx"
        + dim 2: "vector component idx [0->x, 1->y]"
        + units: "Å"
        + pattern type: "rectangular grid" | "jittered rectangular grid" | "no underlying rectangular grid"
        + grid dimensions in units of probe shifts: <ordered pair> | "N/A"

      * defocii: <HDF5 1D dataset>
        
        + dim 1: "defocus idx"
        + units: "Å"

      * output_layer_depths: <HDF5 1D dataset>
        
        + dim 1: "output layer idx"
        + units: "Å"

      * k_x: <HDF5 1D dataset>
        
        + dim 1: "k_x idx"
        + units: "1/Å"

      * k_y: <HDF5 1D dataset>
        
        + dim 1: "k_y idx"
        + units: "1/Å"

    - data: <HDF5 group>

      - 4D_STEM: <HDF5 group>

        * complex_valued_DPs: <HDF5 6D dataset>

          + dim 1: "output layer idx"
          + dim 2: "atomic config idx"
          + dim 3: "defocus idx"
          + dim 4: "probe idx"
          + dim 5: "k_y idx"
          + dim 6: "k_x idx"
          + units: "dimensionless"

    ``prismatique`` can also optionally save the "potential slice" [i.e.
    Eq. :eq:`coarse_grained_potential_1`] data for each subset of frozen phonon
    configurations into separate HDF5 output files by extracting said data from
    the original ``prismatic`` output files. For the ``i`` th subset, the
    corresponding potential slice data is saved to an HDF5 output file with
    basename ``"potential_slices_"+str(i)+".h5"``. Unlike the output data in the
    file ``"stem_simulation_output.h5"``, the layout of the potential slice data
    is kept the same as that found in the original file [i.e. the same HDF5
    paths are used].  The same layout needs to be used in order for
    ``prismatic`` to be able to successfully import/load pre-calculated
    potential slice data for future simulations.

    Moreover if used/generated, ``prismatique`` can also optionally save the
    :math:`S`-matrix [i.e.  Eq. :eq:`S_tilde_matrix`] data into a separate HDF5
    output files by extracting said data from the original ``prismatic`` output
    files. For the ``i`` th subset, the corresponding :math:`S`-matrix data is
    saved to an HDF5 output file with basename ``"S_matrix_"+str(i)+".h5"``. The
    layout of the :math:`S`-matrix data is kept the same as that found in the
    original file for the same reason mentioned for the potential slice
    data. Note that :math:`S`-matrix data is only generated/used when STEM
    simulations are performed using the PRISM algorithm [see the documentation
    for the subpackage :class:`prismatique.stem` for details on this algorithm].

    It is beyond the scope of the documentation to describe the structure of the
    potential slice and :math:`S`-matrix output files. Users can analyze the
    data in these output files with the help of the tools found in the module
    :mod:`prismatique.load`.

    The last file that is always generated after running a STEM simulation is a
    JSON file that contains, in a serialized format, the simulation parameters
    used. This file has the basename ``"stem_sim_params.json"``.

    Parameters
    ----------
    base_params : :class:`prismatique.stem.output.base.Params` | `None`, optional
        The base output parameters for the STEM simulation. See the
        documentation for the class :class:`prismatique.stem.output.base.Params`
        for a discussion on said parameters. If ``base_params`` is set to `None`
        [i.e. the default value], then the aforementioned parameters are set to
        default values.
    alg_specific_params : :class:`prismatique.stem.output.multislice.Params` | :class:`prismatique.stem.output.prism.Params` | `None`, optional
        The output parameters that are applicable only to the particular
        algorithm used perform the STEM simulation. If ``alg_specific_params``
        is set to an instance of the class
        :class:`prismatique.stem.output.multislice.Params`, then the multislice
        algorithm is used to perform the STEM simulation. See the documentation
        for the class :class:`prismatique.stem.output.multislice.Params` for a
        discussion on the output parameters specific to the multislice
        algorithm. If ``alg_specific_params`` is set to an instance of the class
        :class:`prismatique.stem.output.prism.Params`, then the PRISM algorithm
        is used to perform the STEM simulation. See the documentation for the
        class :class:`prismatique.stem.output.prism.Params` for a discussion on
        the output parameters specific to the PRISM algorithm. If
        ``alg_specific_params`` is set to `None` [i.e. the default value], then
        the multislice algorithm is used to perform the STEM simulation with the
        aforementioned output parameters set to default values.

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
        {"base_params": _check_and_convert_base_params,
         "alg_specific_params": _check_and_convert_alg_specific_params}

    _pre_serialization_funcs = \
        {"base_params": _pre_serialize_base_params,
         "alg_specific_params": _pre_serialize_alg_specific_params}

    _de_pre_serialization_funcs = \
        {"base_params": _de_pre_serialize_base_params,
         "alg_specific_params": _de_pre_serialize_alg_specific_params}
    
    def __init__(self,
                 base_params=None,
                 alg_specific_params=None):
        ctor_params = {"base_params": base_params,
                       "alg_specific_params": alg_specific_params}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_output_params(ctor_params):
    output_params = copy.deepcopy(ctor_params["output_params"])
    if output_params is None:
        output_params = Params()
        
    kwargs = {"obj": output_params,
              "obj_name": "output_params",
              "accepted_types": (Params, type(None))}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return output_params



def _pre_serialize_output_params(output_params):
    serializable_rep = output_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_output_params(serializable_rep):
    output_params = Params.de_pre_serialize(serializable_rep)

    return output_params



def layer_depths(sample_specification, alg_specific_params=None):
    r"""Determine the output layer depths from STEM simulation parameter subset.

    For a discussion on output layer depths, see the documentation for the class
    :class:`prismatique.stem.output.multislice.Params`, particularly the
    description of the constructor parameters ``num_slices_per_output`` and
    ``z_start_output``.

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
        slices and frozen phonon configuration subsets respectively. Note that
        of parameters stored in ``sample_specification``, only the following are
        used:

        - sample_specification

          * atomic_coords_filename

          * unit_cell_tiling

          * discretization_params

            + interpolation_factors

            + num_slices

        Otherwise, if ``sample_specification`` is an instance of the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs`, the class
        :class:`prismatique.sample.SMatrixSubsetIDs`, or the class
        :class:`prismatique.sample.PotentialSliceAndSMatrixSubsetIDs`, then
        ``sample_specification`` specifies a set of files, where each file
        stores either the pre-calculated potential slices or :math:`S`-matrices
        for a frozen phonon configuration subset. See the documentation for the
        aforementioned classes for further discussions on specifying
        pre-calculated objects. See the documentation for the subpackage
        :mod:`prismatique.stem` for a discussion on :math:`S`-matrices. 

        Moreover, note that ``sample_specification`` must be an instance of the
        class :class:`prismatique.sample.ModelParams`, or the class
        :class:`prismatique.sample.PotentialSliceSubsetIDs` when the parameter
        ``alg_specific_params`` below is set to `None`, or is an instance of the
        class :class:`prismatique.stem.output.multislice.Params`.
    alg_specific_params : :class:`prismatique.stem.output.multislice.Params` | :class:`prismatique.stem.output.prism.Params` | `None`, optional
        The output parameters that are applicable only to the particular
        algorithm used perform the STEM simulation. If ``alg_specific_params``
        is set to an instance of the class
        :class:`prismatique.stem.output.multislice.Params`, then the multislice
        algorithm is used to perform the STEM simulation. See the documentation
        for the class :class:`prismatique.stem.output.multislice.Params` for a
        discussion on the output parameters specific to the multislice
        algorithm. If ``alg_specific_params`` is set to an instance of the class
        :class:`prismatique.stem.output.prism.Params`, then the PRISM algorithm
        is used to perform the STEM simulation. See the documentation for the
        class :class:`prismatique.stem.output.prism.Params` for a discussion on
        the output parameters specific to the PRISM algorithm. If
        ``alg_specific_params`` is set to `None` [i.e. the default value], then
        the multislice algorithm is used to perform the STEM simulation with the
        aforementioned output parameters set to default values.

    Returns
    -------
    output_layer_depths : `array_like` (`float`, ndim=1)
        The output layer depths, in units of angstroms, in ascending order. Note
        that for STEM simulations that use the PRISM algorithm,
        ``output_layer_depths`` is a single-element array, where the single
        element is the the :math:`z`-dimension of the sample's supercell in
        units of angstroms [see the documentation for the class
        :class:`prismatique.discretization.Params` for a discussion on the
        sample's supercell].

    """
    accepted_types = (prismatique.sample.ModelParams,
                      prismatique.sample.PotentialSliceSubsetIDs,
                      prismatique.sample.SMatrixSubsetIDs,
                      prismatique.sample.PotentialSliceAndSMatrixSubsetIDs)
    prismatique.sample._check_sample_specification(sample_specification,
                                                   accepted_types)

    temp_ctor_params = \
        {"alg_specific_params": alg_specific_params}
    alg_specific_params = \
        _check_and_convert_alg_specific_params(temp_ctor_params)

    _check_sample_specification_wrt_alg_specific_params(sample_specification,
                                                        alg_specific_params)

    output_layer_depths = _layer_depths(sample_specification,
                                        alg_specific_params)
    
    return output_layer_depths



def _check_sample_specification_wrt_alg_specific_params(sample_specification,
                                                        alg_specific_params):
    if "z_start_output" in alg_specific_params.core_attrs:
        if "thermal_params" in sample_specification.core_attrs:
            discretization_params = \
                sample_specification.core_attrs["discretization_params"]
            f_x, f_y = \
                discretization_params.core_attrs["interpolation_factors"]
        elif "interpolation_factors" in sample_specification.core_attrs:
            f_x, f_y = sample_specification.core_attrs["interpolation_factors"]
        else:
            err_msg = \
                _check_sample_specification_wrt_alg_specific_params_err_msg_1
            raise TypeError(err_msg)

        if ((f_x != 1) or (f_y != 1)):
            err_msg = \
                _check_sample_specification_wrt_alg_specific_params_err_msg_2
            raise ValueError(err_msg)
    
    return None



def _layer_depths(sample_specification, alg_specific_params):
    multislice_alg_was_specified = \
        ("z_start_output" in alg_specific_params.core_attrs)
    
    if multislice_alg_was_specified:
        output_layer_depths = \
            _output_layer_depths_for_multislice_alg(sample_specification,
                                                    alg_specific_params)
    else:
        output_layer_depths = \
            (prismatique.sample._supercell_dims(sample_specification)[2],)

    return output_layer_depths



def _output_layer_depths_for_multislice_alg(sample_specification,
                                            alg_specific_params):
    _, _, Delta_Z = \
        prismatique.sample._supercell_dims(sample_specification)
    sample_supercell_slice_thickness = \
        prismatique.sample._supercell_slice_thickness(sample_specification)

    total_num_slices = _total_num_slices(sample_specification)

    func_alias = prismatique.sample._z_start_output_passed_to_pyprismatic
    z_start_output = func_alias(sample_specification, alg_specific_params)
    
    z_slice_idx_start_output = \
        max(np.ceil(z_start_output / sample_supercell_slice_thickness)-1, 0)
    num_slices_per_output = \
        alg_specific_params.core_attrs["num_slices_per_output"]

    if num_slices_per_output == 0:
        output_layer_depths = (Delta_Z,)
        return output_layer_depths
        
    total_num_output_layers = \
        int(total_num_slices / num_slices_per_output)
    total_num_output_layers += \
        int(total_num_slices % num_slices_per_output != 0)
    if z_start_output > 0:
        total_num_output_layers += \
            int((z_slice_idx_start_output+1) % num_slices_per_output == 0)
        total_num_output_layers -= \
            int((z_slice_idx_start_output+1) / num_slices_per_output)

        output_layer_idx_offset = \
            int((z_slice_idx_start_output+1) / num_slices_per_output)
        output_layer_idx_offset += \
            int((z_slice_idx_start_output+1) % num_slices_per_output != 0)
    else:
        output_layer_idx_offset = 1

    distance_between_output_layers = \
        num_slices_per_output * sample_supercell_slice_thickness
    depth_of_first_output_layer = \
        output_layer_idx_offset * distance_between_output_layers

    output_layer_depths = [depth_of_first_output_layer]*total_num_output_layers
    for output_layer_idx in range(1, total_num_output_layers-1):
        output_layer_depths[output_layer_idx] += \
            output_layer_idx * distance_between_output_layers
    output_layer_depths[-1] = Delta_Z
    output_layer_depths = tuple(output_layer_depths)

    return output_layer_depths



def _total_num_slices(sample_specification):
    if "discretization_params" in sample_specification.core_attrs:
        discretization_params = \
            sample_specification.core_attrs["discretization_params"]
        total_num_slices = \
            discretization_params.core_attrs["num_slices"]
    else:
        _, _, Delta_Z = \
            prismatique.sample._supercell_dims(sample_specification)
        sample_supercell_slice_thickness = \
            prismatique.sample._supercell_slice_thickness(sample_specification)
        total_num_slices = \
            int(np.round(Delta_Z / sample_supercell_slice_thickness))

    return total_num_slices



def data_size(stem_system_model_params, output_params=None):
    r"""Calculate the data size of the STEM simulation output that one could
    generate according to a given STEM system model, and output parameter set.

    Note that data size due to HDF5 file overhead and metadata are not taken
    into account.

    Parameters
    ----------
    stem_system_model_params : :class:`prismatique.stem.system.ModelParams`
        The simulation parameters related to the modelling of the STEM
        system. See the documentation for the class
        :class:`prismatique.stem.system.ModelParams` for a discussion on said
        parameters.
    output_params : :class:`prismatique.stem.output.Params` | `None`, optional
        The output parameters for the STEM simulation. See the documentation for
        the class :class:`prismatique.stem.output.Params` for a discussion on
        said parameters. If ``output_params`` is set to `None` [i.e. the default
        value], then the aforementioned simulation parameters are set to default
        values.

    Returns
    -------
    output_data_size : `int`
        The data size in units of bytes.

    """
    stem_system_model_params = \
        _check_stem_system_model_params(stem_system_model_params)
    temp_ctor_params = \
        {"output_params": output_params}
    output_params = \
        _check_and_convert_output_params(temp_ctor_params)

    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]

    kwargs = {"sample_specification": sample_specification,
              "output_params": output_params,
              "format_arg": "sample_specification"}
    _check_sample_specification_wrt_output_params(**kwargs)
        
    output_data_size = _data_size(stem_system_model_params, output_params)
        
    return output_data_size



def _check_stem_system_model_params(stem_system_model_params):
    temp_ctor_params = \
        {"stem_system_model_params": stem_system_model_params}
    check_and_convert_stem_system_model_params = \
        prismatique.stem.system._check_and_convert_stem_system_model_params
    stem_system_model_params = \
        check_and_convert_stem_system_model_params(temp_ctor_params)

    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]
    probe_model_params = \
        stem_system_model_params.core_attrs["probe_model_params"]
    scan_config = \
        stem_system_model_params.core_attrs["scan_config"]

    mod_alias = \
        prismatique.sample
    func_alias = \
        mod_alias._check_and_convert_sample_specification_and_probe_model_params
    
    temp_ctor_params = {"sample_specification": sample_specification,
                        "probe_model_params": probe_model_params}
    _, _ = func_alias(temp_ctor_params)

    kwargs = {"sample_specification": sample_specification,
              "scan_config": scan_config,
              "filename": None}
    _ = prismatique.scan._generate_probe_positions(**kwargs)

    return stem_system_model_params



def _check_sample_specification_wrt_output_params(sample_specification,
                                                  output_params,
                                                  format_arg):
    sample_specification_specifies_a_precalculated_S_matrix_subset = \
        (not (("thermal_params" in sample_specification.core_attrs)
              or ("interpolation_factors" in sample_specification.core_attrs)))
    
    base_params = output_params.core_attrs["base_params"]
    save_potential_slices = base_params.core_attrs["save_potential_slices"]
    if (save_potential_slices
        and sample_specification_specifies_a_precalculated_S_matrix_subset):
        unformatted_err_msg = \
            _check_sample_specification_wrt_output_params_err_msg_1
        err_msg = \
            unformatted_err_msg.format(format_arg)
        raise ValueError(err_msg)

    alg_specific_params = output_params.core_attrs["alg_specific_params"]
    if "z_start_output" in alg_specific_params.core_attrs:
        try:
            kwargs = {"sample_specification": sample_specification,
                      "alg_specific_params": alg_specific_params}
            _check_sample_specification_wrt_alg_specific_params(**kwargs)
        except TypeError:
            unformatted_err_msg = \
                _check_sample_specification_wrt_output_params_err_msg_2
            err_msg = \
                unformatted_err_msg.format(format_arg)
            raise TypeError(err_msg)
        except ValueError:
            unformatted_err_msg = \
                _check_sample_specification_wrt_output_params_err_msg_3
            err_msg = \
                unformatted_err_msg.format(format_arg)
            raise ValueError(err_msg)
        except BaseException as err:
            raise err
    
    return None



def _data_size(stem_system_model_params, output_params):    
    base_output_params = output_params.core_attrs["base_params"]
    alg_specific_output_params = output_params.core_attrs["alg_specific_params"]

    output_data_size = 0

    kwargs = {"stem_system_model_params": stem_system_model_params,
              "output_params": output_params}
    if _intensity_output_is_to_be_saved(base_output_params):
        output_data_size += _data_size_of_intensity_output(**kwargs)
    if _wavefunction_output_is_to_be_saved(base_output_params):
        output_data_size += _data_size_of_wavefunction_output(**kwargs)

    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]    
    kwargs = \
        {"sample_specification": sample_specification}

    if base_output_params.core_attrs["save_potential_slices"]:
        output_data_size += \
            prismatique.sample._potential_slice_set_data_size(**kwargs)
    if "save_S_matrices" in alg_specific_output_params.core_attrs:
        if alg_specific_output_params.core_attrs["save_S_matrices"]:
            kwargs["probe_model_params"] = \
                stem_system_model_params.core_attrs["probe_model_params"]
            output_data_size += \
                prismatique.sample._S_matrix_set_data_size(**kwargs)
            
    return output_data_size



def _intensity_output_is_to_be_saved(base_output_params):
    radial_step_size_for_3d_stem = \
        base_output_params.core_attrs["radial_step_size_for_3d_stem"]
    radial_range_for_2d_stem = \
        base_output_params.core_attrs["radial_range_for_2d_stem"]
    cbed_params = \
        base_output_params.core_attrs["cbed_params"]

    if (cbed_params.core_attrs["save_final_intensity"]
        or base_output_params.core_attrs["save_com"]
        or (radial_step_size_for_3d_stem > 0)
        or (radial_range_for_2d_stem[0] != radial_range_for_2d_stem[1])):
        intensity_output_is_to_be_saved = True
    else:
        intensity_output_is_to_be_saved = False

    return intensity_output_is_to_be_saved



def _wavefunction_output_is_to_be_saved(base_output_params):
    cbed_params = \
        base_output_params.core_attrs["cbed_params"]
    wavefunction_output_is_to_be_saved = \
        cbed_params.core_attrs["save_wavefunctions"]

    return wavefunction_output_is_to_be_saved



def _data_size_of_intensity_output(stem_system_model_params, output_params):
    data_size_of_2d_stem_intensity_output = \
        _data_size_of_2d_stem_intensity_output(stem_system_model_params,
                                               output_params)
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]

    base_output_params = \
        output_params.core_attrs["base_params"]
    radial_step_size_for_3d_stem = \
        base_output_params.core_attrs["radial_step_size_for_3d_stem"]
    radial_range_for_2d_stem = \
        base_output_params.core_attrs["radial_range_for_2d_stem"]
    cbed_params = \
        base_output_params.core_attrs["cbed_params"]

    data_size_of_intensity_output = 0
    if radial_range_for_2d_stem[0] != radial_range_for_2d_stem[1]:
        data_size_of_intensity_output += data_size_of_2d_stem_intensity_output
    if base_output_params.core_attrs["save_com"]:
        data_size_of_intensity_output += 2*data_size_of_2d_stem_intensity_output
    if cbed_params.core_attrs["save_final_intensity"]:
        kwargs = \
            {"sample_specification": sample_specification,
             "signal_is_cbed_pattern_set": True,
             "postprocessing_seq": cbed_params.core_attrs["postprocessing_seq"]}
        num_pixels_in_postprocessed_2d_signal_space = \
            prismatique._signal._num_pixels_in_postprocessed_2d_signal_space
        num_pixels_in_postprocessed_dp = \
            num_pixels_in_postprocessed_2d_signal_space(**kwargs)
        data_size_of_intensity_output += \
            (data_size_of_2d_stem_intensity_output
             * num_pixels_in_postprocessed_dp)
    if radial_step_size_for_3d_stem > 0:
        kwargs = \
            {"stem_system_model_params": stem_system_model_params,
             "output_params": output_params}
        num_pixels_in_3d_stem_integrated_dp = \
            _num_pixels_in_3d_stem_integrated_dp(**kwargs)
        data_size_of_intensity_output += \
            (data_size_of_2d_stem_intensity_output
             * num_pixels_in_3d_stem_integrated_dp)

    return data_size_of_intensity_output



def _data_size_of_2d_stem_intensity_output(stem_system_model_params,
                                           output_params):
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]
    scan_config = \
        stem_system_model_params.core_attrs["scan_config"]

    alg_specific_params = output_params.core_attrs["alg_specific_params"]
    output_layer_depths = _layer_depths(sample_specification,
                                        alg_specific_params)

    kwargs = {"sample_specification": sample_specification,
              "scan_config": scan_config,
              "filename": None}
    probe_positions = prismatique.scan._generate_probe_positions(**kwargs)

    size_of_single = 4  # In bytes.

    data_size_of_2d_stem_intensity_output = (len(output_layer_depths)
                                             * len(probe_positions)
                                             * size_of_single)

    return data_size_of_2d_stem_intensity_output



def _num_pixels_in_3d_stem_integrated_dp(stem_system_model_params,
                                         output_params):
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]

    base_output_params = output_params.core_attrs["base_params"]
    cbed_params = base_output_params.core_attrs["cbed_params"]
    postprocessing_seq = cbed_params.core_attrs["postprocessing_seq"]

    kwargs = \
        {"sample_specification": sample_specification,
         "postprocessing_seq": postprocessing_seq,
         "navigation_dims": tuple(),
         "signal_is_cbed_pattern_set": True,
         "signal_dtype": "float"}
    blank_postprocessed_dp_signal = \
        prismatique._signal._blank_postprocessed_2d_signal(**kwargs)

    num_pixels_in_3d_stem_integrated_dp = \
        _num_bins_in_3d_stem_integration(stem_system_model_params,
                                         output_params,
                                         blank_postprocessed_dp_signal)
    
    return num_pixels_in_3d_stem_integrated_dp



def _num_bins_in_3d_stem_integration(stem_system_model_params,
                                     output_params,
                                     input_signal):
    probe_model_params = \
        stem_system_model_params.core_attrs["probe_model_params"]
    gun_model_params = \
        probe_model_params.core_attrs["gun_model_params"]
    mean_beam_energy = \
        gun_model_params.core_attrs["mean_beam_energy"]
    wavelength = \
        embeam.wavelength(mean_beam_energy)
    
    base_output_params = \
        output_params.core_attrs["base_params"]
    radial_angular_step_size = \
        base_output_params.core_attrs["radial_step_size_for_3d_stem"]
    radial_k_step_size = \
        (radial_angular_step_size / 1000) / wavelength

    # The following code block calculates the radial range from zero to the
    # largest radial distance within the signal space boundaries.
    kwargs = \
        {"radial_range": None,
         "signal": input_signal,
         "center": (0, 0)}
    radial_k_range = \
        np.array(empix._check_and_convert_radial_range_v2(**kwargs)) / 1000

    num_bins = int(np.floor(2 * radial_k_range[1] / radial_k_step_size))

    return num_bins



def _data_size_of_wavefunction_output(stem_system_model_params, output_params):
    sample_specification = \
        stem_system_model_params.core_attrs["sample_specification"]
    probe_model_params = \
        stem_system_model_params.core_attrs["probe_model_params"]
    scan_config = \
        stem_system_model_params.core_attrs["scan_config"]

    defocal_offset_supersampling = \
        probe_model_params.core_attrs["defocal_offset_supersampling"]

    alg_specific_params = output_params.core_attrs["alg_specific_params"]
    output_layer_depths = _layer_depths(sample_specification,
                                        alg_specific_params)

    kwargs = \
        {"sample_specification": sample_specification}
    total_num_frozen_phonon_configs = \
        prismatique.sample._total_num_frozen_phonon_configs(**kwargs)

    kwargs = {"sample_specification": sample_specification,
              "scan_config": scan_config,
              "filename": None}
    probe_positions = prismatique.scan._generate_probe_positions(**kwargs)

    kwargs = \
        {"sample_specification": sample_specification,
         "signal_is_cbed_pattern_set": True}
    num_pixels_in_unprocessed_dp = \
        prismatique._signal._num_pixels_in_unprocessed_2d_signal_space(**kwargs)

    size_of_single = 4  # In bytes.
    size_of_complex_single = 2 * size_of_single  # In bytes.

    data_size_of_wavefunction_output = (total_num_frozen_phonon_configs
                                        * defocal_offset_supersampling
                                        * len(output_layer_depths)
                                        * len(probe_positions)
                                        * num_pixels_in_unprocessed_dp
                                        * size_of_complex_single)

    return data_size_of_wavefunction_output



###########################
## Define error messages ##
###########################

_check_and_convert_alg_specific_params_err_msg_1 = \
    ("The object ``alg_specific_params`` must be an instance of one of the "
     "following classes: (`prismatique.stem.output.multislice.Params`, "
     "`prismatique.stem.output.prism.Params`, `NoneType`).")

_check_sample_specification_wrt_alg_specific_params_err_msg_1 = \
    ("The object ``sample_specification`` must be an instance of the class "
     "`prismatique.sample.ModelParams`, or the class "
     "`prismatique.sample.PotentialSliceSubsetIDs` when the parameter "
     "``alg_specific_params`` is set to `None`, or is an instance of the class "
     ":class:`prismatique.stem.output.multislice.Params`.")
_check_sample_specification_wrt_alg_specific_params_err_msg_2 = \
    ("The object "
     "``sample_specification`` must specify that the interpolation factors "
     "be set to unity when the parameter ``alg_specific_params`` is set to "
     "`None`, or is an instance of the class "
     ":class:`prismatique.stem.output.multislice.Params`.")

_check_sample_specification_wrt_output_params_err_msg_1 = \
    ("The object ``{}`` must be an instance of the class "
     "`prismatique.sample.ModelParams`, or the class "
     "`prismatique.sample.PotentialSliceSubsetIDs` when specifying that the "
     "potential slices to be used in the simulation are to be saved as output.")
_check_sample_specification_wrt_output_params_err_msg_2 = \
    ("The object ``{}`` must be an instance of the class "
     "`prismatique.sample.ModelParams`, or the class "
     "`prismatique.sample.PotentialSliceSubsetIDs` when specifying that the "
     "multislice algorithm is to be used, as is implied by the object "
     "``output_params``.")
_check_sample_specification_wrt_output_params_err_msg_3 = \
    ("The object ``{}`` must specify that the interpolation factors be set to "
     "unity when also specifying that the multislice algorithm is to be used, "
     "as is implied by the object ``output_params``.")
