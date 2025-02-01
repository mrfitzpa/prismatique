r"""For specifying simulation parameters related to the discretization of 
real-space and Fourier/:math:`k`-space.
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



def _check_and_convert_z_supersampling(ctor_params):
    kwargs = {"obj": ctor_params["z_supersampling"],
              "obj_name": "z_supersampling"}
    z_supersampling = czekitout.convert.to_nonnegative_int(**kwargs)
    
    return z_supersampling



def _pre_serialize_z_supersampling(z_supersampling):
    serializable_rep = z_supersampling

    return serializable_rep



def _de_pre_serialize_z_supersampling(serializable_rep):
    z_supersampling = serializable_rep

    return z_supersampling



def _check_and_convert_sample_supercell_reduced_xy_dims_in_pixels(ctor_params):
    kwargs = {"obj": ctor_params["sample_supercell_reduced_xy_dims_in_pixels"],
              "obj_name": "sample_supercell_reduced_xy_dims_in_pixels"}
    sample_supercell_reduced_xy_dims_in_pixels = \
        czekitout.convert.to_pair_of_positive_ints(**kwargs)
    
    return sample_supercell_reduced_xy_dims_in_pixels



def _pre_serialize_sample_supercell_reduced_xy_dims_in_pixels(
        sample_supercell_reduced_xy_dims_in_pixels):
    serializable_rep = sample_supercell_reduced_xy_dims_in_pixels

    return serializable_rep



def _de_pre_serialize_sample_supercell_reduced_xy_dims_in_pixels(
        serializable_rep):
    sample_supercell_reduced_xy_dims_in_pixels = serializable_rep

    return sample_supercell_reduced_xy_dims_in_pixels



def _check_and_convert_interpolation_factors(ctor_params):
    kwargs = {"obj": ctor_params["interpolation_factors"],
              "obj_name": "interpolation_factors"}
    interpolation_factors = czekitout.convert.to_pair_of_positive_ints(**kwargs)
    
    return interpolation_factors



def _pre_serialize_interpolation_factors(interpolation_factors):
    serializable_rep = interpolation_factors

    return serializable_rep



def _de_pre_serialize_interpolation_factors(serializable_rep):
    interpolation_factors = serializable_rep

    return interpolation_factors



def _check_and_convert_num_slices(ctor_params):
    kwargs = {"obj": ctor_params["num_slices"],
              "obj_name": "num_slices"}
    num_slices = czekitout.convert.to_positive_int(**kwargs)
    
    return num_slices



def _pre_serialize_num_slices(num_slices):
    serializable_rep = num_slices

    return serializable_rep



def _de_pre_serialize_num_slices(serializable_rep):
    num_slices = serializable_rep

    return num_slices




class Params(fancytypes.PreSerializableAndUpdatable):
    r"""The simulation parameters related to the discretization of real-space 
    and Fourier/:math:`k`-space.

    In ``prismatic``, a sample is constructed by defining an orthorhombic unit
    cell of atoms, which may be optionally tiled a finite number of times in the
    :math:`x`-, :math:`y`-, and :math:`z`-directions to form a supercell, upon
    which periodic boundary conditions are imposed on the supercell in the
    :math:`x`- and :math:`y`-directions. In the case of crystal samples that do
    not possess primitive cells with orthorhombic symmetry, one can typically
    find a larger unit cell within the crystal with the required orthorhombic
    symmetry. 

    Both the PRISM and the standard multislice algorithms make certain
    assumptions regarding the total potential within the sample. First, the
    total potential is treated classically. Second, the isolated atom
    approximation is assumed, wherein the total potential within
    :math:`V\left(\mathbf{r}\right)` the sample is modelled by:

    .. math ::
        V\left(\mathbf{r}\right)=\sum_{j=1}^{N}V_{j}^{\left(\text{atom}\right)}
        \left(\left|\mathbf{r}-\mathbf{r}_{j}\right|\right),
        :label: isolated_atom_approx

    where

    .. math ::
        \mathbf{r}=x\hat{\mathbf{x}}+y\hat{\mathbf{y}}+z\hat{\mathbf{z}},
        :label: r_vector_in_potential_params

    .. math ::
        \mathbf{r}_{j}=x_{j}\hat{\mathbf{x}}+y_{j}\hat{\mathbf{y}}
        +z_{j}\hat{\mathbf{z}},
        :label: r_j_vector

    :math:`N` is the number of atoms within the sample,
    :math:`V_{j}^{\left(\text{atom}\right)}
    \left(\left|\mathbf{r}-\mathbf{r}_{j}\right|\right)` is the effective
    isolated potential of the :math:`j^{\text{th}}` atom at the position
    :math:`\mathbf{r}_{j}` [we will refer to :math:`j` as the atom label]. In
    the isolated atom approximation bond effects are neglected. Corrections due
    to the bond effects should be small in most cases since electron scattering
    is mainly from the nucleus with the core and valence electrons screening the
    nucleus [Kirkland1]_. The errors due to the isolated atom approximation are
    most prominent for low angle scattering, which may cause problems with phase
    contrast bright field images, whereas high angle annular dark field STEM
    images should be more accurate [Kirkland1]_. In any case, due to the
    simplicity and cheap computational costs, ``prismatic`` adopts the isolated
    atom approximation described above. The effective isolated potentials are
    parameterized based on single-atom scattering factors determined through
    Hartree-Fock calculations [Kirkland1]_, [DaCosta1]_.  Third, the beam
    electrons are assumed to interact with an effective potential
    :math:`\tilde{V}\left(\mathbf{r}\right)` within the sample, which is
    obtained by coarse-graining :math:`V\left(\mathbf{r}\right)` along the
    :math:`z`-axis:

    .. math ::
        \tilde{V}\left(\mathbf{r}\right)=\sum_{n=0}^{N_{\text{slices}}-1}
        \delta\left(z-z_{n}\right)
        \tilde{V}_{n}\left(x,y\right),
        :label: multislice_approx

    where :math:`N_{\text{slices}}` is a positive integer;

    .. math ::
        z_{n}=n\delta z,
        :label: z_n_in_potential_params

    with

    .. math ::
        \delta z = \frac{\Delta Z}{N_{\text{slices}}},
        :label: slice_thickness_in_potential_params

    :math:`\Delta Z` is the :math:`z`-dimension of the sample's supercell in
    units of length; and

    .. math ::
        \tilde{V}_{n}\left(x,y\right) & =\int_{z_{n}}^{z_{n+1}}dz\,
        V\left(\mathbf{r}\right)\\
        & =\sum_{j=1}^{N}\int_{z_{n}}^{z_{n+1}}dz\,
        V_{j}^{\left(\text{atom}\right)}
        \left(\left|\mathbf{r}-\mathbf{r}_{j}\right|\right).
        :label: coarse_grained_potential_1

    In other words, the sample is partitioned into :math:`N_{\text{slices}}`
    slices of thickness :math:`\delta z` along the :math:`z`-axis and
    :math:`V\left(\mathbf{r}\right)` is coarse-grained over these slices to
    obtain :math:`\tilde{V}\left(\mathbf{r}\right)`. Note that :math:`z=z_0` and
    :math:`z=z_{N_{\text{slices}}}` are the :math:`z`-coordinates of the
    entrance and exit surfaces of the sample's supercell respectively.  We will
    refer to :math:`\tilde{V}_{n}\left(x,y\right)` as the :math:`n^{\text{th}}`
    potential slice.

    At this point, the user can choose one of two approaches to complete the
    calculation of the effective potential
    :math:`\tilde{V}_{n}\left(x,y\right)`.  In the simpler approach, we
    approximate Eq. :eq:`coarse_grained_potential_1` by

    .. math ::
        \tilde{V}_{n}\left(x,y\right)=\sum_{j\in A_{n}}
        \tilde{V}_{j}^{\left(\text{atom}\right)}\left(x-x_{j},y-y_{j}\right),
        :label: coarse_grained_potential_2

    where the :math:`\tilde{V}_{j}^{\left(\text{atom}\right)}
    \left(x-x_{j},y-y_{j}\right)` are the 2D projected isolated atomic
    potentials:

    .. math ::
        \tilde{V}_{j}^{\left(\text{atom}\right)}\left(x-x_{j},y-y_{j}\right)
        =\int_{-\infty}^{\infty}dz\,V_{j}^{\left(\text{atom}\right)}
        \left(\left|\mathbf{r}-\mathbf{r}_{j}\right|\right),
        :label: projected_potentials

    and :math:`A_{n}` is the set of atom labels corresponding to the atoms with
    :math:`z`-coordinates :math:`z_{j}\in\left[z_{n},z_{n+1}\right)`. The
    integrals in Eq. :eq:`projected_potentials` can be evaluated analytically,
    thus avoiding costly numerical integrations. This simpler approach should
    only yield accurate results when the slice thickness :math:`\delta z` is
    large compared to the effective extent of each isolated atomic potential in
    the sample and that the :math:`z_{j}\sim\delta z\left(1+2l\right)/2`, where
    :math:`0\le l<N_{\text{slices}}`. In the second approach, the limits of
    integration in Eq. :eq:`coarse_grained_potential_1` are not extended to
    infinity, and the integrals are evaluated numerically. This approach is more
    accurate but is more computationally expensive.

    Both real-space and momentum/Fourier/angular-space need to be discretized in
    order to handle wavefunctions, as well as the potential slices
    numerically. Each potential slice is discretized with pixel sizes in the
    :math:`x`- and :math:`y`-directions given by:

    .. math ::
        \Delta x=\frac{\Delta X}{N_{x}},\quad\Delta y=\frac{\Delta Y}{N_{y}},
        :label: potential_slice_pixel_size

    where :math:`\Delta x` and :math:`\Delta y` are the potential slice pixel
    sizes along the :math:`x`- and :math:`y`-directions respectively,
    :math:`\Delta X` and :math:`\Delta Y` are the :math:`x`- and
    :math:`y`-dimensions of the sample's supercell in units of length
    respectively, and :math:`N_x` and :math:`N_y` are the :math:`x`- and
    :math:`y`-dimensions of the sample's supercell in units of potential slice
    pixels respectively. We will refer to :math:`\left(N_{x},N_{y}\right)` as
    the "pixelated potential slice :math:`xy`-dimensions". For practical
    reasons, it is convenient to refer to the quantity:

    .. math ::
        \left(\tilde{N}_{x},\tilde{N}_{y}\right)=\left(\frac{N_{x}}{4f_{x}},
        \frac{N_{y}}{4f_{y}}\right),
        :label: reduced_supercell_xy_dims_in_pixels_in_potential_params

    as the "sample supercell's reduced :math:`xy`-dimensions in units of
    pixels", where :math:`f_{x}` and :math:`f_{y}` are positive integers called
    the :math:`x` and :math:`y` interpolation factors, which are discussed
    further in the documentation for the subpackage :mod:`prismatique.stem` and
    the class :class:`prismatique.tilt.Params`. Note that in STEM simulations
    that use the multislice algorithm, :math:`f_{x}=f_{y}=1`. In order for
    :math:`\left(\tilde{N}_{x},\tilde{N}_{y}\right)` to be a pair of integers,
    we constrain :math:`N_{x}` and :math:`N_{y}` to be divisible by
    :math:`4f_{x}` and :math:`4f_{y}` respectively. In ``prismatique``, the user
    specifies indirectly :math:`\Delta x` and :math:`\Delta y` by specifying
    directly :math:`\left(\tilde{N}_{x},\tilde{N}_{y}\right)`.

    Parameters
    ----------
    z_supersampling : `int`, optional
        See the above discussion for context. If ``z_supersampling`` is set to
        ``0``, then the first of the two approaches discussed above is used to
        calculate :math:`\tilde{V}_{n}\left(x,y\right)`,
        i.e. Eq. :eq:`coarse_grained_potential_2` is used. Otherwise, if
        ``z_supersampling`` is set to a positive integer, then
        Eq. :eq:`coarse_grained_potential_1` is used to calculate 
        :math:`\tilde{V}_{n}\left(x,y\right)`, and ``z_supersampling`` specifies
        the number of evenly spaced quadrature points to use in numerically
        evaluating the integral. Note that ``z_supersampling`` must be a
        non-negative integer.
    sample_supercell_reduced_xy_dims_in_pixels : `array_like` (`int`, shape=(``2``,)), optional
        ``sample_supercell_reduced_xy_dims_in_pixels[0]`` and
        ``sample_supercell_reduced_xy_dims_in_pixels[1]`` specify
        :math:`\tilde{N}_x` and :math:`\tilde{N}_y` respectively, where
        :math:`\left(\tilde{N}_{x},\tilde{N}_{y}\right)` are the sample 
        supercell's reduced :math:`xy`-dimensions in units of pixels.
    interpolation_factors : `array_like` (`int`, shape=(``2``,)), optional
        ``interpolation_factors[0]`` and ``interpolation_factors[1]`` are the
        interpolation factors :math:`f_x` and :math:`f_y` respectively, which
        are discussed above. Note that ``interpolation_factors`` must be set to
        ``(1, 1)`` if using the multislice algorithm for STEM simulations.
    num_slices : `int`, optional
        The number of slices :math:`N_{\text{slices}}` used to partition the 
        sample.

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
        {"z_supersampling": _check_and_convert_z_supersampling,
         "sample_supercell_reduced_xy_dims_in_pixels": \
         _check_and_convert_sample_supercell_reduced_xy_dims_in_pixels,
         "interpolation_factors": _check_and_convert_interpolation_factors,
         "num_slices": _check_and_convert_num_slices}

    _pre_serialization_funcs = \
        {"z_supersampling": _pre_serialize_z_supersampling,
         "sample_supercell_reduced_xy_dims_in_pixels": \
         _pre_serialize_sample_supercell_reduced_xy_dims_in_pixels,
         "interpolation_factors": _pre_serialize_interpolation_factors,
         "num_slices": _pre_serialize_num_slices}

    _de_pre_serialization_funcs = \
        {"z_supersampling": _de_pre_serialize_z_supersampling,
         "sample_supercell_reduced_xy_dims_in_pixels": \
         _de_pre_serialize_sample_supercell_reduced_xy_dims_in_pixels,
         "interpolation_factors": _de_pre_serialize_interpolation_factors,
         "num_slices": _de_pre_serialize_num_slices}
    
    def __init__(self,
                 z_supersampling=16,
                 sample_supercell_reduced_xy_dims_in_pixels=(64, 64),
                 interpolation_factors=(1, 1),
                 num_slices=25):
        ctor_params = {"z_supersampling": z_supersampling,
                       "sample_supercell_reduced_xy_dims_in_pixels": \
                       sample_supercell_reduced_xy_dims_in_pixels,
                       "interpolation_factors": interpolation_factors,
                       "num_slices": num_slices}
        fancytypes.PreSerializableAndUpdatable.__init__(self, ctor_params)

        return None



def _check_and_convert_discretization_params(ctor_params):
    discretization_params = copy.deepcopy(ctor_params["discretization_params"])
    if discretization_params is None:
        discretization_params = Params()
    
    kwargs = {"obj": discretization_params,
              "obj_name": "discretization_params",
              "accepted_types": (Params, type(None))}
    czekitout.check.if_instance_of_any_accepted_types(**kwargs)

    return discretization_params



def _pre_serialize_discretization_params(discretization_params):
    serializable_rep = discretization_params.pre_serialize()

    return serializable_rep



def _de_pre_serialize_discretization_params(serializable_rep):
    discretization_params = Params.de_pre_serialize(serializable_rep)

    return discretization_params
    


###########################
## Define error messages ##
###########################