.. _atomic_coord_generator_example_sec:

Example of atomic coordinate generation
=======================================



Prerequisites for running this example
--------------------------------------

Prior to running any scripts or Jupyter notebooks in the directory
``<root>/examples``, where ``<root>`` is the root of the ``prismatique``
repository, a set of Python libraries need to be installed in the Python
environment within which any such scripts or Jupyter notebooks are to be
executed. See :ref:`this page <examples_prerequisites_for_running_examples_sec>`
for instructions on how to do so.



Example description
-------------------

The first step to running a STEM or HRTEM simulation is to generate the
zero-temperature expectation value of the atomic coordinates, as well as the
root-mean-square :math:`x`-displacement, for each atom in the unit cell of your
sample of interest. This first step does not actually involve the use of any
functions from either the ``prismatic`` or ``prismatique``
libraries. Nonetheless, for completeness we show an example of atomic coordinate
generation below.

In this example, we generate the zero-temperature expectation value of the
atomic coordinates for each atom in the unit cell of a sample of bilayer
:math:`\text{MoS}_2`, where the :math:`c`-axis of the :math:`\text{MoS}_2`
lattice is parallel to the optic axis. See the documentation for the class
:class:`prismatique.discretization.Params` for a discussion on sample unit cells
and supercells. We make a distinction here between *sample* unit cells and
*lattice* unit cells: the latter is a smaller unit cell which can be tiled to
construct the former. Lattice unit cells are convenient when the sample unit
cell possesses translational symmetry over a significant portion of the volume
it takes up, e.g. a crystalline nanoparticle.

It is easiest to work with a sample unit cell with square dimensions in the
:math:`xy`-plane when possible. Strictly speaking, in the case of the
:math:`\text{MoS}_2` sample mentioned above, it is impossible to construct a
square supercell that is also subject to periodic boundary conditions. That
being said, one can construct a sample supercell with nearly square dimensions
in the :math:`xy`-plane by tiling appropriately the chosen orthorhombic
lattice unit cell. We choose the following real-space lattice vectors to define
the orthorhombic lattice unit cell of our :math:`\text{MoS}_2` sample:


.. math ::
    \begin{align*}
    \mathbf{a}_{\text{MoS}_2;1}&=a_{\text{MoS}_2}\hat{\mathbf{x}},\\
    \mathbf{a}_{\text{MoS}_2;2}&=\sqrt{3}a_{\text{MoS}_2}\hat{\mathbf{y}},\\
    \mathbf{a}_{\text{MoS}_2;3}&=c_{\text{MoS}_2}\hat{\mathbf{z}},
    \end{align*}
    :label: MoS2_rspace_lattice_vectors

where :math:`a_{\text{MoS}_2}=3.1604\,\text{Å}` and
:math:`c_{\text{MoS}_2}=12.295\,\text{Å}`. Since
:math:`\left|\mathbf{a}_{\text{MoS}_2;2}\right| /
\left|\mathbf{a}_{\text{MoS}_2;1}\right|` equals an irrational number, one can
never construct a supercell with square dimensions in the
:math:`xy`-dimensions. If we tile the lattice unit cell 7 times in the
:math:`x`-direction, and 4 times in the :math:`y`-direction, then we obtain a
supercell with nearly square dimensions in the :math:`xy`-plane. The resulting
:math:`x`- and :math:`y`-dimensions, :math:`\Delta X` and :math:`\Delta Y`, of
the supercell are:

.. math ::
    \Delta X=7 a_{\text{MoS}_2},\quad\Delta Y=6.93 a_{\text{MoS}_2},
    :label: MoS2_supercell_dims

In this example, we use the aforementioned tiling to construct our supercell.

The positions of the S atoms in the lattice unit cell are:

.. math ::
    \begin{align*}
    \delta_{\text{S};1}&=\frac{1}{2}\mathbf{a}_{\text{MoS}_{2};1}
    +\frac{1}{6}\mathbf{a}_{\text{MoS}_{2};2}
    +\left(u-\frac{1}{2}\right)\mathbf{a}_{\text{MoS}_{2};3},\\
    \delta_{\text{S};2}&=\frac{1}{2}\mathbf{a}_{\text{MoS}_{2};1}
    +\frac{1}{6}\mathbf{a}_{\text{MoS}_{2};2}
    +\left(1-u\right)\mathbf{a}_{\text{MoS}_{2};3},\\
    \delta_{\text{S};3}&=\frac{2}{3}\mathbf{a}_{\text{MoS}_{2};2}
    +\left(-\frac{1}{2}+u\right)\mathbf{a}_{\text{MoS}_{2};3},\\
    \delta_{\text{S};4}&=\frac{2}{3}\mathbf{a}_{\text{MoS}_{2};2}
    +\left(1-u\right)\mathbf{a}_{\text{MoS}_{2};3},\\
    \delta_{\text{S};5}&=\frac{1}{3}\mathbf{a}_{\text{MoS}_{2};2}
    +u\mathbf{a}_{\text{MoS}_{2};3},\\
    \delta_{\text{S};6}&=\frac{1}{3}\mathbf{a}_{\text{MoS}_{2};2}
    +\left(\frac{3}{2}-u\right)\mathbf{a}_{\text{MoS}_{2};3},\\
    \delta_{\text{S};7}&=\frac{1}{2}\mathbf{a}_{\text{MoS}_{2};1}
    +\frac{5}{6}\mathbf{a}_{\text{MoS}_{2};2}+u\mathbf{a}_{\text{MoS}_{2};3},\\
    \delta_{\text{S};8}&=\frac{1}{2}\mathbf{a}_{\text{MoS}_{2};1}
    +\frac{5}{6}\mathbf{a}_{\text{MoS}_{2};2}
    +\left(\frac{3}{2}-u\right)\mathbf{a}_{\text{MoS}_{2};3},
    \end{align*}
    :label: S_atoms_in_lattice_unit_cell

where :math:`u=0.612`. The positions of the Mo atoms in the lattice unit cell
are:

.. math ::
    \begin{align*}
    \delta_{\text{Mo};1}&=\frac{1}{3}\mathbf{a}_{\text{MoS}_{2};2}
    +\frac{1}{4}\mathbf{a}_{\text{MoS}_{2};3},\\
    \delta_{\text{Mo};2}&=\frac{1}{2}\mathbf{a}_{\text{MoS}_{2};1}
    +\frac{5}{6}\mathbf{a}_{\text{MoS}_{2};2}
    +\frac{1}{4}\mathbf{a}_{\text{MoS}_{2};3},\\
    \delta_{\text{Mo};3}&=\frac{1}{2}\mathbf{a}_{\text{MoS}_{2};1}
    +\frac{1}{6}\mathbf{a}_{\text{MoS}_{2};2}
    +\frac{3}{4}\mathbf{a}_{\text{MoS}_{2};3},\\
    \delta_{\text{Mo};4}&=\frac{1}{2}\mathbf{a}_{\text{MoS}_{2};1}
    +\frac{5}{6}\mathbf{a}_{\text{MoS}_{2};2}
    +\frac{3}{4}\mathbf{a}_{\text{MoS}_{2};3}.
    \end{align*}
    :label: Mo_atoms_in_lattice_unit_cell

Note that a single lattice unit cell tiled laterally forms two atomic layers of
:math:`\text{MoS}_2`.


      
Code
----

Below is the code that generates the atomic coordinates described above. You can
also find the same code in the file
``<root>/examples/atomic_coord_generator/generate.py`` of the repository, where
``<root>`` is the root of the ``prismatique`` repository. To run the script from
the terminal, change into the directory containing said script, and then issue
the following command::

  python generate.py

The atomic coordinates are saved to the file
``<root>/examples/data/atomic_coords.xyz``.

If you would like to modify this script for your own work, it is recommended
that you copy the original script and save it elsewhere outside of the git
repository so that the changes made are not tracked by git.

.. literalinclude:: ../../../examples/atomic_coord_generator/generate.py
