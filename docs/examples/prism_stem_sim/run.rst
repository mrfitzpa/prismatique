.. _prism_stem_sim_example_sec:

Example of running a STEM simulation using the PRISM algorithm
==============================================================

Example description
-------------------

In this example, we simulate a STEM experiment involving the bilayer
:math:`\text{MoS}_2` sample that we defined in the
:ref:`atomic_coord_generator_example_sec` page. The STEM simulation is
implemented using the PRISM algorithm. See the documentation for the subpackage
:mod:`prismatique.stem` for a discussion on the PRISM algorithm.

Code
----

Below is the code that performs the aforementioned STEM simulation. You can also
find the same code in the file ``<root>/examples/prism_stem_sim/run.py`` of the
repository, where ``<root>`` is the root of the ``prismatique`` repository. To
run the script from the terminal, change into the directory containing said
script, and then issue the following command::

  python run.py

The output files generated by this script are saved in the directory
``<root>/examples/data/prism_stem_sim_output``, where ``<root>`` is the root of
the ``prismatique`` repository. To analyze the output, use the Jupyter notebook
``<root>/examples/output_data_analysis_notebooks/analyzing_prism_stem_sim_output.ipynb``.

If you would like to modify this script for your own work, it is recommended
that you copy the original script and save it elsewhere outside of the git
repository so that the changes made are not tracked by git.

.. literalinclude:: ../../../examples/prism_stem_sim/run.py