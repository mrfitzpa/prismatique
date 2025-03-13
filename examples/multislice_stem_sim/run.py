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
"""An example of running a STEM simulation that uses the multislice algorithm.

See the link
https://mrfitzpa.gitlab.io/prismatique/examples/multislice-stem-sim/run.html for
a description of the example.

A NOTE BEFORE STARTING
----------------------

To run this script from the terminal as is, i.e. without modifications, change
into the directory containing said script, and then issue the following
command::

  python run.py

The output files generated by this script are saved in the directory
``<root>/examples/data/multislice-stem-sim-output``, where ``<root>`` is the
root of the ``prismatique`` repository. To analyze the output, use the Jupyter
notebook
``<root>/examples/output-data-analysis-notebooks/analyzing-multislice-stem-sim-output.ipynb``.

If you would like to modify this script for your own work, it is recommended
that you copy the original script and save it elsewhere outside of the git
repository so that the changes made are not tracked by git.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For making directories and creating path objects.
import pathlib

# For spawning shell processes and getting the path to current script.
import os



# For running the STEM simulation.
import prismatique



###############################################
## Define classes, functions, and contstants ##
###############################################



###########################
## Define error messages ##
###########################



#########################
## Main body of script ##
#########################

msg = "Setting up the STEM simulation for the MoS2 sample..."
print(msg)
print()



# In another example script, we generate all the simulation parameters required
# to generate the potential slices, the S-matrices, and for simulating STEM and
# HRTEM experiments involving the MoS2 sample. Rather than set explicitly again
# all the simulation parameters required for the STEM simulation, we can run the
# other script to generate the parameters if they do not already exist, and then
# load these parameters in this current script.
path_to_current_script = pathlib.Path(os.path.realpath(__file__))
path_to_examples_dir = str(path_to_current_script.parent.parent)
path_to_data_dir = path_to_examples_dir + "/data"
path_to_serialized_params = (path_to_data_dir
                             + "/sim_param_generator_output"
                             + "/multislice_stem_sim_params.json")
for try_count in range(2):
    try:
        multislice_stem_sim_params = \
            prismatique.stem.sim.Params.load(path_to_serialized_params)
        stem_system_model_params = \
            multislice_stem_sim_params.core_attrs["stem_system_model_params"]
        sample_specification = \
            stem_system_model_params.core_attrs["sample_specification"]
        atomic_coords_filename = \
            sample_specification.core_attrs["atomic_coords_filename"]

        module_alias = prismatique.sample
        module_alias.check_atomic_coords_file_format(atomic_coords_filename)

        break  # No exceptions raised, hence we can exit the for-loop early.

    except:
        msg = ("Need to generate the simulation parameter sets for the various "
               "simulated experiments involving the MoS2 sample...")
        print(msg)
        print()
        script_to_execute = (path_to_examples_dir
                             + "/sim_param_generator/generate.py")
        os.system("python " + script_to_execute)
        msg = ("Resuming the setup of the STEM simulation for the MoS2 "
               "sample...")
        print(msg)
        print()

sample_model_params = sample_specification



# In another example script, we generate the potential slices for the MoS2
# sample, which can be used in turn as the sample model to perform the STEM
# simulation, rather than the sample model parameters we loaded above. Below we
# check whether the potential slices have been generated: if they exists then we
# load and use them in this current script to perform the STEM simulation,
# otherwise we use the sample model parameters that we loaded above to perform
# the STEM simulation from scratch.
path_to_potential_slices = (path_to_data_dir
                            + "/potential_slice_generator_output")
try:
    thermal_params = \
        sample_model_params.core_attrs["thermal_params"]
    num_frozen_phonon_configs_per_subset = \
        thermal_params.core_attrs["num_frozen_phonon_configs_per_subset"]
    num_subsets = \
        thermal_params.core_attrs["num_subsets"]

    filenames = tuple((path_to_potential_slices
                       + "/potential_slices_of_subset_{}.h5".format(idx))
                      for idx in range(num_subsets))

    # See the documentation for the class
    # :class:`prismatique.sample.PotentialSliceSubsetIDs` for a description of
    # its construction parameters and how they relate to specifying
    # precalculated potential slices.
    discretization_params = \
        sample_model_params.core_attrs["discretization_params"]
    interpolation_factors = \
        discretization_params.core_attrs["interpolation_factors"]
    kwargs = \
        {"filenames": filenames,
         "interpolation_factors": interpolation_factors,
         "max_num_frozen_phonon_configs_per_subset": \
         num_frozen_phonon_configs_per_subset}
    sample_specification = \
        prismatique.sample.PotentialSliceSubsetIDs(**kwargs)


        
    # Lastly, we use the potential slices and the parameters loaded above to
    # perform the STEM simulation. See the documentation for the function
    # :func:`prismatique.stem.sim.run` for a discussion on the input parameters
    # of said function.
    new_core_attr_subset_candidate = {"sample_specification": \
                                      sample_specification}
    stem_system_model_params.update(new_core_attr_subset_candidate)
    new_core_attr_subset_candidate = {"stem_system_model_params": \
                                      stem_system_model_params}
    multislice_stem_sim_params.update(new_core_attr_subset_candidate)

    msg = "Running the STEM simulation for the MoS2 sample..."
    print(msg)
    print()    
    prismatique.stem.sim.run(sim_params=multislice_stem_sim_params)
        
except:
    # Could not find valid potential slices to load, hence we perform the STEM
    # simulation from scratch using the sample model parameters.
    new_core_attr_subset_candidate = {"sample_specification": \
                                      sample_model_params}
    stem_system_model_params.update(new_core_attr_subset_candidate)
    new_core_attr_subset_candidate = {"stem_system_model_params": \
                                      stem_system_model_params}
    multislice_stem_sim_params.update(new_core_attr_subset_candidate)
    prismatique.stem.sim.run(sim_params=multislice_stem_sim_params)

    

output_params = multislice_stem_sim_params.core_attrs["output_params"]
base_output_params = output_params.core_attrs["base_params"]
output_dirname = base_output_params.core_attrs["output_dirname"]
msg = ("Finished the STEM simulation for the MoS2 sample: all output files "
       "were saved in the directory {}.".format(output_dirname))
print(msg)
print()
