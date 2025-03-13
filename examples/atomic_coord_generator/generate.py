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
"""An example of atomic coordinate generation.

See the link
https://mrfitzpa.github.io/prismatique/examples/atomic-coord-generator/generate.html
for a description of the example.

A NOTE BEFORE STARTING
----------------------

To run this script from the terminal as is, i.e. without modifications, change
into the directory containing said script, and then issue the following
command::

  python generate.py

The atomic coordinates are saved to the file
``<root>/examples/data/atomic-coords.xyz``, where ``<root>`` is the root of the
``prismatique`` repository.

If you would like to modify this script for your own work, it is recommended
that you copy the original script and save it elsewhere outside of the git
repository so that the changes made are not tracked by git.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For making directories and creating path objects.
import pathlib

# For getting the path to the current script.
import os



# For general array handling.
import numpy as np



###############################################
## Define classes, functions, and contstants ##
###############################################



###########################
## Define error messages ##
###########################



#########################
## Main body of script ##
#########################

msg = "Generating the atomic coordinates for the bilayer MoS2 sample..."
print(msg)
print()



# "a" and "c" lattice parameters of MoS2 in Å.
a_MoS2 = 3.1604
c_MoS2 = 12.295

# The atomic potential extent in Å. See the documentation for the core attribute
# ``atomic_potential_extent`` of the class
# :class:`prismatique.sample.ModelParams` for a description of this quantity.
atomic_potential_extent = 3

# MoS2 unit-cell lattice vectors. These vectors will yield a hexagonal lattice.
a_MoS2_1 = a_MoS2 * np.array([1.0, 0.0, 0.0])
a_MoS2_2 = a_MoS2 * np.array([0.0, np.sqrt(3), 0.0])
a_MoS2_3 = c_MoS2 * np.array([0.0, 0.0, 1.0])

# Used to define the positions of atoms in unit cell.
u = 0.612

# Positions of S atoms in unit cell.
delta_S_1 = (1/2)*a_MoS2_1 + (1/6)*a_MoS2_2 + (u-1/2)*a_MoS2_3
delta_S_2 = (1/2)*a_MoS2_1 + (1/6)*a_MoS2_2 + (1-u)*a_MoS2_3
delta_S_3 = (2/3)*a_MoS2_2 + (-1/2+u)*a_MoS2_3
delta_S_4 = (2/3)*a_MoS2_2 + (1-u)*a_MoS2_3
delta_S_5 = (1/3)*a_MoS2_2 + u*a_MoS2_3
delta_S_6 = (1/3)*a_MoS2_2 + (3/2-u)*a_MoS2_3
delta_S_7 = (1/2)*a_MoS2_1 + (5/6)*a_MoS2_2 + u*a_MoS2_3
delta_S_8 = (1/2)*a_MoS2_1 + (5/6)*a_MoS2_2 + (3/2-u)*a_MoS2_3

# Positions of Mo atoms in unit cell.
delta_Mo_1 = (1/3)*a_MoS2_2 + (1/4)*a_MoS2_3
delta_Mo_2 = (1/2)*a_MoS2_1 + (5/6)*a_MoS2_2 + (1/4)*a_MoS2_3
delta_Mo_3 = (1/2)*a_MoS2_1 + (1/6)*a_MoS2_2 + (3/4)*a_MoS2_3
delta_Mo_4 = (2/3)*a_MoS2_2 + (3/4)*a_MoS2_3

# Define the orthorhombic MoS2 unit cell in two parts.
S_unit_cell = np.array([delta_S_1, delta_S_2, delta_S_3, delta_S_4,
                        delta_S_5, delta_S_6, delta_S_7, delta_S_8])
Mo_unit_cell = np.array([delta_Mo_1, delta_Mo_2, delta_Mo_3, delta_Mo_4])

# Specify the number of atomic layers.
num_atomic_layers = 2



# We tile the Mo and S unit cells.
num_y_tiles = 4
num_x_tiles = int(np.round(num_y_tiles * np.sqrt(3)))
x_tiling_indices = range(0, num_x_tiles)
y_tiling_indices = range(0, num_y_tiles)
z_tiling_indices = range(0, (num_atomic_layers//2)+(num_atomic_layers%2))

Mo_sample_unit_cell = []
S_sample_unit_cell = []

for x_tiling_idx in x_tiling_indices:
    for y_tiling_idx in y_tiling_indices:
        for z_tiling_idx in z_tiling_indices:
            shift = (x_tiling_idx*a_MoS2_1
                     + y_tiling_idx*a_MoS2_2
                     + z_tiling_idx*a_MoS2_3)

            if ((z_tiling_idx == max(z_tiling_indices))
                and (num_atomic_layers%2 == 1)):
                current_Mo_cell = np.array([delta_Mo+shift
                                            for delta_Mo
                                            in Mo_unit_cell[:2]])
                current_S_cell = np.array([delta_S+shift
                                           for delta_S
                                           in S_unit_cell[:4]])
            else:
                current_Mo_cell = np.array([delta_Mo+shift
                                            for delta_Mo
                                            in Mo_unit_cell])
                current_S_cell = np.array([delta_S+shift
                                           for delta_S
                                           in S_unit_cell])

            for position_of_current_atom in current_Mo_cell:
                x, y, z = position_of_current_atom
                Mo_sample_unit_cell.append((x, y, z))
                
            for position_of_current_atom in current_S_cell:
                x, y, z = position_of_current_atom
                S_sample_unit_cell.append((x, y, z))

Mo_sample_unit_cell = np.array(Mo_sample_unit_cell)
S_sample_unit_cell = np.array(S_sample_unit_cell)



# Find minimum and maximum z-coordinates of the sample.
single_species_sample_unit_cells = (Mo_sample_unit_cell, S_sample_unit_cell)

min_z = np.inf
max_z = -np.inf

sample_unit_cells = single_species_sample_unit_cells
for single_species_sample_unit_cell in sample_unit_cells:
    candidate_min_z = np.amin(single_species_sample_unit_cell[:, 2])
    candidate_max_z = np.amax(single_species_sample_unit_cell[:, 2])
    min_z = min(candidate_min_z, min_z)
    max_z = max(candidate_max_z, max_z)


                
# Determine sample unit cell dimensions, keeping in mind we want a buffer of 3 Å
# [i.e. the value of ``atomic_potential_extent``] from the edges of the sample
# unit cell along the z-axis.
Delta_X = num_x_tiles * np.linalg.norm(a_MoS2_1)
Delta_Y = num_y_tiles * np.linalg.norm(a_MoS2_2)
Delta_Z = (max_z-min_z) + 2*atomic_potential_extent
sample_unit_cell_dims = (Delta_X, Delta_Y, Delta_Z)



# Apply a global shift to the atomic coordinates along the z-axis.
global_shift = np.array([0, 0, min_z-atomic_potential_extent])

for idx, _ in enumerate(Mo_sample_unit_cell):
    Mo_sample_unit_cell[idx] -= global_shift
for idx, _ in enumerate(S_sample_unit_cell):
    S_sample_unit_cell[idx] -= global_shift


            
# Write atomic coordinates of sample unit cell to file.
path_to_current_script = pathlib.Path(os.path.realpath(__file__))
output_dirname = str(path_to_current_script.parent.parent) + "/data"
output_filename = output_dirname + "/atomic_coords.xyz"

pathlib.Path(output_dirname).mkdir(parents=True, exist_ok=True)

# See the documentation for the core attribute ``atomic_coords_filename`` of the
# class :class:`prismatique.sample.ModelParams` for a description of the
# required format of the file in which to save the atomic coordinates.
with open(output_filename, "w") as file_obj:
    unformatted_line = "Bilayer MoS2 Sample; atomic_potential_extent={}\n"
    formatted_line = unformatted_line.format(atomic_potential_extent)
    file_obj.write(formatted_line)

    unformatted_line = "\t{:18.14f}\t{:18.14f}\t{:18.14f}\n"
    formatted_line = unformatted_line.format(*sample_unit_cell_dims)
    file_obj.write(formatted_line)

    # Number not used except to write file with expected format.
    occ = 1
    
    Z = 42  # Atomic number of Mo.

    # The RMS x-displacement of Mo atoms at room temperature. Value was taken
    # from experimental data for the RMS of the in-plane displacement in
    # Schönfeld et al., Acta Cryst. B39, 404-407 (1983).
    u_x_rms = 0.069

    for position_of_atom in Mo_sample_unit_cell:
        x, y, z = position_of_atom
        unformatted_line = ("{}\t{:18.14f}\t{:18.14f}"
                            "\t{:18.14f}\t{:18.14f}\t{:18.14f}\n")
        formatted_line = unformatted_line.format(Z, x, y, z, occ, u_x_rms)
        file_obj.write(formatted_line)

    Z = 16  # Atomic number of S.

    # The RMS x-displacement of S atoms at room temperature. Value was taken
    # from experimental data for the RMS of the in-plane displacement in
    # Schönfeld et al., Acta Cryst. B39, 404-407 (1983).
    u_x_rms = 0.062

    for position_of_atom in S_sample_unit_cell:
        x, y, z = position_of_atom
        unformatted_line = ("{}\t{:18.14f}\t{:18.14f}"
                            "\t{:18.14f}\t{:18.14f}\t{:18.14f}\n")
        formatted_line = unformatted_line.format(Z, x, y, z, occ, u_x_rms)
        file_obj.write(formatted_line)

    file_obj.write("-1")



msg = ("Finished generating the atomic coordinates for the bilayer MoS2 "
       "sample: the atomic coordinates were written to the file {}"
       ".".format(output_filename))
print(msg)
print()
