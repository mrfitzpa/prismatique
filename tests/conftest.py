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
r"""Contains fixtures to be used for various tests.

"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For removing directories.
import shutil

# For performing operations on file and directory paths.
import pathlib



# For operations related to unit tests.
import pytest

# For general array handling.
import numpy as np

# For modelling beams and lenses in electron microscopy.
import embeam

# For postprocessing diffraction patterns (DPs) and images.
import empix

# For saving objects to HDF5 files.
import h5pywrappers



# For simulating electron microscopy experiments using Multislice algorithms.
import prismatique



##################################
## Define classes and functions ##
##################################



class Helpers():
    @staticmethod
    def remove_output_files():
        filenames = (Helpers.generate_atomic_coords_filename_1(),
                     Helpers.generate_atomic_coords_filename_2(),
                     Helpers.generate_atomic_coords_filename_3(),
                     Helpers.generate_atomic_coords_filename_4(),
                     Helpers.generate_atomic_coords_filename_5(),
                     Helpers.generate_scan_config_filename_1(),
                     Helpers.generate_scan_config_filename_2(),
                     Helpers.generate_scan_config_filename_3())
        for filename in filenames:
            pathlib.Path(filename).unlink(missing_ok=True)

        output_dirname = Helpers.generate_output_dirname_1()
        shutil.rmtree(output_dirname, ignore_errors=True)

        return None



    @staticmethod
    def generate_atomic_coords_file_1():
        filename = \
            Helpers.generate_atomic_coords_filename_1()

        method_alias = \
            Helpers.generate_C_sample_unit_cell_and_C_sample_unit_cell_dims
        C_sample_unit_cell, C_sample_unit_cell_dims = \
            method_alias(Delta_X_eq_Delta_Y=False)

        with open(filename, "w") as file_obj:
            line = "Graphene Sample 1\n"
            file_obj.write(line)

            Helpers.write_atomic_coord_data_to_file(file_obj,
                                                    C_sample_unit_cell,
                                                    C_sample_unit_cell_dims)

        return None



    @staticmethod
    def generate_atomic_coords_filename_1():
        atomic_coords_filename = "atomic_coords_file_1.xyz"

        return atomic_coords_filename



    @staticmethod
    def generate_C_sample_unit_cell_and_C_sample_unit_cell_dims(
            Delta_X_eq_Delta_Y):
        a_1, a_2, a_3 = \
            Helpers.generate_graphene_unit_cell_lattice_vectors()
        C_unit_cell = \
            Helpers.generate_C_unit_cell()
        x_tiling_indices, y_tiling_indices, z_tiling_indices = \
            Helpers.generate_x_y_and_z_tiling_indices()
        atomic_potential_extent = \
            Helpers.generate_atomic_potential_extent_1()

        C_sample_unit_cell = tuple()

        for x_tiling_idx in x_tiling_indices:
            for y_tiling_idx in y_tiling_indices:
                for z_tiling_idx in z_tiling_indices:
                    shift = (x_tiling_idx*a_1
                             + y_tiling_idx*a_2
                             + z_tiling_idx*a_3)

                    current_C_cell = np.array([delta_C+shift
                                               for delta_C
                                               in C_unit_cell])

                    for position_of_current_atom in current_C_cell:
                        x, y, z_prime = position_of_current_atom
                        C_sample_unit_cell += ((x, y, z_prime),)

        C_sample_unit_cell = np.array(C_sample_unit_cell)

        min_coords = tuple()
        max_coords = tuple()
        for axis in range(3):
            min_coords += (np.amin(C_sample_unit_cell[:, axis]),)
            max_coords += (np.amax(C_sample_unit_cell[:, axis]),)

        Delta_X = (max_coords[0]-min_coords[0])
        Delta_Y = (Delta_X
                   if Delta_X_eq_Delta_Y
                   else (max_coords[1]-min_coords[1]))
        Delta_Z = (max_coords[2]-min_coords[2]) + 2*atomic_potential_extent
        C_sample_unit_cell_dims = (Delta_X, Delta_Y, Delta_Z)

        global_shift = np.array([min_coords[0],
                                 min_coords[1],
                                 min_coords[2]-atomic_potential_extent])
        C_sample_unit_cell[:] -= global_shift

        return C_sample_unit_cell, C_sample_unit_cell_dims



    @staticmethod
    def generate_C_unit_cell():
        a_1, a_2, a_3 = Helpers.generate_graphene_unit_cell_lattice_vectors()

        delta_C_1 = (0/6)*a_1 + (0/6)*a_2 + (0/6)*a_3
        delta_C_2 = (0/6)*a_1 + (2/6)*a_2 + (0/6)*a_3
        delta_C_3 = (3/6)*a_1 + (3/6)*a_2 + (0/6)*a_3
        delta_C_4 = (3/6)*a_1 + (5/6)*a_2 + (0/6)*a_3

        C_unit_cell = np.array([delta_C_1, delta_C_2, delta_C_3, delta_C_4])

        return C_unit_cell



    @staticmethod
    def generate_graphene_unit_cell_lattice_vectors():
        a = 2.46

        a_1 = a * np.array([1.0, 0.0, 0.0])
        a_2 = a * np.array([0.0, np.sqrt(3), 0.0])
        a_3 = np.array([0.0, 0.0, 1.0])

        graphene_unit_cell_lattice_vectors = (a_1, a_2, a_3)

        return graphene_unit_cell_lattice_vectors



    @staticmethod
    def generate_x_y_and_z_tiling_indices():
        num_y_tiles = 3
        num_x_tiles = int(np.round(num_y_tiles * np.sqrt(3)))
        
        x_tiling_indices = \
            range(-(num_x_tiles//2), -(num_x_tiles//2)+num_x_tiles)
        y_tiling_indices = \
            range(-(num_y_tiles//2), -(num_y_tiles//2)+num_y_tiles)
        z_tiling_indices = \
            range(0, 1)

        return x_tiling_indices, y_tiling_indices, z_tiling_indices



    @staticmethod
    def generate_atomic_potential_extent_1():
        atomic_potential_extent = 3

        return atomic_potential_extent



    @staticmethod
    def write_atomic_coord_data_to_file(file_obj,
                                        C_sample_unit_cell,
                                        C_sample_unit_cell_dims):
        unformatted_line = "\t{:18.14f}\t{:18.14f}\t{:18.14f}\n"
        formatted_line = unformatted_line.format(*C_sample_unit_cell_dims)
        file_obj.write(formatted_line)

        occ = 1
        Z_of_C = 6
        u_x_rms_of_C = 0.104

        single_species_sample_unit_cells = (C_sample_unit_cell,)
        Z_set = (Z_of_C,)
        u_x_rms_set = (u_x_rms_of_C,)
        zip_obj = zip(single_species_sample_unit_cells, Z_set, u_x_rms_set)
        
        for triplet in zip_obj:
            single_species_sample_unit_cell, Z, u_x_rms = triplet
            for position_of_atom in single_species_sample_unit_cell:
                x, y, z_prime = position_of_atom
                unformatted_line = ("{}\t{:18.14f}\t{:18.14f}"
                                    "\t{:18.14f}\t{:18.14f}\t{:18.14f}\n")
                args = (Z, x, y, z_prime, occ, u_x_rms)
                formatted_line = unformatted_line.format(*args)
                file_obj.write(formatted_line)

        file_obj.write("-1")

        return None



    @staticmethod
    def generate_atomic_coords_file_2():
        filename = \
            Helpers.generate_atomic_coords_filename_2()

        method_alias = \
            Helpers.generate_C_sample_unit_cell_and_C_sample_unit_cell_dims
        C_sample_unit_cell, C_sample_unit_cell_dims = \
            method_alias(Delta_X_eq_Delta_Y=True)

        with open(filename, "w") as file_obj:
            line = "Graphene Sample 2\n"
            file_obj.write(line)

            Helpers.write_atomic_coord_data_to_file(file_obj,
                                                    C_sample_unit_cell,
                                                    C_sample_unit_cell_dims)

        return None



    @staticmethod
    def generate_atomic_coords_filename_2():
        atomic_coords_filename = "atomic_coords_file_2.xyz"

        return atomic_coords_filename



    @staticmethod
    def generate_scan_config_file_1():
        filename = Helpers.generate_scan_config_filename_1()

        with open(filename, "w") as file_obj:
            lines = ("Probe Positions 1\n",
                     "1.0 1.0\n",
                     "1.0 1.1\n",
                     "-1")
            for line in lines:
                file_obj.write(line)

        return None



    @staticmethod
    def generate_scan_config_filename_1():
        scan_config_filename = "scan_config_file_1.txt"

        return scan_config_filename



    @staticmethod
    def generate_stem_sim_params_1():
        kwargs = {"stem_system_model_params": \
                  Helpers.generate_stem_system_model_params_1(),
                  "output_params": \
                  Helpers.generate_stem_sim_output_params_1(),
                  "worker_params": \
                  Helpers.generate_worker_params_1(),
                  "skip_validation_and_conversion": \
                  False}
        sim_params = prismatique.stem.sim.Params(**kwargs)

        return sim_params


    
    @staticmethod
    def generate_stem_system_model_params_1():
        kwargs = {"sample_specification": \
                  Helpers.generate_sample_specification_1(),
                  "probe_model_params": \
                  Helpers.generate_probe_model_params_2(),
                  "specimen_tilt": \
                  (0, 0),
                  "scan_config": \
                  Helpers.generate_scan_config_1(),
                  "skip_validation_and_conversion": \
                  False}
        stem_system_model_params = prismatique.stem.system.ModelParams(**kwargs)

        return stem_system_model_params



    @staticmethod
    def generate_sample_specification_1():
        kwargs = {"atomic_coords_filename": \
                  Helpers.generate_atomic_coords_filename_1(),
                  "unit_cell_tiling": \
                  Helpers.generate_unit_cell_tiling_1(),
                  "discretization_params": \
                  Helpers.generate_discretization_params_1(),
                  "atomic_potential_extent": \
                  Helpers.generate_atomic_potential_extent_1(),
                  "thermal_params": \
                  Helpers.generate_thermal_params_1(),
                  "skip_validation_and_conversion": \
                  False}
        sample_specification = prismatique.sample.ModelParams(**kwargs)

        return sample_specification



    @staticmethod
    def generate_unit_cell_tiling_1():
        unit_cell_tiling = (1, 1, 2)

        return unit_cell_tiling



    @staticmethod
    def generate_discretization_params_1():
        kwargs = {"z_supersampling": \
                  8,
                  "sample_supercell_reduced_xy_dims_in_pixels": \
                  (16, 16),
                  "interpolation_factors": \
                  (1, 1),
                  "num_slices": \
                  2,
                  "skip_validation_and_conversion": \
                  False}
        discretization_params = prismatique.discretization.Params(**kwargs)

        return discretization_params



    @staticmethod
    def generate_thermal_params_1():
        kwargs = {"enable_thermal_effects": \
                  True,
                  "num_frozen_phonon_configs_per_subset": \
                  1,
                  "num_subsets": \
                  2,
                  "rng_seed": \
                  None,
                  "skip_validation_and_conversion": \
                  False}
        thermal_params = prismatique.thermal.Params(**kwargs)

        return thermal_params



    @staticmethod
    def generate_probe_model_params_1():
        kwargs = {"lens_model_params": Helpers.generate_lens_model_params_1(),
                  "defocal_offset_supersampling": 1,
                  "convergence_semiangle": 15,
                  "skip_validation_and_conversion": False}
        probe_model_params = embeam.stem.probe.ModelParams(**kwargs)
        
        return probe_model_params



    @staticmethod
    def generate_lens_model_params_1():
        coherent_aberrations = Helpers.generate_coherent_aberrations_1()
        
        kwargs = {"coherent_aberrations": coherent_aberrations,
                  "chromatic_aberration_coef": 0,
                  "mean_current": 1,
                  "std_dev_current": 0,
                  "skip_validation_and_conversion": False}
        lens_model_params = embeam.lens.ModelParams(**kwargs)
        
        return lens_model_params



    @staticmethod
    def generate_coherent_aberrations_1():
        coherent_aberrations = (Helpers.generate_defocus_aberration_1(),
                                Helpers.generate_spherical_aberration_1(),)
        
        return coherent_aberrations



    @staticmethod
    def generate_spherical_aberration_1():
        spherical_aberration = embeam.coherent.Aberration(m=4, 
                                                          n=0, 
                                                          C_mag=0.1, 
                                                          C_ang=0)
        
        return spherical_aberration



    @staticmethod
    def generate_defocus_aberration_1():
        defocus_aberration = embeam.coherent.Aberration(m=2, 
                                                        n=0, 
                                                        C_mag=0.1, 
                                                        C_ang=0)
        
        return defocus_aberration


    
    @staticmethod
    def generate_scan_config_1():
        scan_config = ((1.0, 1.0),
                       (1.0, 1.1))

        return scan_config



    @staticmethod
    def generate_stem_sim_output_params_1():
        kwargs = {"base_params": \
                  Helpers.generate_base_output_params_1(),
                  "alg_specific_params": \
                  Helpers.generate_alg_specific_params_1(),
                  "skip_validation_and_conversion": \
                  False}
        stem_sim_output_params = prismatique.stem.output.Params(**kwargs)

        return stem_sim_output_params



    @staticmethod
    def generate_base_output_params_1():
        kwargs = {"output_dirname": \
                  Helpers.generate_output_dirname_1(),
                  "max_data_size": \
                  2*10**9,
                  "cbed_params": \
                  Helpers.generate_cbed_params_1(),
                  "radial_step_size_for_3d_stem": \
                  1,
                  "radial_range_for_2d_stem": \
                  (0, 1),
                  "save_com": \
                  True,
                  "save_potential_slices": \
                  False,
                  "skip_validation_and_conversion": \
                  False}
        base_output_params = prismatique.stem.output.base.Params(**kwargs)

        return base_output_params



    @staticmethod
    def generate_output_dirname_1():
        output_dirname = "sim_output_files"

        return output_dirname



    @staticmethod
    def generate_cbed_params_1():
        kwargs = {"postprocessing_seq": \
                  Helpers.generate_postprocessing_seq_1(),
                  "avg_num_electrons_per_postprocessed_dp": \
                  1,
                  "apply_shot_noise": \
                  True,
                  "save_wavefunctions": \
                  True,
                  "save_final_intensity": \
                  True,
                  "skip_validation_and_conversion": \
                  False}
        cbed_params = prismatique.cbed.Params(**kwargs)

        return cbed_params



    @staticmethod
    def generate_postprocessing_seq_1():
        kwargs = {"block_dims": (1, 1),
                  "padding_const": 0,
                  "downsample_mode": "mean",
                  "skip_validation_and_conversion": False}
        postprocessing_step_1 = empix.OptionalDownsamplingParams(**kwargs)

        kwargs = {"center": None, "window_dims": None}
        postprocessing_step_2 = empix.OptionalCroppingParams(**kwargs)

        kwargs = {"new_signal_space_sizes": None,
                  "new_signal_space_scales": None,
                  "new_signal_space_offsets": None}
        postprocessing_step_3 = empix.OptionalResamplingParams(**kwargs)
        
        postprocessing_seq = (postprocessing_step_1,
                              postprocessing_step_2,
                              postprocessing_step_3)

        return postprocessing_seq



    @staticmethod
    def generate_alg_specific_params_1():
        kwargs = {"enable_S_matrix_refocus": \
                  False,
                  "save_S_matrices": \
                  True,
                  "skip_validation_and_conversion": \
                  False}
        alg_specific_params = prismatique.stem.output.prism.Params(**kwargs)

        return alg_specific_params



    @staticmethod
    def generate_worker_params_1():
        worker_params = prismatique.worker.Params()

        return worker_params



    @staticmethod
    def generate_dummy_sim_output_files(sim_params):
        output_dirname = Helpers.get_output_dirname_from_sim_params(sim_params)
        pathlib.Path(output_dirname).mkdir(parents=True, exist_ok=True)

        module_alias = (prismatique.stem.sim
                        if isinstance(sim_params, prismatique.stem.sim.Params)
                        else prismatique.hrtem.sim)
        module_alias._initialize_output_files(sim_params)

        Helpers.generate_dummy_sample_specification_output_files(sim_params)

        filename = (output_dirname + "/stem_sim_params.json"
                    if isinstance(sim_params, prismatique.stem.sim.Params)
                    else output_dirname + "/hrtem_sim_params.json")

        kwargs = {"filename": filename, "overwrite": True}
        sim_params.dump(**kwargs)

        sample_specification = \
            Helpers.get_sample_specification_from_sim_params(sim_params)

        kwargs = {"filename": output_dirname + "/sample_model_params.json",
                  "overwrite": True}
        sample_specification.dump(**kwargs)
        
        return None



    @staticmethod
    def generate_dummy_sample_specification_output_files(sim_params):
        sample_specification = \
            Helpers.get_sample_specification_from_sim_params(sim_params)
        output_dirname = \
            Helpers.get_output_dirname_from_sim_params(sim_params)

        kwargs = \
            {"sample_specification": sample_specification}
        sample_specification_output_metadata = \
            Helpers.generate_sample_specification_output_metadata(**kwargs)
        potential_slice_output_data = \
            Helpers.generate_dummy_potential_slice_output_data(**kwargs)

        output_datasets = (potential_slice_output_data,)
        
        if isinstance(sim_params, prismatique.stem.sim.Params):
            kwargs["probe_model_params"] = \
                Helpers.get_probe_model_params_from_sim_params(sim_params)
            S_matrix_output_data = \
                Helpers.generate_dummy_S_matrix_output_data(**kwargs)
            
            output_datasets += (S_matrix_output_data,)

        for sample_specification_output_data in output_datasets:
            method_name = "write_sample_specification_output_metadata_to_files"
            method_alias = getattr(Helpers, method_name)
            kwargs = {"sample_specification_output_data": \
                      sample_specification_output_data,
                      "output_dirname": \
                      output_dirname,
                      "sample_specification_output_metadata": \
                      sample_specification_output_metadata}
            method_alias(**kwargs)

            method_name = "write_sample_specification_output_data_to_files"
            method_alias = getattr(Helpers, method_name)
            del kwargs["sample_specification_output_metadata"]
            method_alias(**kwargs)

        return None



    @staticmethod
    def get_sample_specification_from_sim_params(sim_params):
        if isinstance(sim_params, prismatique.stem.sim.Params):
            stem_system_model_params = \
                sim_params.core_attrs["stem_system_model_params"]
            sample_specification = \
                stem_system_model_params.core_attrs["sample_specification"]
        else:
            hrtem_system_model_params = \
                sim_params.core_attrs["hrtem_system_model_params"]
            sample_specification = \
                hrtem_system_model_params.core_attrs["sample_specification"]

        return sample_specification



    @staticmethod
    def get_output_dirname_from_sim_params(sim_params):
        if isinstance(sim_params, prismatique.stem.sim.Params):
            output_params = sim_params.core_attrs["output_params"]
            base_output_params = output_params.core_attrs["base_params"]
            output_dirname = base_output_params.core_attrs["output_dirname"]
        else:
            output_params = sim_params.core_attrs["output_params"]
            output_dirname = output_params.core_attrs["output_dirname"]

        return output_dirname



    @staticmethod
    def generate_sample_specification_output_metadata(sample_specification):
        unit_cell_tiling = sample_specification.core_attrs["unit_cell_tiling"]

        discretization_params = \
            sample_specification.core_attrs["discretization_params"]
        interpolation_factors = \
            discretization_params.core_attrs["interpolation_factors"]

        func_alias = prismatique.sample._supercell_dims
        sample_supercell_dims = func_alias(sample_specification)
        
        sample_unit_cell_dims = tuple((np.array(sample_supercell_dims)
                                       / np.array(unit_cell_tiling)).tolist())

        func_alias = prismatique.sample._supercell_lateral_pixel_size
        sample_supercell_lateral_pixel_size = func_alias(sample_specification)

        func_alias = prismatique.sample._supercell_slice_thickness
        sample_supercell_slice_thickness = func_alias(sample_specification)

        sample_specification_output_metadata = \
            {"c": sample_unit_cell_dims,
             "t": unit_cell_tiling,
             "px": sample_supercell_lateral_pixel_size[0],
             "py": sample_supercell_lateral_pixel_size[1],
             "s": sample_supercell_slice_thickness,
             "fx": interpolation_factors[0],
             "fy": interpolation_factors[1]}

        return sample_specification_output_metadata



    @staticmethod
    def generate_dummy_potential_slice_output_data(sample_specification):
        thermal_params = \
            sample_specification.core_attrs["thermal_params"]
        num_frozen_phonon_configs_per_subset = \
            thermal_params.core_attrs["num_frozen_phonon_configs_per_subset"]
        num_subsets = \
            thermal_params.core_attrs["num_subsets"]

        func_alias = prismatique.sample._supercell_xy_dims_in_pixels
        sample_supercell_xy_dims_in_pixels = func_alias(sample_specification)

        func_alias = prismatique.sample._num_slices
        num_slices = func_alias(sample_specification)

        potential_slice_output_data_shape = \
            (num_subsets,
             num_frozen_phonon_configs_per_subset,
             sample_supercell_xy_dims_in_pixels[0],
             sample_supercell_xy_dims_in_pixels[1],
             num_slices)
        potential_slice_output_data = \
            np.zeros(potential_slice_output_data_shape, dtype=np.float32)

        return potential_slice_output_data



    @staticmethod
    def get_probe_model_params_from_sim_params(sim_params):
        stem_system_model_params = \
            sim_params.core_attrs["stem_system_model_params"]
        probe_model_params = \
            stem_system_model_params.core_attrs["probe_model_params"]

        return probe_model_params



    @staticmethod
    def generate_dummy_S_matrix_output_data(sample_specification,
                                            probe_model_params):
        thermal_params = \
            sample_specification.core_attrs["thermal_params"]
        num_frozen_phonon_configs_per_subset = \
            thermal_params.core_attrs["num_frozen_phonon_configs_per_subset"]
        num_subsets = \
            thermal_params.core_attrs["num_subsets"]

        discretization_params = \
            sample_specification.core_attrs["discretization_params"]
        interpolation_factors = \
            discretization_params.core_attrs["interpolation_factors"]

        func_alias = prismatique.sample._supercell_xy_dims_in_pixels
        sample_supercell_xy_dims_in_pixels = func_alias(sample_specification)

        func_alias = prismatique.sample._S_matrix_k_xy_vectors
        S_matrix_k_xy_vectors = func_alias(sample_specification,
                                           probe_model_params)

        S_matrix_output_data_shape = \
            (num_subsets,
             num_frozen_phonon_configs_per_subset,
             sample_supercell_xy_dims_in_pixels[0]//interpolation_factors[0]//2,
             sample_supercell_xy_dims_in_pixels[1]//interpolation_factors[1]//2,
             len(S_matrix_k_xy_vectors))
        S_matrix_output_data = \
            np.zeros(S_matrix_output_data_shape, dtype=np.complex64)

        return S_matrix_output_data



    @staticmethod
    def write_sample_specification_output_metadata_to_files(
            sample_specification_output_data,
            output_dirname,
            sample_specification_output_metadata):
        num_subsets = \
            sample_specification_output_data.shape[0]
        unformatted_basename = \
            ("potential_slices_of_subset_{}.h5"
             if (sample_specification_output_data.dtype == np.float32)
             else "S_matrices_of_subset_{}.h5")

        for subset_idx in range(num_subsets):
            unformatted_filename = output_dirname + "/" + unformatted_basename
            filename = unformatted_filename.format(subset_idx)
            
            path_in_file = ("/4DSTEM_simulation/metadata/metadata_0/original"
                            "/simulation_parameters")

            kwargs = {"filename": filename, "path_in_file": path_in_file}
            obj_id = h5pywrappers.obj.ID(**kwargs)

            kwargs = {"group": None, "group_id": obj_id, "write_mode": "a"}
            h5pywrappers.group.save(**kwargs)

            for hdf5_attr_name in sample_specification_output_metadata:
                hdf5_attr = sample_specification_output_metadata[hdf5_attr_name]

                kwargs = {"obj_id": obj_id, "attr_name": hdf5_attr_name}
                attr_id = h5pywrappers.attr.ID(**kwargs)

                kwargs = {"attr": hdf5_attr,
                          "attr_id": attr_id,
                          "write_mode": "a"}
                h5pywrappers.attr.save(**kwargs)

        return None



    @staticmethod
    def write_sample_specification_output_data_to_files(
            output_dirname,
            sample_specification_output_data):
        num_subsets = \
            sample_specification_output_data.shape[0]
        num_atomic_configs_per_subset = \
            sample_specification_output_data.shape[1]
        unformatted_basename = \
            ("potential_slices_of_subset_{}.h5"
             if (sample_specification_output_data.dtype == np.float32)
             else "S_matrices_of_subset_{}.h5")
        unformatted_path_in_file = \
            ("/4DSTEM_simulation/data/realslices/ppotential_fp{}/data"
             if (sample_specification_output_data.dtype == np.float32)
             else "/4DSTEM_simulation/data/realslices/smatrix_fp{}/data")

        for subset_idx in range(num_subsets):
            unformatted_filename = output_dirname + "/" + unformatted_basename
            filename = unformatted_filename.format(subset_idx)
            
            for atomic_config_idx in range(num_atomic_configs_per_subset):
                args = (str(atomic_config_idx).rjust(4, "0"),)
                path_in_file = unformatted_path_in_file.format(*args)

                kwargs = {"filename": filename, "path_in_file": path_in_file}
                dataset_id = h5pywrappers.obj.ID(**kwargs)

                multi_dim_slice = (subset_idx, atomic_config_idx)
                dataset = sample_specification_output_data[multi_dim_slice]

                kwargs = {"dataset": dataset,
                          "dataset_id": dataset_id,
                          "write_mode": "a"}
                h5pywrappers.dataset.save(**kwargs)

        return None



    @staticmethod
    def generate_invalid_sim_output_file_1():
        output_dirname = Helpers.generate_output_dirname_1()
        pathlib.Path(output_dirname).mkdir(parents=True, exist_ok=True)

        filename = Helpers.generate_invalid_sim_output_filename_1()

        kwargs = {"filename": filename, "path_in_file": "/data/foobar"}
        obj_id = h5pywrappers.obj.ID(**kwargs)

        kwargs = {"group": None, "group_id": obj_id, "write_mode": "a"}
        h5pywrappers.group.save(**kwargs)

        path_in_files = ("/metadata/k_x",
                         "/metadata/k_y",
                         "/data/4D_STEM/intensity_DPs")
        datasets = (np.array([0.0, 1.0]),
                    np.array([0.0, 1.0]),
                    np.zeros((3, 3)))
        zip_obj = zip(path_in_files, datasets)

        for path_in_file, dataset in zip_obj:
            kwargs = {"filename": filename, "path_in_file": path_in_file}
            dataset_id = h5pywrappers.obj.ID(**kwargs)

            kwargs = {"dataset": dataset,
                      "dataset_id": dataset_id,
                      "write_mode": "a"}
            h5pywrappers.dataset.save(**kwargs)
        
        return None



    @staticmethod
    def generate_invalid_sim_output_filename_1():
        output_dirname = Helpers.generate_output_dirname_1()
        filename = output_dirname + "/invalid_output_file_1.h5"
        
        return filename
    


    @staticmethod
    def generate_stem_sim_params_2():
        kwargs = {"stem_system_model_params": \
                  Helpers.generate_stem_system_model_params_2(),
                  "output_params": \
                  Helpers.generate_stem_sim_output_params_2(),
                  "worker_params": \
                  Helpers.generate_worker_params_1(),
                  "skip_validation_and_conversion": \
                  False}
        sim_params = prismatique.stem.sim.Params(**kwargs)

        return sim_params



    @staticmethod
    def generate_stem_system_model_params_2():
        kwargs = {"sample_specification": \
                  Helpers.generate_sample_specification_2(),
                  "probe_model_params": \
                  Helpers.generate_probe_model_params_1(),
                  "specimen_tilt": \
                  (0, 0),
                  "scan_config": \
                  Helpers.generate_scan_config_2(),
                  "skip_validation_and_conversion": \
                  False}
        stem_system_model_params = prismatique.stem.system.ModelParams(**kwargs)

        return stem_system_model_params



    @staticmethod
    def generate_sample_specification_2():
        kwargs = {"atomic_coords_filename": \
                  Helpers.generate_atomic_coords_filename_2(),
                  "unit_cell_tiling": \
                  Helpers.generate_unit_cell_tiling_2(),
                  "discretization_params": \
                  Helpers.generate_discretization_params_1(),
                  "atomic_potential_extent": \
                  Helpers.generate_atomic_potential_extent_1(),
                  "thermal_params": \
                  Helpers.generate_thermal_params_2(),
                  "skip_validation_and_conversion": \
                  False}
        sample_specification = prismatique.sample.ModelParams(**kwargs)

        return sample_specification



    @staticmethod
    def generate_unit_cell_tiling_2():
        unit_cell_tiling = (1, 1, 1)

        return unit_cell_tiling



    @staticmethod
    def generate_thermal_params_2():
        kwargs = {"enable_thermal_effects": \
                  True,
                  "num_frozen_phonon_configs_per_subset": \
                  2,
                  "num_subsets": \
                  1,
                  "rng_seed": \
                  1,
                  "skip_validation_and_conversion": \
                  False}
        thermal_params = prismatique.thermal.Params(**kwargs)

        return thermal_params



    @staticmethod
    def generate_scan_config_2():
        kwargs = {"step_size": (1.0, 1.0),
                  "window": (0.499, 0.5, 0.499, 0.5),
                  "skip_validation_and_conversion": False}
        scan_config = prismatique.scan.rectangular.Params(**kwargs)

        return scan_config



    @staticmethod
    def generate_stem_sim_output_params_2():
        kwargs = {"base_params": \
                  Helpers.generate_base_output_params_2(),
                  "alg_specific_params": \
                  Helpers.generate_alg_specific_params_2(),
                  "skip_validation_and_conversion": \
                  False}
        stem_sim_output_params = prismatique.stem.output.Params(**kwargs)

        return stem_sim_output_params



    @staticmethod
    def generate_base_output_params_2():
        kwargs = {"output_dirname": \
                  Helpers.generate_output_dirname_1(),
                  "max_data_size": \
                  2*10**9,
                  "cbed_params": \
                  Helpers.generate_cbed_params_2(),
                  "radial_step_size_for_3d_stem": \
                  0,
                  "radial_range_for_2d_stem": \
                  (0, 0),
                  "save_com": \
                  False,
                  "save_potential_slices": \
                  True,
                  "skip_validation_and_conversion": \
                  False}
        base_output_params = prismatique.stem.output.base.Params(**kwargs)

        return base_output_params



    @staticmethod
    def generate_cbed_params_2():
        kwargs = {"postprocessing_seq": \
                  tuple(),
                  "avg_num_electrons_per_postprocessed_dp": \
                  1,
                  "apply_shot_noise": \
                  False,
                  "save_wavefunctions": \
                  False,
                  "save_final_intensity": \
                  False,
                  "skip_validation_and_conversion": \
                  False}
        cbed_params = prismatique.cbed.Params(**kwargs)

        return cbed_params



    @staticmethod
    def generate_alg_specific_params_2():
        kwargs = {"enable_S_matrix_refocus": \
                  False,
                  "save_S_matrices": \
                  False,
                  "skip_validation_and_conversion": \
                  False}
        alg_specific_params = prismatique.stem.output.prism.Params(**kwargs)

        return alg_specific_params



    @staticmethod
    def generate_stem_sim_params_3():
        kwargs = {"stem_system_model_params": \
                  Helpers.generate_stem_system_model_params_3(),
                  "output_params": \
                  Helpers.generate_stem_sim_output_params_3(),
                  "worker_params": \
                  Helpers.generate_worker_params_1(),
                  "skip_validation_and_conversion": \
                  False}
        sim_params = prismatique.stem.sim.Params(**kwargs)

        return sim_params



    @staticmethod
    def generate_stem_system_model_params_3():
        kwargs = {"sample_specification": \
                  Helpers.generate_sample_specification_2(),
                  "probe_model_params": \
                  Helpers.generate_probe_model_params_1(),
                  "specimen_tilt": \
                  (0, 0),
                  "scan_config": \
                  Helpers.generate_scan_config_3(),
                  "skip_validation_and_conversion": \
                  False}
        stem_system_model_params = prismatique.stem.system.ModelParams(**kwargs)

        return stem_system_model_params



    @staticmethod
    def generate_scan_config_3():
        scan_config = Helpers.generate_scan_config_filename_1()

        return scan_config



    @staticmethod
    def generate_stem_sim_output_params_3():
        kwargs = {"base_params": \
                  Helpers.generate_base_output_params_3(),
                  "alg_specific_params": \
                  Helpers.generate_alg_specific_params_1(),
                  "skip_validation_and_conversion": \
                  False}
        stem_sim_output_params = prismatique.stem.output.Params(**kwargs)

        return stem_sim_output_params



    @staticmethod
    def generate_base_output_params_3():
        kwargs = {"output_dirname": \
                  Helpers.generate_output_dirname_1(),
                  "max_data_size": \
                  1,
                  "cbed_params": \
                  Helpers.generate_cbed_params_2(),
                  "radial_step_size_for_3d_stem": \
                  0,
                  "radial_range_for_2d_stem": \
                  (0, 0),
                  "save_com": \
                  False,
                  "save_potential_slices": \
                  False,
                  "skip_validation_and_conversion": \
                  False}
        base_output_params = prismatique.stem.output.base.Params(**kwargs)

        return base_output_params



    @staticmethod
    def generate_stem_sim_params_4():
        kwargs = {"stem_system_model_params": \
                  Helpers.generate_stem_system_model_params_4(),
                  "output_params": \
                  Helpers.generate_stem_sim_output_params_4(),
                  "worker_params": \
                  Helpers.generate_worker_params_1(),
                  "skip_validation_and_conversion": \
                  False}
        sim_params = prismatique.stem.sim.Params(**kwargs)

        return sim_params



    @staticmethod
    def generate_stem_system_model_params_4():
        kwargs = {"sample_specification": \
                  Helpers.generate_sample_specification_2(),
                  "probe_model_params": \
                  Helpers.generate_probe_model_params_2(),
                  "specimen_tilt": \
                  (0, 0),
                  "scan_config": \
                  Helpers.generate_scan_config_2(),
                  "skip_validation_and_conversion": \
                  False}
        stem_system_model_params = prismatique.stem.system.ModelParams(**kwargs)

        return stem_system_model_params



    @staticmethod
    def generate_probe_model_params_2():
        kwargs = {"lens_model_params": Helpers.generate_lens_model_params_2(),
                  "defocal_offset_supersampling": 2,
                  "convergence_semiangle": 20,
                  "skip_validation_and_conversion": False}
        probe_model_params = embeam.stem.probe.ModelParams(**kwargs)
        
        return probe_model_params



    @staticmethod
    def generate_lens_model_params_2():
        coherent_aberrations = Helpers.generate_coherent_aberrations_2()
        
        kwargs = {"coherent_aberrations": coherent_aberrations,
                  "chromatic_aberration_coef": 1,
                  "mean_current": 1,
                  "std_dev_current": 0.1,
                  "skip_validation_and_conversion": False}
        lens_model_params = embeam.lens.ModelParams(**kwargs)
        
        return lens_model_params



    @staticmethod
    def generate_coherent_aberrations_2():
        coherent_aberrations = (Helpers.generate_spherical_aberration_1(),)
        
        return coherent_aberrations



    @staticmethod
    def generate_stem_sim_output_params_4():
        kwargs = {"base_params": \
                  Helpers.generate_base_output_params_4(),
                  "alg_specific_params": \
                  Helpers.generate_alg_specific_params_3(),
                  "skip_validation_and_conversion": \
                  False}
        stem_sim_output_params = prismatique.stem.output.Params(**kwargs)

        return stem_sim_output_params



    @staticmethod
    def generate_base_output_params_4():
        kwargs = {"output_dirname": \
                  Helpers.generate_output_dirname_1(),
                  "max_data_size": \
                  2*10**9,
                  "cbed_params": \
                  Helpers.generate_cbed_params_3(),
                  "radial_step_size_for_3d_stem": \
                  1,
                  "radial_range_for_2d_stem": \
                  (0, 1),
                  "save_com": \
                  True,
                  "save_potential_slices": \
                  False,
                  "skip_validation_and_conversion": \
                  False}
        base_output_params = prismatique.stem.output.base.Params(**kwargs)

        return base_output_params



    @staticmethod
    def generate_cbed_params_3():
        kwargs = {"postprocessing_seq": \
                  tuple(),
                  "avg_num_electrons_per_postprocessed_dp": \
                  1,
                  "apply_shot_noise": \
                  False,
                  "save_wavefunctions": \
                  False,
                  "save_final_intensity": \
                  True,
                  "skip_validation_and_conversion": \
                  False}
        cbed_params = prismatique.cbed.Params(**kwargs)

        return cbed_params



    @staticmethod
    def generate_alg_specific_params_3():
        kwargs = {"num_slices_per_output": \
                  1,
                  "z_start_output": \
                  float("inf"),
                  "skip_validation_and_conversion": \
                  False}
        alg_specific_params = prismatique.stem.output.multislice.Params(**kwargs)

        return alg_specific_params



    @staticmethod
    def generate_stem_sim_params_5():
        kwargs = {"stem_system_model_params": \
                  Helpers.generate_stem_system_model_params_5(),
                  "output_params": \
                  Helpers.generate_stem_sim_output_params_1(),
                  "worker_params": \
                  Helpers.generate_worker_params_1(),
                  "skip_validation_and_conversion": \
                  False}
        sim_params = prismatique.stem.sim.Params(**kwargs)

        return sim_params



    @staticmethod
    def generate_stem_system_model_params_5():
        kwargs = {"sample_specification": \
                  Helpers.generate_sample_specification_3(),
                  "probe_model_params": \
                  Helpers.generate_probe_model_params_2(),
                  "specimen_tilt": \
                  (0, 0),
                  "scan_config": \
                  Helpers.generate_scan_config_2(),
                  "skip_validation_and_conversion": \
                  False}
        stem_system_model_params = prismatique.stem.system.ModelParams(**kwargs)

        return stem_system_model_params



    @staticmethod
    def generate_sample_specification_3():
        output_dirname = Helpers.generate_output_dirname_2()

        filename = output_dirname + "/S_matrices_of_subset_0.h5"

        kwargs = {"filenames": (filename,),
                  "skip_validation_and_conversion": False}
        sample_specification = prismatique.sample.SMatrixSubsetIDs(**kwargs)

        return sample_specification



    @staticmethod
    def generate_output_dirname_2():
        output_dirname = "extra_sim_output_files"

        return output_dirname



    @staticmethod
    def generate_sample_specification_4():
        kwargs = {"atomic_coords_filename": \
                  Helpers.generate_atomic_coords_filename_1(),
                  "unit_cell_tiling": \
                  Helpers.generate_unit_cell_tiling_1(),
                  "discretization_params": \
                  Helpers.generate_discretization_params_2(),
                  "atomic_potential_extent": \
                  Helpers.generate_atomic_potential_extent_1(),
                  "thermal_params": \
                  Helpers.generate_thermal_params_1(),
                  "skip_validation_and_conversion": \
                  False}
        sample_specification = prismatique.sample.ModelParams(**kwargs)

        return sample_specification



    @staticmethod
    def generate_discretization_params_2():
        kwargs = {"z_supersampling": \
                  8,
                  "sample_supercell_reduced_xy_dims_in_pixels": \
                  (16, 16),
                  "interpolation_factors": \
                  (1, 1),
                  "num_slices": \
                  10,
                  "skip_validation_and_conversion": \
                  False}
        discretization_params = prismatique.discretization.Params(**kwargs)

        return discretization_params



    @staticmethod
    def generate_hrtem_sim_params_1():
        kwargs ={"hrtem_system_model_params": \
                  Helpers.generate_hrtem_system_model_params_1(),
                  "output_params": \
                  Helpers.generate_hrtem_sim_output_params_1(),
                  "worker_params": \
                  Helpers.generate_worker_params_1(),
                  "skip_validation_and_conversion": \
                  False}
        sim_params = prismatique.hrtem.sim.Params(**kwargs)

        return sim_params


    
    @staticmethod
    def generate_hrtem_system_model_params_1():
        kwargs = \
            {"sample_specification": \
             Helpers.generate_sample_specification_1(),
             "gun_model_params": \
             None,
             "lens_model_params": \
             Helpers.generate_lens_model_params_1(),
             "tilt_params": \
             Helpers.generate_tilt_params_1(),
             "objective_aperture_params": \
             Helpers.generate_objective_aperture_params_1(),
             "defocal_offset_supersampling": \
             1,
             "skip_validation_and_conversion": \
             False}
        hrtem_system_model_params = \
            prismatique.hrtem.system.ModelParams(**kwargs)

        return hrtem_system_model_params



    @staticmethod
    def generate_tilt_params_1():
        tilt_params = None

        return tilt_params



    @staticmethod
    def generate_objective_aperture_params_1():
        objective_aperture_params = None

        return objective_aperture_params



    @staticmethod
    def generate_hrtem_sim_output_params_1():
        kwargs = {"output_dirname": \
                  Helpers.generate_output_dirname_1(),
                  "max_data_size": \
                  2*10**9,
                  "image_params": \
                  Helpers.generate_image_params_1(),
                  "save_potential_slices": \
                  False,
                  "skip_validation_and_conversion": \
                  False}
        hrtem_sim_output_params = prismatique.hrtem.output.Params(**kwargs)

        return hrtem_sim_output_params



    @staticmethod
    def generate_image_params_1():
        kwargs = {"postprocessing_seq": \
                  Helpers.generate_postprocessing_seq_1(),
                  "avg_num_electrons_per_postprocessed_image": \
                  1,
                  "apply_shot_noise": \
                  True,
                  "save_wavefunctions": \
                  True,
                  "save_final_intensity": \
                  True,
                  "skip_validation_and_conversion": \
                  False}
        image_params = prismatique.hrtem.image.Params(**kwargs)

        return image_params



    @staticmethod
    def generate_hrtem_sim_params_2():
        kwargs = {"hrtem_system_model_params": \
                  Helpers.generate_hrtem_system_model_params_2(),
                  "output_params": \
                  Helpers.generate_hrtem_sim_output_params_2(),
                  "worker_params": \
                  Helpers.generate_worker_params_1(),
                  "skip_validation_and_conversion": \
                  False}
        sim_params = prismatique.hrtem.sim.Params(**kwargs)

        return sim_params



    @staticmethod
    def generate_hrtem_system_model_params_2():
        kwargs = \
            {"sample_specification": \
             Helpers.generate_sample_specification_2(),
             "gun_model_params": \
             None,
             "lens_model_params": \
             Helpers.generate_lens_model_params_1(),
             "tilt_params": \
             Helpers.generate_tilt_params_2(),
             "objective_aperture_params": \
             Helpers.generate_objective_aperture_params_1(),
             "defocal_offset_supersampling": \
             1,
             "skip_validation_and_conversion": \
             False}
        hrtem_system_model_params = \
            prismatique.hrtem.system.ModelParams(**kwargs)

        return hrtem_system_model_params



    @staticmethod
    def generate_tilt_params_2():
        kwargs = {"offset": (0, 0),
                  "window": (0, 0, 0, 0),
                  "spread": 0.0001,
                  "skip_validation_and_conversion": False}
        tilt_params = prismatique.tilt.Params(**kwargs)

        return tilt_params



    @staticmethod
    def generate_hrtem_sim_output_params_2():
        kwargs = {"output_dirname": \
                  Helpers.generate_output_dirname_1(),
                  "max_data_size": \
                  2*10**9,
                  "image_params": \
                  Helpers.generate_image_params_2(),
                  "save_potential_slices": \
                  False,
                  "skip_validation_and_conversion": \
                  False}
        hrtem_sim_output_params = prismatique.hrtem.output.Params(**kwargs)

        return hrtem_sim_output_params



    @staticmethod
    def generate_image_params_2():
        kwargs = {"postprocessing_seq": \
                  tuple(),
                  "avg_num_electrons_per_postprocessed_image": \
                  1,
                  "apply_shot_noise": \
                  False,
                  "save_wavefunctions": \
                  False,
                  "save_final_intensity": \
                  False,
                  "skip_validation_and_conversion": \
                  False}
        image_params = prismatique.hrtem.image.Params(**kwargs)

        return image_params



    @staticmethod
    def generate_hrtem_sim_params_3():
        kwargs = {"hrtem_system_model_params": \
                  Helpers.generate_hrtem_system_model_params_2(),
                  "output_params": \
                  Helpers.generate_hrtem_sim_output_params_3(),
                  "worker_params": \
                  Helpers.generate_worker_params_1(),
                  "skip_validation_and_conversion": \
                  False}
        sim_params = prismatique.hrtem.sim.Params(**kwargs)

        return sim_params



    @staticmethod
    def generate_hrtem_sim_output_params_3():
        kwargs = {"output_dirname": \
                  Helpers.generate_output_dirname_1(),
                  "max_data_size": \
                  1,
                  "image_params": \
                  Helpers.generate_image_params_2(),
                  "save_potential_slices": \
                  True,
                  "skip_validation_and_conversion": \
                  False}
        hrtem_sim_output_params = prismatique.hrtem.output.Params(**kwargs)

        return hrtem_sim_output_params



    @staticmethod
    def generate_hrtem_sim_params_4():
        kwargs = {"hrtem_system_model_params": \
                  Helpers.generate_hrtem_system_model_params_3(),
                  "output_params": \
                  Helpers.generate_hrtem_sim_output_params_4(),
                  "worker_params": \
                  Helpers.generate_worker_params_1(),
                  "skip_validation_and_conversion": \
                  False}
        sim_params = prismatique.hrtem.sim.Params(**kwargs)

        return sim_params



    @staticmethod
    def generate_hrtem_system_model_params_3():
        kwargs = \
            {"sample_specification": \
             Helpers.generate_sample_specification_2(),
             "gun_model_params": \
             None,
             "lens_model_params": \
             Helpers.generate_lens_model_params_2(),
             "tilt_params": \
             Helpers.generate_tilt_params_2(),
             "objective_aperture_params": \
             Helpers.generate_objective_aperture_params_1(),
             "defocal_offset_supersampling": \
             2,
             "skip_validation_and_conversion": \
             False}
        hrtem_system_model_params = \
            prismatique.hrtem.system.ModelParams(**kwargs)

        return hrtem_system_model_params



    @staticmethod
    def generate_hrtem_sim_output_params_4():
        kwargs = {"output_dirname": \
                  Helpers.generate_output_dirname_1(),
                  "max_data_size": \
                  2*10**9,
                  "image_params": \
                  Helpers.generate_image_params_3(),
                  "save_potential_slices": \
                  False,
                  "skip_validation_and_conversion": \
                  False}
        hrtem_sim_output_params = prismatique.hrtem.output.Params(**kwargs)

        return hrtem_sim_output_params



    @staticmethod
    def generate_image_params_3():
        kwargs = {"postprocessing_seq": \
                  tuple(),
                  "avg_num_electrons_per_postprocessed_image": \
                  1,
                  "apply_shot_noise": \
                  False,
                  "save_wavefunctions": \
                  False,
                  "save_final_intensity": \
                  True,
                  "skip_validation_and_conversion": \
                  False}
        image_params = prismatique.hrtem.image.Params(**kwargs)

        return image_params



    @staticmethod
    def generate_invalid_atomic_coords_file_1():
        filename = Helpers.generate_atomic_coords_filename_3()

        with open(filename, "w") as file_obj:
            lines = ("Invalid Atomic Coords File 1\n"
                     "5 5 5\n",
                     "8 1.0 1.0 1.0 1.0 0.0\n",
                     "-2")
            for line in lines:
                file_obj.write(line)
        
        return None



    @staticmethod
    def generate_atomic_coords_filename_3():
        atomic_coords_filename = "atomic_coords_file_3.xyz"
        
        return atomic_coords_filename



    @staticmethod
    def generate_invalid_atomic_coords_file_2():
        filename = Helpers.generate_atomic_coords_filename_4()

        with open(filename, "w") as file_obj:
            lines = ("Invalid Atomic Coords File 2\n"
                     "5 5 -5\n",
                     "8 1.0 1.0 1.0 1.0 0.0\n",
                     "-1")
            for line in lines:
                file_obj.write(line)
        
        return None



    @staticmethod
    def generate_atomic_coords_filename_4():
        atomic_coords_filename = "atomic_coords_file_4.xyz"
        
        return atomic_coords_filename



    @staticmethod
    def generate_invalid_atomic_coords_file_3():
        filename = Helpers.generate_atomic_coords_filename_5()

        with open(filename, "w") as file_obj:
            lines = ("Invalid Atomic Coords File 2\n"
                     "5 5 5\n",
                     "200 1.0 1.0 1.0 1.0 0.0\n",
                     "-1")
            for line in lines:
                file_obj.write(line)
        
        return None



    @staticmethod
    def generate_atomic_coords_filename_5():
        atomic_coords_filename = "atomic_coords_file_5.xyz"
        
        return atomic_coords_filename



    @staticmethod
    def generate_invalid_scan_config_file_1():
        filename = Helpers.generate_scan_config_filename_2()

        with open(filename, "w") as file_obj:
            lines = ("Probe Positions 1\n",
                     "1.0 1.0\n",
                     "1.0 1.1\n",
                     "-2")
            for line in lines:
                file_obj.write(line)

        return None



    @staticmethod
    def generate_scan_config_filename_2():
        scan_config_filename = "scan_config_file_2.txt"

        return scan_config_filename



    @staticmethod
    def generate_invalid_scan_config_file_2():
        filename = Helpers.generate_scan_config_filename_3()

        with open(filename, "w") as file_obj:
            lines = ("Probe Positions 1\n",
                     "1.0 1.0\n",
                     "1.0 ABC\n",
                     "-1")
            for line in lines:
                file_obj.write(line)

        return None



    @staticmethod
    def generate_scan_config_filename_3():
        scan_config_filename = "scan_config_file_3.txt"

        return scan_config_filename



@pytest.fixture
def helpers():
    return Helpers
