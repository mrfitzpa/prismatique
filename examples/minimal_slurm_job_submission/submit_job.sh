#!/bin/bash
#SBATCH --job-name=some_job_name
#SBATCH --account=<account code, e.g. def-ausername>
#SBATCH --nodes=1
#SBATCH --gpus-per-node=v100l:1  # Number of GPU nodes per core and GPU type.
#SBATCH --cpus-per-task=2         # CPU cores/threads
#SBATCH --mem=1G               # memory per node
#SBATCH --time=0-00:30            # time (DD-HH:MM)
#SBATCH --mail-user=<email address to receive job notifications and updates>
#SBATCH --mail-type=ALL

# Path to python script to run, relative to this current file. Path should
# include the ``.py`` extension as well.
rel_path_to_python_script="./job.py"

# Path to output directory.
output_dirname="../data/minimal_slurm_job_submission_output"

mkdir -p ${output_dirname}
rm -f ${output_dirname}/gpu_usage.txt
nvidia-smi dmon -i 0 -s mu -d 5 -o TD > ${output_dirname}/gpu_usage.txt &

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

( /cvmfs/soft.computecanada.ca/gentoo/2020/usr/bin/time -v \
    python ${rel_path_to_python_script} ) \
    2>${output_dirname}/time_and_cpu_mem_stats.txt
