.. _minimal_slurm_job_submission_example_sec:

Digital Research Alliance of Canada job submission example
==========================================================

A minimal example of submitting a ``prismatique`` job on Digital Research
Alliance of Canada (DRAC) can be found in the directory
``<root>/examples/minimal_slurm_job_submission``, wherein the file ``job.py`` is
the job script, and the file ``submit_job.sh`` is the submission script. To
submit the example on DRAC, change into the directory containing the submission
script, and then issue the following command::

  sbatch submit_job.sh

Note that the job script itself is trivial. The purpose of this example is to
simply demonstrate how to submit jobs, hence the more important file in the
aforementioned directory is the submission script ``submit_job.sh``. Statistics
about GPU and CPU usage are saved to the files
``<root>/examples/data/minimal_slurm_job_submission_output/gpu_usage.txt`` and
``<root>/examples/data/minimal_slurm_job_submission_output/time_and_cpu_mem_stats.txt``
respectively, where ``<root>`` is the root of the ``prismatique`` repository.

If you would like to modify the job script and/or the submission script for your
own work, it is recommended that you copy the original scripts and save them
elsewhere outside of the git repository so that the changes made are not tracked
by git.

Contents of the file ``job.py``:

.. literalinclude:: ../../examples/minimal_slurm_job_submission/job.py

Contents of the file ``submit_job.sh``:

.. literalinclude:: ../../examples/minimal_slurm_job_submission/submit_job.sh