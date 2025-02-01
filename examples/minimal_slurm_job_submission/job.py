"""An example job script to be submitted on Digital Research Alliance of Canada
machines.

See the link
https://mrfitzpa.gitlab.io/prismatique/examples/simulation/minimal-slurm-job-submission.html
for a description of the example.

A NOTE BEFORE STARTING
----------------------

To submit the job, change into the directory containing this job script, and
run::

  sbatch submit-job.sh

Note that the job script itself is trivial. The purpose of this example is to
simply demonstrate how to submit jobs, hence the more important file in the
directory containing this current job script is the submission script
``submit_job.sh``. Statistics about GPU and CPU usage are saved to the files
``<root>/examples/data/minimal-slurm-job-submission-output/gpu-usage.txt`` and
``<root>/examples/data/minimal-slurm-job-submission-output/time-and-cpu-mem-stats.txt``
respectively, where ``<root>`` is the root of the ``prismatique`` repository.

If you would like to modify this script for your own work, it is recommended
that you copy the original script and save it elsewhere outside of the git
repository so that the changes made are not tracked by git.

"""



#####################################
## Load libraries/packages/modules ##
#####################################



############################
## Authorship information ##
############################

__author__     = "Matthew Fitzpatrick"
__copyright__  = "Copyright 2023"
__credits__    = ["Matthew Fitzpatrick"]
__maintainer__ = "Matthew Fitzpatrick"
__email__      = "mrfitzpa@uvic.ca"
__status__     = "Development"



#########################
## Main body of script ##
#########################

print("Hello world.")