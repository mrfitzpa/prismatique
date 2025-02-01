"""This script creates a virtual environment and install the minimal 
dependencies of ``prismatique``.

See the installation instructions in the file ``docs/INSTALL.rst`` for further
details.
"""



#####################################
## Load libraries/packages/modules ##
#####################################

# For parsing command line arguments.
import argparse

# For removing temporary files.
import os

# For determining the OS.
import sys

# For running terminals commands within this current python script.
import subprocess

# For determining the machine architecture.
import platform

# For determining the hostname.
import socket



# For running conda commands within this current python script.
import conda.cli



############################
## Authorship information ##
############################

__author__       = "Matthew Fitzpatrick"
__copyright__    = "Copyright 2023"
__credits__      = ["Matthew Fitzpatrick"]
__maintainer__   = "Matthew Fitzpatrick"
__email__        = "mrfitzpa@uvic.ca"
__status__       = "Development"



###########################
## Define error messages ##
###########################

_err_msg_1 = \
    ("\n"
     "The correct form of the command should be:\n"
     "\n"
     "    python setup_env.py "
     "--python=<version> --name=<name>\n"
     "\n"
     "where ``<version>`` is the python version number indicating the version "
     "of python that the user wishes to install in the conda environment to be "
     "created, and ``<name>`` is the name of the new conda environment to be "
     "installed.")
_err_msg_2 = \
    ("An invalid python version was specified. Note that the major version "
     "must be equal to 3.")
_err_msg_3 = \
    ("Virtual environment name cannot contain whitespace.")
_err_msg_4 = \
    ("The current operating system of this machine is not supported by "
     "``prismatique``.")



#########################
## Main body of script ##
#########################

if __name__ == "__main__":
    # Parse command line arguments.
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument('--python', default="")
        parser.add_argument('--name', default="prismatique")
        args = parser.parse_args()
    except:
        raise SystemExit(_err_msg_1)
        
        

    # Check whether a valid python version has been specified.
    try:
        python_version_str = args.python
        if python_version_str != "":
            [int(string) for string in python_version_str.split(".")]
    except:
        raise ValueError(_err_msg_2)
        
        

    # Check whether a valid conda environment name has been specified.
    env_name = args.name
    if " " in env_name:
        raise ValueError(_err_msg_3)
        
        

    # Determine operating system of current machine in use.
    if sys.platform.startswith("linux"):
        operating_system = "linux"
    elif sys.platform.startswith("darwin"):
        operating_system = "macOS"
    elif sys.platform.startswith("win"):
        operating_system = "windows"
    else:
        raise SystemExit(_err_msg_4)



    # Check if CUDA compatible GPUs are available. If they are available, then
    # check the CUDA version as well.
    if sys.platform.startswith("linux"):
        cmds = ("nvidia-smi",)
    elif sys.platform.startswith("darwin"):
        cmds = ("",)
    elif sys.platform.startswith("win"):
        cmd_1 = (r"C:\Windows\System32\DriverStore"
                 r"\FileRepository\nvdm*\nvidia-smi.exe")
        cmd_2 = r"C:\Program Files\NVIDIA Corporation\NVSMI\nvidia-smi.exe"
        cmds = (cmd_1, cmd_2)
    for cmd in cmds:
        proc = subprocess.run(cmd, text=True, capture_output=True, shell=True)
        strings = proc.stdout.split()
        if len(strings) > 0:
            idx = [idx for idx in range(len(strings))
                   if strings[idx] == "CUDA"][0]
            cuda_version_str = strings[idx+2]
            if float(cuda_version_str) < 10.2:
                pyprismatic_spec_str = "pyprismatic=2.*=cpu*"
            else:
                pyprismatic_spec_str = "pyprismatic=2.*=gpu*"
                cuda_version_str += ".*"
                break
        else:
            pyprismatic_spec_str = "pyprismatic=2.*=cpu*"



    # The dependencies of ``prismatique`` are installed in subsets to increase
    # the chances of successfully creating the virtual environment with all the
    # necessary dependencies.



    # Create temporary file that lists a subset of ``prismatique`` dependencies.
    lines = []
    if python_version_str != "":
        lines.append("python=="+python_version_str)
    lines.append("ipympl")
    lines.append("pyqt")
    lines.append("conda")
    lines.append("pip")

    temp_txt_file = "_requirements.txt"
    with open(temp_txt_file, "w") as file_obj:
        file_obj.write("\n".join(lines))
    
    
    
    # Create new conda environment with the above ``prismatique`` dependency
    # subset installed. Remove aforementioned temporary file upon completion.
    # conda_cmd = conda.cli.python_api.Commands.CREATE
    try:
        args = ("create",
                "-n", env_name,
                "-c", "conda-forge",
                "-y",
                "--file", temp_txt_file)
        conda.cli.main(*args)
        if os.path.isfile(temp_txt_file):
            os.remove(temp_txt_file)
    except BaseException as err:
        if os.path.isfile(temp_txt_file):
            os.remove(temp_txt_file)
        raise err



    # Install the dependencies from github and gitlab.
    # conda_cmd = conda.cli.python_api.Commands.RUN
    with open("requirements.txt", "r") as file_obj:
        lines = file_obj.readlines()
    for line in lines:
        if "git" in line:
            _, _, url = line.split()
            args = ("run",
                    "-n", env_name,
                    "pip", "install", url)
            conda.cli.main(*args)


        
    # Install the remaining dependencies.
    libraries = [pyprismatic_spec_str]
    if pyprismatic_spec_str == "pyprismatic=2.*=gpu*":
        libraries.append("cudatoolkit=="+cuda_version_str)
    args = ("install",
            "-n", env_name,
            "-c", "conda-forge",
            "-y") + tuple(libraries)
    conda.cli.main(*args)
