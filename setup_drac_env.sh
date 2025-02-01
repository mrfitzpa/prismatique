#!/bin/bash

source ./load_drac_modules.sh
virtualenv --no-download ~/$1
source ~/$1/bin/activate
pip install --no-index --upgrade pip
pip install --no-index hyperspy h5py pyFAI pyprismatic-gpu pytest ipympl
pip install git+https://gitlab.com/mrfitzpa/embeam.git
pip install git+https://github.com/mrfitzpa/h5pywrappers.git
