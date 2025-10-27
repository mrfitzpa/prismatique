# Prismatique

[![Test library](https://github.com/mrfitzpa/prismatique/actions/workflows/test_library.yml/badge.svg)](https://github.com/mrfitzpa/prismatique/actions/workflows/test_library.yml)
[![Code Coverage](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/mrfitzpa/75a9dc87bd52856433f20af68b721f9a/raw/prismatique_coverage_badge.json)](https://github.com/mrfitzpa/prismatique/actions/workflows/measure_code_coverage.yml)
[![Documentation](https://img.shields.io/badge/docs-read-brightgreen)](https://mrfitzpa.github.io/prismatique)
[![PyPi Version](https://img.shields.io/pypi/v/prismatique.svg)](https://pypi.org/project/prismatique)
[![Conda-Forge Version](https://img.shields.io/conda/vn/conda-forge/prismatique.svg)](https://anaconda.org/conda-forge/prismatique)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

`prismatique` is a Python library that functions essentially as a wrapper to the
Python library `pyprismatic`, which itself is a thin wrapper to `prismatic`, a
CUDA/C++ package for fast image simulations in scanning transmission electron
microscopy and high-resolution transmission electron microscopy. You can find
more information about `pyprismatic` and `prismatic`
[here](https://prism-em.com).

Visit the `prismatique` [website](https://mrfitzpa.github.io/prismatique) for a
web version of the installation instructions, the reference guide, and the
examples archive.

The source code can be found in the [`prismatique` GitHub
repository](https://github.com/mrfitzpa/prismatique).



## Table of contents

- [Instructions for installing and uninstalling
  `prismatique`](#instructions-for-installing-and-uninstalling-prismatique)
  - [Installing `prismatique`](#installing-prismatique)
    - [Installing `prismatique` using
      `pip`](#installing-prismatique-using-pip)
    - [Installing `prismatique` using
      `conda`](#installing-prismatique-using-conda)
  - [Uninstalling `prismatique`](#uninstalling-prismatique)
- [Learning how to use `prismatique`](#learning-how-to-use-prismatique)



## Instructions for installing and uninstalling `prismatique`



### Installing `prismatique`

For all installation scenarios, first open up the appropriate command line
interface. On Unix-based systems, you could open e.g. a terminal. On Windows
systems you could open e.g. an Anaconda Prompt as an administrator.

GPU acceleration is available for `prismatique` installed on Linux and Windows
machines that have NVIDIA GPUs. You will need to make sure that you have a
NVIDIA driver installed with CUDA version 10.2.89 or greater.



#### Installing `prismatique` using `pip` and `conda` together

The easiest way to install `prismatique` involves using both the conda package
manager and `pip`. While it is possible to install `prismatique` without the use
of the conda package manager, it is more difficult. Because of this, we discuss
only the simplest installation procedure below.

Of course, to use the conda package manager, one must install either `anaconda3`
or `miniconda3`. For installation instructions for `anaconda3` click
[here](https://docs.anaconda.com/anaconda/install/index.html); for installation
instructions for `miniconda3` click
[here](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/macos.html).

First, open up the appropriate command line interface. On Unix-based systems,
you would open a terminal. On Windows systems you would open an Anaconda Prompt
as an administrator.

Next, you can optionally update your conda package manager by issuing the
following command:

    conda update conda

It is recommended that you install `prismatique` and its dependencies in a
virtual environment: click
[here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
for a discussion on the creation and management of conda virtual
environments. The remaining instructions assumes that you activate the conda
(virtual) environment in which you intend to install `prismatique` and its
dependencies.

The first dependency that we need to install is `pyprismatic`. GPU acceleration
is available for `pyprismatic` (and thus `prismatique`) if the following
conditions are met:

1. You are using a Linux or Windows machine that has NVIDIA GPUs.
2. A NVIDIA driver is installed with CUDA version 10.2.89 or greater.

If the above conditions have been met, and you would like to be able to use GPUs
with `prismatique`, run the following command:

    conda install -c conda-forge pyprismatic=2.\*=gpu\* cudatoolkit==<X>.<Y>.\*

where `<X>` and `<Y>` are the major and minor versions of CUDA installed on your
machine, e.g. CUDA version 10.2.89 has a major version of `10`, and a minor
version of `2`. Users can omit `cudatoolkit==<X>.<Y>.\*` if they do not require
a specific version of `cudatoolkit`, which should apply to most scenarios. For
CPU support only, run the following command instead:

    conda install -c conda-forge pyprismatic=2.\*=cpu\*

The easiest way to install the remaining dependencies, along with `prismatique`
is to use `pip` by running the following command:

    pip install prismatique

The above command will install the latest stable version of
`prismatique`. Another option is to use `conda`:

    conda install -c conda-forge prismatique

As yet another option, you can install the latest development version of
`prismatique` from the main branch of the [prismatique GitHub
repository](https://github.com/mrfitzpa/prismatique). To do so, one must first
clone the repository by running the following command:

    git clone https://github.com/mrfitzpa/prismatique.git

then subsequently change into the root of the cloned repository, and then run
the following command:

    pip install .

Note that you must include the period as well.

Optionally, for additional features in `prismatique`, one can install additional
dependencies upon installing `prismatique` via `pip`. To install a subset of
additional dependencies (along with the standard installation), run the
following command from the root of the repository:

    pip install .[<selector>]

where `<selector>` can be one of the following:

* `tests`: to install the dependencies necessary for running unit tests;
* `examples`: to install the dependencies necessary for executing files stored
  in `<root>/examples`, where `<root>` is the root of the repository;
* `docs`: to install the dependencies necessary for documentation generation;
* `all`: to install all of the above optional dependencies.

Alternatively, one can run:

    pip install prismatique[<selector>]

elsewhere in order to install the latest stable version of `prismatique`, along
with the subset of additional dependencies specified by `<selector>`.



### Uninstalling `prismatique`

If `prismatique` was installed using `pip`, then to uninstall, run the following
command:

    pip uninstall prismatique

If `prismatique` was installed using `conda`, then to uninstall, run the
following command:

    conda remove prismatique



## Learning how to use `prismatique`

For those new to the `prismatique` library, it is recommended that they take a
look at the [Examples](https://mrfitzpa.github.io/prismatique/examples.html)
page, which contain code examples that show how one can use the `prismatique`
library. While going through the examples, readers can consult the [prismatique
reference
guide](https://mrfitzpa.github.io/prismatique/_autosummary/prismatique.html) to
understand what each line of code is doing.