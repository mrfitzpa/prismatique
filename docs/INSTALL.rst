.. _installation_instructions_sec:

Installation instructions
=========================

This page covers how to install ``prismatique`` on both non--Digital Research
Alliance of Canada (DRAC) machines and DRAC machines.

In both cases, the installation instructions assume that ``git`` has been
installed and that git commands can be issued from the appropriate command line
interface. On Unix-based systems (e.g. Linux and Mac) this is the terminal, and
on Windows systems this is the Anaconda Prompt ran as an Administrator.

GPU acceleration is available for ``prismatique`` installed on machines that
have NVIDIA GPUs. You will need to make sure that you have a NVIDIA driver
installed with CUDA version 10.2.89 or greater. DRAC servers also has NVIDIA
GPUs available, see their website for more details.

Typical installation times range from 3-15 minutes.

Installation on non--Digital Research Alliance of Canada machines
-----------------------------------------------------------------

**Note**: these instructions have been tested on ``Ubuntu 20.04`` and ``Windows
10``.

The easiest way to install ``prismatique`` involves using both the conda package
manager and ``pip``. While it is possible to install ``prismatique`` without the
use of the conda package manager, it is more difficult. Because of this, we
discuss only the simplest installation procedure below.

Of course, to use the conda package manager, one must install either
``anaconda3`` or ``miniconda3``. For installation instructions for ``anaconda3``
click `here <https://docs.anaconda.com/anaconda/install/index.html>`__; for
installation instructions for ``miniconda3`` click `here
<https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/macos.html>`__.

First, open up the appropriate command line interface. On Unix-based systems,
you would open a terminal. On Windows systems you would open an Anaconda Prompt
as an administrator.

Next, you can optionally update your conda package manager by issuing the
following command (WARNING: IF YOU ARE USING A MACHINE IN WHICH YOU ARE NOT
SUPPOSE TO UPDATE THE CONDA PACKAGE, DO NOT ISSUE THE FOLLOWING COMMAND)::

  conda update conda

Next, assuming that you have downloaded/cloned the ``prismatique`` git
repository, change into the root directory of said repository.

To install the minimal requirements for ``prismatique`` and ``python3.<Y>`` ---
where ``<Y>`` is the minor version of Python of interest (e.g.  ``7`` or ``9``)
--- into a new virtual environment with the name ``<name>`` (e.g.
``prismatique`` or ``stem_sim``), issue the following commands (may take several
minutes)::

  python setup_env.py --python=3.<Y> --name=<name>
  conda activate <name>
  
This will create a new virtual environment named ``<name>`` with the
aforementioned minimal requirements installed, and activate the new virtual
environment. The minimal requirements include ``numpy``, ``scipy``, ``sympy``,
``scikit-image``, ``scikit-learn``, ``hyperspy``, ``h5py``, ``matplotlib``,
``numba``, ``panda``, ``pyFAI``, ``pyopencl``, ``pyprismatic``, ``pytest``,
``emconstants``, ``czekitout``, ``fancytypes``, ``empix``, and ``embeam``. The
``--python`` argument can be omitted in which case the latest possible version
of python will be installed. The ``--name`` argument can also be omitted, in
which case the virtual environment name will default to ``prismatique``.

The above installation script will automatically determine whether your machine
supports GPU acceleration for ``prismatique``: if GPU acceleration is supported
then the GPU version of ``prismatique`` will be automatically installed,
otherwise the CPU-only version of ``prismatique`` will be automatically
installed.

Once you have installed the required packages, you can install ``prismatique``
by issuing the following command within your virtual environment::

  pip install .

Note that you must include the period as well.

Optionally, for additional features in ``prismatique``, one can install
additional dependencies upon installing ``prismatique``. To install a subset of
additional dependencies, run the following command from the root of the
repository::

  pip install .[<selector>]

where ``<selector>`` can be one of the following:

* ``doc``: to install the dependencies necessary for documentation generation;
* ``examples``: to install the dependencies necessary for running any example
  notebooks;
* ``all``: to install all additional dependencies.

Installation on Digital Research Alliance of Canada
---------------------------------------------------

First, SSH into the DRAC server that you intend to run ``prismatique`` jobs
on. Next, assuming that you have downloaded/cloned the ``prismatique`` git
repository, change into the root directory of said repository.

To install the minimal requirements for ``prismatique`` into a new virtual
environment with the name ``<name>`` (e.g. ``prismatique`` or ``stem_sim``),
issue the following command (may take several minutes)::

  source setup_drac_env.sh <name>

This will first load a subset of modules pre-installed on DRAC that are required
to install ``prismatique``, then create a new virtual environment named
``<name>`` with the remaining minimal requirements to install ``prismatique``,
and then activate said virtual environment. Note that the above script will
create a virtual environment with ``python3.9``. Moreover, the above script will
only try to install the GPU-version of ``prismatique``.

Once you have installed the required packages, you can install ``prismatique``
by issuing the following command within your virtual environment::

  pip install .

Note that you must include the period as well.

Every time you SSH into DRAC, you must load the same set of modules
pre-installed on DRAC that are required to run ``prismatique``. This can be done
by running the following command from the root of the ``prismatique``
repository::

  source load_drac_modules.sh

Then you must activate your virtual environment by running the following command
from anywhere::

  source ~/<name>/bin/activate

where ``<name>`` is the name of the virtual environment.

In the case that you are done running ``prismatique`` jobs and you would like to
e.g. work in another environment that requires a different set of modules
pre-installed on DRAC, you can unload the set you previously loaded for
``prismatique`` by running the following command from the root of the
``prismatique`` repository::

  source unload_drac_modules.sh

Update prismatique
------------------

If you, or someone else has made changes to this library, you can reinstall it
by issuing the following command from the root of the repository (within your
virtual environment)::
  
  pip install .

or the command::

  pip install .[<selector>]

where ``<selector>`` was described in the first section.

Uninstall prismatique
---------------------

To uninstall ``prismatique``, run the following command from the root of the
repository (within your virtual environment)::

  pip uninstall prismatique

Exploring examples of using prismatique
---------------------------------------

Examples of using ``prismatique`` can be found in a set of scripts and notebooks
in the directory ``<root>/examples``, where ``<root>`` is the root of the
repository. The dependencies required for running these example scripts and
notebooks can be installed by running the following command from the root of the
repository (within your virtual environment)::

  pip install .[examples]

or the command::

  pip install .[all]

Note that the latter command will install all extra dependencies of
``prismatique``.

Since the repository tracks the notebooks under their original basenames, we
recommend that you copy whatever original notebook of interest and rename it to
whatever other basename before executing any cells. This way you can explore any
notebook by executing and modifying cells without changing the originals, which
are being tracked by git.

Generating documention files on non--Digital Research Alliance of Canada machines
---------------------------------------------------------------------------------

To generate documentation in html format from source files, you will also need
to install several other packages. This can be done by running the following
command from the root of the repository::

  pip install .[doc]

or the command::

  pip install .[all]

Note that the latter command will install all extra dependencies of
``prismatique``.

Assuming that you are in the root of the repository, that you have installed all
the prerequisite packages, and that ``prismatique`` has been installed, you can
generate the ``prismatique`` documentation html files by issuing the following
commands within your virtual environment::

  cd docs
  make html

This will generate a set of html files in ``./_build/html`` containing the
documentation of ``prismatique``. You can then open any of the files using your
favorite web browser.

If ``prismatique`` has been updated, the documentation has most likely changed
as well. To update the documentation simply run::

  make html

again to generate the new documentation.
