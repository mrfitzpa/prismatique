.. _examples_prerequisites_for_running_examples_sec:

Prerequisites for running example scripts or Jupyter notebooks
==============================================================

Prior to running any scripts or Jupyter notebooks in the directory
``<root>/examples``, where ``<root>`` is the root of the ``prismatique``
repository, a set of Python libraries need to be installed in the Python
environment within which any such scripts or Jupyter notebooks are to be
executed.

The Python libraries that need to be installed in said Python environment are::

  pyprismatic>=2.0
  jupyter
  ipympl
  prismatique

The easiest way to install these libraries is within a ``conda`` virtual
environment. It is recommended that ``pyprismatic`` be installed separately,
prior to installing any of the other Python libraries listed above. GPU
acceleration is available for ``pyprismatic`` (and thus ``prismatique``) if the
following conditions are met:

1. You are using a Linux or Windows machine that has NVIDIA GPUs.
2. A NVIDIA driver is installed with CUDA version 10.2.89 or greater.

If the above conditions have been met, and you would like to be able to use GPUs
with ``prismatique``, run the following command::

  conda install -c conda-forge pyprismatic=2.*=gpu* cudatoolkit==<X>.<Y>.*

where ``<X>`` and ``<Y>`` are the major and minor versions of CUDA installed on
your machine, e.g. CUDA version 10.2.89 has a major version of ``10``, and a
minor version of ``2``. Users can omit ``cudatoolkit==<X>.<Y>.*`` if they do not
require a specific version of ``cudatoolkit``, which should apply to most
scenarios. For CPU support only, run the following command instead::

  conda install -c conda-forge pyprismatic=2.*=cpu*

Once ``pyprismatic`` is installed, you can install the remaining libraries via
``pip`` by running the following command::

  pip install prismatique[examples]

or alternatively, via ``conda`` by running the following command::

  conda install -c conda-forge prismatique ipympl jupyter
