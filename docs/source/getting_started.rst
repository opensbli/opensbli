Getting Started
===============

Dependencies
------------

You should first ensure that all the core dependencies listed in the ``README.md`` file are satisfied. Many of these packages can be installed either via a package manager such as ``apt`` or via the `Python package manager (pip) <https://pypi.python.org/pypi/pip>`_.

Obtaining OPS
-------------

In order to target and compile the generated OPSC code, you will need to have OPS available. First, clone the `OPS GitHub <https://github.com/gihanmudalige/OPS>`_ repository using

.. code-block:: bash

    git clone https://github.com/gihanmudalige/OPS.git

You will then need to set up your OPS-related environment variables, listed below. Note that the values given here are system-dependent and may need to be adapted depending on where the MPI or HDF5 libraries are installed. Furthermore, it is assumed that the OPS GitHub repository has been cloned in your home (~) directory.

.. code-block:: bash

    export OPS_INSTALL_PATH=~/OPS/ops
    export OPS_COMPILER=gnu
    export MPI_INSTALL_PATH=/usr/
    export HDF5_INSTALL_PATH=/usr/

You can include these export commands in your ``~/.bashrc`` file to save typing them out each time you open up a new terminal.

Installing AutoFD
-----------------

You can install AutoFD using

.. code-block:: bash

    sudo make install

from within the base directory of AutoFD. Alternatively, particularly for developers of AutoFD, you can simply point your ``PYTHONPATH`` environment variable to the AutoFD base directory using, for example,

.. code-block:: bash

    export PYTHONPATH=$PYTHONPATH:~/autofd

Compiling the generated code
----------------------------
In the directory where the generated code resides, first copy the generated code into a new ``src`` directory:

.. code-block:: bash

    mkdir src
    cp OPSC_nssolver.cpp src/
    cp auto_kernel.h src/
    cd src

In the copied source, insert the simulation parameters in ``OPSC_nssolver.cpp``, and then translate this code with the OPS translator using:

.. code-block:: bash

    python ~/OPS/translator/python/c/ops.py OPSC_nssolver.cpp

This will give you the code targetted towards different backends, e.g. CUDA, MPI, OpenMP, etc.

Finally, copy across the ``Makefile`` from one of the OPS ``apps``, and modify it so that it will compile the source for your simulation setup.
