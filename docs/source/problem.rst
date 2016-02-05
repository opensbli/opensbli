Defining and Running a Problem
==============================
Essentially, AutoFD comprises the following classes and modules (emboldened below), which define the abstraction employed:

* A **Problem** defines the physical problem's dimension, the equations that must be solved, and any accompanying formulas, constants, etc.
* This Problem comprises many **Equations**, and an **Algorithm**.
* After expanding the Equations, the Algorithm is applied to the expanded equations and formulas to prepare a computational **System**.
* With this System, and the help of ``codegen_utils.py``, the **OPSC** code can be generated.
* All LaTeX writing is handled by the **LatexWriter** class.

AutoFD will expect an ``equations`` file, and an ``algorithm`` file containing all the problem-specific settings and configurations (the governing equations, any constitutive equations for e.g. temperature-dependent viscosity, what time-stepping scheme is to be used, the boundary conditions, etc.). There are several examples of these files provided in the applications (``apps``) directory of the AutoFD package.

Generating code
---------------
Once defined, AutoFD can then read in these files and generate the code for this particular problem. To proceed with the code generation step, use the ``bin/autofd-generate`` script:

.. code-block:: bash

    python bin/autofd-generate /path/to/directory/containing/problem/files

For the OPSC format, AutoFD's code generator will create two files: ``auto_kernel.h``, ``OPSC_nssolver.cpp``, and the translated file ``OPSC_nssolver_ops.cpp``.

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
