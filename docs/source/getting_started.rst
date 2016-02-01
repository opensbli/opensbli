Getting Started
===============

Environment variables
---------------------
In order to compile the generated code, you will need to set up your OPS installation using, for example,

```
export OPS_INSTALL_PATH=~/OPS/ops
export OPS_COMPILER=gnu
```

Compiling the generated code
----------------------------
In the directory where the generated code resides, first copy the generated code into a new `src` directory:

```
mkdir src
cp OPSC_nssolver.cpp src/
cp auto_kernel.h src/
cd src
```

In the copied source, insert the simulation parameters in `OPSC_nssolver.cpp`, and then translate this code with the OPS translator using:

```
python ~/OPS/translator/python/c/ops.py OPSC_nssolver.cpp
```

This will give you the code targetted towards different backends, e.g. CUDA, MPI, OpenMP, etc.

Finally, copy across the `Makefile` from one of the OPS `apps`, and modify it so that it will compile the source for your simulation setup.
