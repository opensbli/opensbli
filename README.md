# OpenSBLI
OpenSBLI is an open-source code-generation system for compressible fluid dynamics (CFD) on heterogeneous computing architectures. Written in Python, OpenSBLI uses explicit high-order finite-difference schemes on structured curvilinear meshes. Shock-capturing is performed by a choice of high-order Weighted Essentially Non-Oscillatory (WENO) or Targeted Essentially Non-Oscillatory (TENO) schemes. OpenSBLI generates a complete CFD solver in the Oxford Parallel Structured (OPS) domain specific language. The OPS library is embedded in C code, enabling massively-parallel execution of the code on a variety of high-performance-computing architectures, including GPUs. 

## How to cite this work
The current reference for OpenSBLI is: 

D.J. Lusher, S.P. Jammy, N.D. Sandham. OpenSBLI: Automated code-generation for heterogeneous computing architectures applied to compressible fluid dynamics on structured grids. Computer Physics Communications Vol. 267, October 2021, 108063.
 
```bash
@article{OpenSBLI_LJS2021,
title = {{OpenSBLI: Automated code-generation for heterogeneous computing architectures applied to compressible fluid dynamics on structured grids}},
journal = {Computer Physics Communications},
volume = {267},
pages = {108063},
year = {2021},
issn = {0010-4655},
doi = {https://doi.org/10.1016/j.cpc.2021.108063},
author = {David J. Lusher and Satya P. Jammy and Neil D. Sandham},
}
```

## Getting started

### Dependencies
First ensure that the following dependencies are satisfied:
#### Core code-generation:
The following dependencies are required for generating a code and running a simulation:

* Python 2.7
* Sympy == 1.1
* Numpy
* Scipy 0.19.1
* OPS (to target the generated OPSC code towards different backends) [OPS project's repository](https://github.com/gihanmudalige/OPS).

#### Testing and documentation:

* pytest (for running the test suite)
* python-flake8 (for linting the code base)
* Sphinx (to build the documentation)

#### Note on previous version

* The version 1.0 of OpenSBLI can be downloaded from [release][https://github.com/opensbli/opensbli/releases]
* No enhancements of version 1.0 are planned

#### Postprocessing:

* Matplotlib for plot scripts
* python-h5py

## Installation

### Development branch

Add OpenSBLI to your `PYTHONPATH` environment variable using

```
export PYTHONPATH=$PYTHONPATH:/path/to/OpenSBLI/base/directory
```

## Contact
If you wish to report a bug with the software, please contact [Satya P. Jammy](mailto:S.P.Jammy@soton.ac.uk) or [David J. Lusher](mailto:D.Lusher@soton.ac.uk)
