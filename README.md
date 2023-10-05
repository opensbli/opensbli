# OpenSBLI
OpenSBLI is an open-source code-generation system for compressible fluid dynamics (CFD) on heterogeneous computing architectures. Written in Python, OpenSBLI uses explicit high-order finite-difference schemes on structured curvilinear meshes. Shock-capturing is performed by a choice of high-order Weighted Essentially Non-Oscillatory (WENO) or Targeted Essentially Non-Oscillatory (TENO) schemes. OpenSBLI generates a complete CFD solver in the Oxford Parallel Structured (OPS) domain specific language. The OPS library is embedded in C code, enabling massively-parallel execution of the code on a variety of high-performance-computing architectures, including GPUs. 

## How to cite this work
The current reference for OpenSBLI is: 

D.J. Lusher, S.P. Jammy, N.D. Sandham. *OpenSBLI: Automated code-generation for heterogeneous computing architectures applied to compressible fluid dynamics on structured grids*. **Computer Physics Communications Vol. 267, 108063 (2021)**.
 
```
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

## Recent applications of OpenSBLI
1. DJ Lusher, M Zauner, A Sansica, A Hashimoto. *Automatic Code-Generation to Enable High-Fidelity Simulations of Multi-Block Airfoils on GPUs*. **AIAA SciTech Forum, 1222 (2023)**.
2. A Gillespie, ND Sandham. *Numerical study of the effect of sidewalls on shock train behaviour*. **Flow 3, E12 (2023)**.
3. DJ Lusher, GN Coleman. *Numerical Study of Compressible Wall-Bounded Turbulence–the Effect of Thermal Wall Conditions on the Turbulent Prandtl Number in the Low-Supersonic Regime*. **International Journal of Computational Fluid Dynamics 36 (9), 797-815 (2022)**.
4. DJ Lusher, GN Coleman. *Numerical Study of the Turbulent Prandtl Number in Supersonic Plane-Channel Flow – the Effect of Thermal Boundary Conditions*. **NASA Technical Memorandum 10483 (NASA/TM–20220010483) (2022)**.
5. A Gillespie, ND Sandham. *Shock train response to high-frequency backpressure forcing*. **AIAA Journal 60 (6), 3736-3748 (2022)**.
6. A Hamzehloo, DJ Lusher, S Laizet, ND Sandham. *Direct numerical simulation of compressible turbulence in a counter-flow channel configuration*. **Physical Review Fluids 6 (9), 094603 (2021)**.
7. DJ Lusher, ND Sandham. *Assessment of low-dissipative shock-capturing schemes for the compressible Taylor–Green vortex*. **AIAA Journal 59 (2), 533-545 (2021)**.
8. A Hamzehloo, DJ Lusher, S Laizet, ND Sandham. *On the performance of WENO/TENO schemes to resolve turbulence in DNS/LES of high‐speed compressible flows*. **International Journal for Numerical Methods in Fluids 93 (1), 176-196 (2021)**.
9. DJ Lusher, ND Sandham. *Shock-wave/boundary-layer interactions in transitional rectangular duct flows*. **Flow, Turbulence and Combustion 105, 649-670 (2020)**.
10. DJ Lusher, ND Sandham. *The effect of flow confinement on laminar shock-wave/boundary-layer interactions*. **Journal of Fluid Mechanics 897, A18 (2020)**.
11. DJ Lusher, SP Jammy, ND Sandham. *Shock-wave/boundary-layer interactions in the automatic source-code generation framework OpenSBLI*. **Computers & Fluids 173, 17-21 (2018)**.

## Contact
If you wish to report a bug with the software, please contact [Satya P. Jammy](mailto:S.P.Jammy@soton.ac.uk) or [David J. Lusher](mailto:lusher.david@jaxa.jp)
