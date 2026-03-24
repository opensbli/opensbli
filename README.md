# OpenSBLI
OpenSBLI is an open-source code-generation system for compressible fluid dynamics (CFD) on heterogeneous computing architectures, including GPUs. Written in Python, OpenSBLI uses explicit high-order finite-difference schemes on structured curvilinear multi-block meshes. Shock-capturing is performed by a choice of high-order Weighted Essentially Non-Oscillatory (WENO) or Targeted Essentially Non-Oscillatory (TENO) schemes. OpenSBLI generates a complete CFD solver in the Oxford Parallel Structured (OPS) domain specific language. The OPS library is embedded in C code, enabling massively-parallel execution of the code on a variety of high-performance-computing architectures, including GPUs. 

## How to cite this work
The latest references for the OpenSBLI solver are: 

D.J. Lusher, A. Sansica, N.D. Sandham, J. Meng, B. Siklósi, A. Hashimoto. *OpenSBLI v3.0: High-fidelity multi-block transonic aerofoil CFD simulations using domain specific languages on GPUs*. **Computer Physics Communications Vol. 307, 109406 (2025)**.

D.J. Lusher, S.P. Jammy, N.D. Sandham. *OpenSBLI: Automated code-generation for heterogeneous computing architectures applied to compressible fluid dynamics on structured grids*. **Computer Physics Communications Vol. 267, 108063 (2021)**.
 
```
@article{LSSMSH_OpenSBLIv3_CPC2025,
title = {{OpenSBLI v3.0: High-fidelity multi-block transonic aerofoil CFD simulations using domain specific languages on GPUs}},
journal = {Computer Physics Communications},
volume = {307},
pages = {109406},
year = {2025},
issn = {0010-4655},
doi = {https://doi.org/10.1016/j.cpc.2024.109406},
author = {David J. Lusher and Andrea Sansica and Neil D. Sandham and Jianping Meng and Bálint Siklósi and Atsushi Hashimoto},
}

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

* Python 3.7
* Sympy == 1.1
* Numpy
* Scipy
* OPS (to target the generated OPSC code towards different backends) [OPS project's repository](https://github.com/gihanmudalige/OPS).

#### Testing and documentation:

* pytest (for running the test suite)
* python-flake8 (for linting the code base)

#### Postprocessing:

* Matplotlib for plot scripts
* python-h5py

## Installation

For a traditional makefile-based installation see the instructions in the user manual (under the docs directory of this branch).

For a simplified installation script (suitable for novice users) that avoids detailed path specifications and uses Cmake for installation and OPS translation, see the readme_script_based_installation.md file in this branch.

<!--
### Development branch

 Add OpenSBLI to your `PYTHONPATH` environment variable using 

```
export PYTHONPATH=$PYTHONPATH:/path/to/OpenSBLI/base/directory
```
-->

## Recent applications of OpenSBLI
1. A Musawi, ND Sandham. Effects of thermal non-equilibrium during transition to turbulence of a high enthalpy free shear layer. **Physics of Fluids 38 (1) (2026)**.
2. B Siklósi, PK Sharma, DJ Lusher, IZ Reguly, ND Sandham. Reduced and mixed precision turbulent flow simulations using explicit finite difference schemes. **Future Generation Computer Systems, 175, 108111 (2026)**.
3. A Musawi, ND Sandham. Efficient viscosity and thermal conductivity formulation for scale-resolved hypersonic flow simulations. **AIAA Journal 63 (12), 5123-5135 (2025)**.
4. M Mauriello, PK Sharma, L Larchevêque, N Sandham. Role of nonlinearities induced by deterministic forcing in the low-frequency dynamics of transitional shock wave/boundary layer interaction. **Journal of Fluid Mechanics 1016, A6 (2025)**.
5. DJ Lusher, A Sansica, A Hashimoto. *Implicit large eddy simulations of three-dimensional turbulent transonic buffet on wide-span infinite wings*. **Journal of Fluid Mechanics (2025)**.
6. JB Chapelier, DJ Lusher, et al. *Comparison of high-order numerical methodologies for the simulation of the supersonic Taylor–Green vortex flow*. **Physics of Fluids 36, 055146 (2024).**
7. DJ Lusher, A Sansica, A Hashimoto. *Effect of Tripping and Domain Width on Transonic Buffet on Periodic NASA-CRM Airfoils*. **AIAA Journal 62 (11), 4411-4430 (2024)**.
8. A Hamzehloo, DJ Lusher, ND Sandham. *Direct numerical simulations and spectral proper orthogonal decomposition analysis of shocklet-containing turbulent channel counter-flows*. **International Journal of Heat and Fluid Flow, 104, 109229 (2023)**.
9. A Gillespie, ND Sandham. *Numerical study of the effect of sidewalls on shock train behaviour*. **Flow 3, E12 (2023)**.
10. DJ Lusher, GN Coleman. *Numerical Study of Compressible Wall-Bounded Turbulence–the Effect of Thermal Wall Conditions on the Turbulent Prandtl Number in the Low-Supersonic Regime*. **International Journal of Computational Fluid Dynamics 36 (9), 797-815 (2022)**.
11. A Gillespie, ND Sandham. *Shock train response to high-frequency backpressure forcing*. **AIAA Journal 60 (6), 3736-3748 (2022)**.
12. A Hamzehloo, DJ Lusher, S Laizet, ND Sandham. *Direct numerical simulation of compressible turbulence in a counter-flow channel configuration*. **Physical Review Fluids 6 (9), 094603 (2021)**.
13. DJ Lusher, ND Sandham. *Assessment of low-dissipative shock-capturing schemes for the compressible Taylor–Green vortex*. **AIAA Journal 59 (2), 533-545 (2021)**.
14. A Hamzehloo, DJ Lusher, S Laizet, ND Sandham. *On the performance of WENO/TENO schemes to resolve turbulence in DNS/LES of high‐speed compressible flows*. **International Journal for Numerical Methods in Fluids 93 (1), 176-196 (2021)**.
15. DJ Lusher, ND Sandham. *Shock-wave/boundary-layer interactions in transitional rectangular duct flows*. **Flow, Turbulence and Combustion 105, 649-670 (2020)**.
16. DJ Lusher, ND Sandham. *The effect of flow confinement on laminar shock-wave/boundary-layer interactions*. **Journal of Fluid Mechanics 897, A18 (2020)**.
17. DJ Lusher, SP Jammy, ND Sandham. *Shock-wave/boundary-layer interactions in the automatic source-code generation framework OpenSBLI*. **Computers & Fluids 173, 17-21 (2018)**.

## Contact
If you wish to report a bug with the software, please contact [Satya P. Jammy](mailto:satya.jammy@ddn.upes.ac.in) or [David J. Lusher](mailto:lusher.david@jaxa.jp)
