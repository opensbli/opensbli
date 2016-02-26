Introduction
============

Overview
--------

`OpenSBLI <https://bitbucket.org/spjammy/autofd>`_ is an automatic code generator that expands a set of equations written in Einstein notation, and automatically generates code (in the OPSC language) which performs the finite difference approximation to obtain a solution. This OPSC code can then be targetted with the `OPS library <http://www.oerc.ox.ac.uk/projects/ops>`_ towards specific hardware backends, such as MPI/OpenMP for execution on CPUs, and CUDA/OpenCL for execution on GPUs.

The main focus of OpenSBLI is on the solution of the compressible Navier-Stokes equations with application to shock-boundary layer interactions (SBLI). However, in principle, any set of equations that can be written in Einstein notation may be solved using the code generation framework. This highlights one of the main advantages of such a high-level, abstract approach to computational model development.

From an implementation perspective, the OpenSBLI codebase is written in the Python (2.7.x) language and depends on the SymPy library to process the Einstein formulation of the equations. The code generator will then write out the model code which performs the finite difference approximations in any of the supported languages (currently only OPSC, although the structure of the codebase is such that other languages can be integrated with minimal effort).

The development of OpenSBLI was supported by the `UK Turbulence Consortium <http://www.turbulence.ac.uk>`_ and the `ExaFLOW project <http://exaflow-project.eu/>`_.

Licensing
---------

OpenSBLI is released under the `GNU General Public License <http://www.gnu.org/licenses/gpl-3.0.en.html>`_. See the file called ``LICENSE`` for more information.

Support
-------

The preferred method of reporting bugs and issues with OpenSBLI is to submit an issue via the repository's issue tracker. One can also email the authors `Satya P Jammy <mailto:S.P.Jammy@soton.ac.uk>`_ and `Christian T. Jacobs <mailto:C.T.Jacobs@soton.ac.uk>`_ directly.

