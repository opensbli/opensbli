Introduction
============

Overview
--------

`OpenSBLI <https://bitbucket.org/spjammy/autofd>`_ is an automatic code generator which expands a set of equations written in Einstein notation, and writes out the corresponding finite difference code in the OPSC language. This OPSC code can then be targetted with the `OPS library <http://www.oerc.ox.ac.uk/projects/ops>`_ towards specific hardware backends, such as MPI/OpenMP for CPUs, and CUDA/OpenCL for GPUs.

The OpenSBLI codebase is written in the Python (2.7.x) language and depends on the SymPy library to process the Einstein formulation of the equations.

Licensing
---------

OpenSBLI is released under the `GNU General Public License <http://www.gnu.org/licenses/gpl-3.0.en.html>`_. See the file called ``LICENSE`` for more information.

Support
-------

The preferred method of reporting bugs and issues with OpenSBLI is to submit an issue via the repository's issue tracker. One can also email the authors `Satya P Jammy <mailto:S.P.Jammy@soton.ac.uk>`_ and `Christian T. Jacobs <C.T.Jacobs@soton.ac.uk>`_ directly.
