Introduction
============

Overview
--------

`AutoFD <https://bitbucket.org/spjammy/codegen>`_ is an automatic code generator which expands a set of equations written in Einstein notation, and writes out the corresponding finite difference code in OPSC format. This OPSC code can then be targetted with the `OPS library <http://www.oerc.ox.ac.uk/projects/ops>`_ towards specific hardware backends, such as MPI/OpenMP for CPUs, and CUDA/OpenCL for GPUs.

The AutoFD codebase is written in the Python (2.7+) language.

Licensing
---------

AutoFD is released under the `GNU General Public License <http://www.gnu.org/licenses/gpl-3.0.en.html>`_. See the file called ``LICENSE`` for more information.
