# OpenSBLI

OpenSBLI is an automatic code generator which expands a set of equations written in Einstein notation, and writes out the finite difference code in the OPSC language.

## Getting started

### Dependencies
First ensure that the following dependencies are satisfied:
#### Core code-generation:
The following dependencies are required for generating a code and running a simulation:

* Python 2.7
* Sympy >= 1.0
* Numpy
* Scipy 0.19.1
* OPS (to target the generated OPSC code towards different backends) [OPS project's repository](https://github.com/gihanmudalige/OPS).

#### Testing and documentation:

* pytest (for running the test suite)
* python-flake8 (for linting the code base)
* Sphinx (to build the documentation)

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

<!-- #### System-wide
OpenSBLI is a Python package and can be installed system-wide by running

```
sudo make install
```

from the OpenSBLI base directory (i.e. the same directory this README file is in).



## Documentation
The documentation for OpenSBLI can be built using Sphinx via the following command:

```
make docs
```

This will build the documentation in HTML format and can be opened in a Web browser.

## Contact
If you wish to report a bug with the software, please contact [Satya P. Jammy](mailto:S.P.Jammy@soton.ac.uk) or [Christian T. Jacobs](mailto:C.T.Jacobs@soton.ac.uk).

## Licence
OpenSBLI is released under the GNU General Public License. See the file called `LICENSE` for more information. -->