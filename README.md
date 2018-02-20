# OpenSBLI

OpenSBLI is an automatic code generator which expands a set of equations written in Einstein notation, and writes out the finite difference code in the OPSC language.

[![Build Status](https://travis-ci.org/opensbli/opensbli.svg?branch=master)](https://travis-ci.org/opensbli/opensbli)
[![Documentation Status](https://readthedocs.org/projects/opensbli/badge/?version=latest)](http://opensbli.readthedocs.io/en/latest/?badge=latest)

## Getting started

### Dependencies
First ensure that the following dependencies are satisfied:

* Python 2.7
* SymPy v1.0 (works only with SymPy v1.0)
* python-h5py
* pytest (for running the test suite)
* python-flake8 (for linting the code base)
* Sphinx (to build the documentation)

If you have pip installed, you can easily do this using `sudo pip install -r requirements.txt`.

* OPS (to target the generated OPSC code towards different backends should be installed manually), specifically revision 178ec4f7c1ccb1917a85b4248820cfecb912ac6f of the [OPS project's repository](https://github.com/gihanmudalige/OPS) or later.

### Installation

#### System-wide
OpenSBLI is a Python package and can be installed system-wide by running

```
sudo make install
```

from the OpenSBLI base directory (i.e. the same directory this `README` file is in).

#### Local builds
Alternatively, you can just add OpenSBLI to your `PYTHONPATH` environment variable using

```
export PYTHONPATH=$PYTHONPATH:/path/to/OpenSBLI/base/directory
```

## Documentation
The documentation for OpenSBLI can be built using Sphinx via the following command:

```
make docs
```

This will build the documentation in HTML format and can be opened in a Web browser. Alternatively, the OpenSBLI documentation can be viewed on [Read the Docs](http://opensbli.readthedocs.io/).

## Citing

If you use OpenSBLI, please consider citing the following paper:

* C. T. Jacobs, S. P. Jammy, N. D. Sandham (2017). **OpenSBLI: A framework for the automated derivation and parallel execution of finite difference solvers on a range of computer architectures**. *Journal of Computational Science*, 18:12-23, DOI: [10.1016/j.jocs.2016.11.001](http://dx.doi.org/10.1016/j.jocs.2016.11.001)

```
@Article{Jacobs_etal_2017,
  Title                    = {{OpenSBLI: A framework for the automated derivation and parallel execution of finite difference solvers on a range of computer architectures}},
  Author                   = {Jacobs, C. T. and Jammy, S. P. and Sandham, N. D.},
  Journal                  = {{Journal of Computational Science}},
  Year                     = {2017},
  Pages                    = {12--23},
  Volume                   = {18},
  Doi                      = {10.1016/j.jocs.2016.11.001}
}
```

## Contact
If you wish to report a bug with the software, please contact [Satya P. Jammy](mailto:S.P.Jammy@soton.ac.uk) or [Christian T. Jacobs](mailto:C.T.Jacobs@soton.ac.uk).

## Licence
OpenSBLI is released as an open-source project under the GNU General Public License. See the file called `LICENSE` for more information.
