# AutoFD

AutoFD is an automatic code generator which expands a set of equations written in Einstein notation, and writes out the finite difference code in either OPSC or Fortran (serial).

## Getting started

### Dependencies
First ensure that the following dependencies are satisfied:

* Python 2.7
* astyle
* Sympy >= 0.7.6.1
* python-h5py
* pytest (for running the test suite)
* python-flake8 (for linting the code base)
* Sphinx (to build the documentation)

### Installation

#### System-wide
AutoFD is a Python package and can be installed system-wide by running

```
sudo make install
```

from the AutoFD base directory (i.e. the same directory this README file is in).

#### Local builds
Alternatively, you can just add AutoFD to your `PYTHONPATH` environment variable using

```
export PYTHONPATH=$PYTHONPATH:/path/to/autofd/base/directory
```

## Documentation
The documentation for AutoFD can be built using Sphinx via the following command:

```
make docs
```

This will build the documentation in HTML format and can be opened in a Web browser.

## Contact
If you wish to report a bug with the software, please contact [Satya P Jammy](mailto:S.P.Jammy@soton.ac.uk) or [Christian T. Jacobs](mailto:C.T.Jacobs@soton.ac.uk).

## Licence
AutoFD is released under the GNU General Public License. See the file called `LICENSE` for more information.

