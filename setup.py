#!/usr/bin/env python

#    OpenSBLI: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs, Neil D. Sandham

#    This file is part of OpenSBLI.

#    OpenSBLI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    OpenSBLI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with OpenSBLI.  If not, see <http://www.gnu.org/licenses/>.

from distutils.core import setup

setup(name='OpenSBLI',
      version='1.0',
      description='An automatic code generator which expands a set of equations written in Einstein notation, and writes out the finite difference code in a supported language.',
      author='Satya P. Jammy, Christian T. Jacobs, Neil D. Sandham',
      url='https://bitbucket.org/spjammy/opensbli',
      packages=['opensbli'],
      package_dir={'opensbli': 'opensbli'},
      scripts=[],
      classifiers=[
      'Environment :: Console',
      'Intended Audience :: Science/Research',
      'Natural Language :: English',
      'Programming Language :: Python :: 2.7',
      'Topic :: Scientific/Engineering :: Physics',
      'Topic :: Software Development :: Code Generators',
      ]
      )
