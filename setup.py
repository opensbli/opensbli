#!/usr/bin/env python

#    AutoFD: An automatic code generator for solving differential equations.
#    Copyright (C) 2016 Satya P. Jammy, Christian T. Jacobs

#    This file is part of AutoFD.

#    AutoFD is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    AutoFD is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with AutoFD.  If not, see <http://www.gnu.org/licenses/>.

from distutils.core import setup

setup(name='AutoFD',
      version='0.1-dev',
      description='An automatic code generator which expands a set of equations written in Einstein notation, and writes out the finite difference code in the OPSC language.',
      author='Satya P. Jammy, Christian T. Jacobs',
      url='https://bitbucket.org/spjammy/codegen',
      packages=['autofd'],
      package_dir={'autofd': 'autofd'},
      scripts=["bin/autofd-clean"],
      classifiers=[
      'Development Status :: 2 - Pre-Alpha',
      'Environment :: Console',
      'Intended Audience :: Science/Research',
      'Natural Language :: English',
      'Programming Language :: Python :: 2.7',
      'Topic :: Scientific/Engineering :: Physics',
      'Topic :: Software Development :: Code Generators',
      ]
      )
