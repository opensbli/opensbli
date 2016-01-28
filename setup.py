#!/usr/bin/env python

from distutils.core import setup

setup(name='AutoFD',
      version='0.1-dev',
      description='An automatic code generator which expands a set of equations written in Einstein notation, and writes out the finite difference code in either OPSC or Fortran.',
      author='Satya P. Jammy, Christian T. Jacobs',
      url='https://bitbucket.org/spjammy/codegen',
      packages=['autofd'],
      package_dir={'autofd': 'autofd'},
      scripts=["bin/autofd-generate", "bin/autofd-clean"],
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
