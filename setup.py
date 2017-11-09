#!/usr/bin/env python

from distutils.core import setup

setup(name='OpenSBLI',
      version='2.b0',
      description='An automatic code generator which expands a set of equations written in Einstein notation, and writes out the finite difference code in a supported language.',
      author='',
      url='https://bitbucket.org/spjammy/opensbli',
      packages=['opensbli'],
      package_dir={'opensbli': 'opensbli'},
      scripts=[],
      classifiers=[
      'Development Status :: 2 - Beta',
      'Environment :: Console',
      'Intended Audience :: Science/Research',
      'Natural Language :: English',
      'Programming Language :: Python :: 2.7',
      'Topic :: Scientific/Engineering :: Physics',
      'Topic :: Software Development :: Code Generators',
      ]
      )
