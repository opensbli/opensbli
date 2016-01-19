#!/usr/bin/env python

from distutils.core import setup

setup(name='autofd',
      version='0.1-dev',
      description='',
      author='Satya P. Jammy, Christian T. Jacobs',
      url='https://bitbucket.org/spjammy/codegen',
      packages=['autofd'],
      package_dir = {'autofd': 'autofd'},
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

