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
      data_files=[]
     )

