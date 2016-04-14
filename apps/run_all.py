#!/usr/bin/env python

""" Run all of the applications. This assumes that each application has its own run.py file. """

import os, os.path
import sys
import subprocess
import glob

apps = [name for name in os.listdir(".") if os.path.isdir(os.path.join(".", name))]

for app in apps:
    subprocess.call("cd %s; python run.py" % (app), shell=True)
