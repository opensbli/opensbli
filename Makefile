#!/bin/sh

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

.PHONY: clean install test lint docs

install:
	@echo ">>> Installing..."
	python setup.py install

clean:
	@echo ">>> Cleaning..."
	python setup.py clean

lint:
	@echo ">>> Linting..."
	flake8 autofd
	flake8 bin

test: lint
	@echo ">>> Running test suite..."
	py.test tests

docs:
	@echo ">>> Building documentation..."
	cd docs; make html; cd ..
