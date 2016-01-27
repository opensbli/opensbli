#!/bin/sh

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
