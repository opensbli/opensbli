#!/bin/sh

.PHONY: clean install test docs

install:
	@echo ">>> Installing..."
	python setup.py install

clean:
	@echo ">>> Cleaning..."
	python setup.py clean
	
test:
	@echo ">>> Running test suite..."
	py.test tests

docs:
	@echo ">>> Building documentation..."
	cd docs; make html; cd ..
