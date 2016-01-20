#!/bin/sh

.PHONY: clean install test

install:
	@echo ">>> Installing..."
	python setup.py install

clean:
	@echo ">>> Cleaning..."
	python setup.py clean
	
test:
	@echo ">>> Running test suite..."
	py.test tests

