#!/bin/sh

EQUATIONS_FILE="equations.tex"

echo ">>> Cleaning up any existing generated files..."
rm -f alg.tex equations.tex *.cpp *.h *.f90
echo ">>> Done."

echo ">>> Generating code for equations..."
./generate.py equations
echo ">>> Done."

if [ -f $EQUATIONS_FILE ];
then
   echo ">>> LaTeX post-processing..."
   astyle -s2 --style=gnu --max-code-length=80 --convert-tabs auto_kernel.h
   # astyle -s2 --style=allman --indent=spaces=2  --max-code-length=80 --break-after-logical auto_kernel.h
   # astyle -s2 --style=allman --indent=spaces=2  --max-code-length=80 --break-after-logical auto_kernel.h
   astyle -s2 --style=gnu --max-code-length=80 OPSC_nssolver.cpp
   pdflatex $EQUATIONS_FILE
   echo ">>> Done."
else
   echo "!!! Warning: No LaTeX file exists. Skipping LaTeX post-processing..."
fi

echo ">>> Cleaning up build files..."
rm *.orig
rm *.pyc
rm *.aux
rm *.log
echo ">>> Done."
