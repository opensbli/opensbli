#!/bin/sh

EQUATIONS_FILE="equations.tex"
ALGORITHM_FILE="alg.tex"

echo ">>> Cleaning all build files..."
./autofd-clean --include-generated-code .
echo ">>> Done."

echo ">>> Generating code for equations..."
./autofd-generate equations
echo ">>> Done."

if [ -f $EQUATIONS_FILE ];
then
   echo ">>> LaTeX post-processing..."
   astyle -s2 --style=gnu --max-code-length=80 --convert-tabs auto_kernel.h
   # astyle -s2 --style=allman --indent=spaces=2  --max-code-length=80 --break-after-logical auto_kernel.h
   # astyle -s2 --style=allman --indent=spaces=2  --max-code-length=80 --break-after-logical auto_kernel.h
   astyle -s2 --style=gnu --max-code-length=80 OPSC_nssolver.cpp
   pdflatex $EQUATIONS_FILE
   pdflatex $ALGORITHM_FILE
   echo ">>> Done."
else
   echo "!!! Warning: No LaTeX file exists. Skipping LaTeX post-processing..."
fi

echo ">>> Cleaning up build files (excluding the generated model code)..."
./autofd-clean .
echo ">>> Done."
