#!/bin/sh

EQUATIONS_FILE="equations.tex"
ALGORITHM_FILE="algorithm.tex"

# Clean up any old build files
echo ">>> Cleaning all build files..."
./autofd-clean --include-generated-code .
echo ">>> Done."

# Code generation step
echo ">>> Generating code for equations..."
./autofd-generate -c equations
echo ">>> Done."

# Make the generated code look pretty.
echo ">>> Stylising the generated source code files..."
astyle -s2 --style=gnu --max-code-length=80 --convert-tabs auto_kernel.h
# astyle -s2 --style=allman --indent=spaces=2  --max-code-length=80 --break-after-logical auto_kernel.h
# astyle -s2 --style=allman --indent=spaces=2  --max-code-length=80 --break-after-logical auto_kernel.h
astyle -s2 --style=gnu --max-code-length=80 OPSC_nssolver.cpp
echo ">>> Done."

# LaTeX post-processing step
echo ">>> LaTeX post-processing..."
if [ -f $EQUATIONS_FILE ];
then
   pdflatex $EQUATIONS_FILE
else
   echo "!!! Warning: No LaTeX file for the equations exists. Skipping..."
fi

if [ -f $ALGORITHM_FILE ];
then
   pdflatex $ALGORITHM_FILE
   echo ">>> Done."
else
   echo "!!! Warning: No LaTeX file for the algorithm exists. Skipping..."
fi
echo ">>> Done."

# Clean up the build files (except the newly-generated code).
echo ">>> Cleaning up build files (excluding the generated model code)..."
./autofd-clean .
echo ">>> Done."
