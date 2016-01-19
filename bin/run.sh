#!/bin/sh

EQUATIONS_FILE="equations.tex"
BUILD_DIR=`printenv AUTOFD_BUILD_DIR`
CURRENT_DIR=`pwd`

echo ">>> Cleaning all build files..."
./autofd-clean --include-generated-code
echo ">>> Done."

echo ">>> Generating code for equations..."
./generate equations
echo ">>> Done."

if [ -f $EQUATIONS_FILE ];
then
   echo ">>> LaTeX post-processing..."
   cd $AUTOFD_BUILD_DIR
   astyle -s2 --style=gnu --max-code-length=80 --convert-tabs auto_kernel.h
   # astyle -s2 --style=allman --indent=spaces=2  --max-code-length=80 --break-after-logical auto_kernel.h
   # astyle -s2 --style=allman --indent=spaces=2  --max-code-length=80 --break-after-logical auto_kernel.h
   astyle -s2 --style=gnu --max-code-length=80 OPSC_nssolver.cpp
   pdflatex $EQUATIONS_FILE
   echo ">>> Done."
else
   echo "!!! Warning: No LaTeX file exists. Skipping LaTeX post-processing..."
fi

echo ">>> Cleaning up build files (excluding the generated model code)..."
./autofd-clean
echo ">>> Done."
