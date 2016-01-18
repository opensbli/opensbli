rm -f alg.tex equations.tex *.cpp *.h *.f90
echo ">>> Expanding equations..."
python expansion.py equations
echo ">>> Done."
# sleep 2s
echo ">>> LaTeX post-processing..."
astyle -s2 --style=gnu   --max-code-length=80 --convert-tabs auto_kernel.h
# astyle -s2 --style=allman --indent=spaces=2  --max-code-length=80 --break-after-logical auto_kernel.h
# astyle -s2 --style=allman --indent=spaces=2  --max-code-length=80 --break-after-logical auto_kernel.h
astyle -s2 --style=gnu  --max-code-length=80  OPSC_nssolver.cpp
pdflatex equations.tex
echo ">>> Done."

echo ">>> Cleaning up build files..."
rm *.orig
rm *.pyc
rm *.aux
rm *.log
echo ">>> Done."
