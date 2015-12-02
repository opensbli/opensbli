rm alg.tex equations.tex *.cpp *.h *.f90
./expansion.py equations
echo "Expanded equations waiting for 10 seconds"
# sleep 2s
echo "doing post stuff"
astyle -s2 --max-code-length=60 auto_kernel.h
astyle -s2 --max-code-length=60 OPSC_nssolver.cpp
pdflatex alg.tex
# pdflatex equations.tex
rm *.orig
rm *.pyc
rm *.aux
rm *.log