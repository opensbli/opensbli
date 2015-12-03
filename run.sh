rm alg.tex equations.tex *.cpp *.h *.f90
./expansion.py equations
echo "Expanded equations waiting for 10 seconds"
# sleep 2s
echo "doing post stuff"
astyle -s2 --style=gnu   --max-code-length=80 --convert-tabs auto_kernel.h
# astyle -s2 --style=allman --indent=spaces=2  --max-code-length=80 --break-after-logical auto_kernel.h
# astyle -s2 --style=allman --indent=spaces=2  --max-code-length=80 --break-after-logical auto_kernel.h
astyle -s2 --style=gnu  --max-code-length=80  OPSC_nssolver.cpp
pdflatex alg.tex
# pdflatex equations.tex
rm *.orig
rm *.pyc
rm *.aux
rm *.log