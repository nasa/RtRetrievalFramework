rm -f *.exe *.mod
#gfortran -c -Wall Chapman_BOA.f90
gfortran -o CF.exe -Wall Chapman_BOA.f90 test_trans.f90

