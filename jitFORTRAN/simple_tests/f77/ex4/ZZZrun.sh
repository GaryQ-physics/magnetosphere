#!/bin/sh
## run in terminal to execute example
## to view the source files in example, use $ cat *

gfortran -c subroutine.f
gfortran -c main.f
gfortran main.o subroutine.o
./a.out
rm a.out
rm *.o

gfortran -c -fPIC subroutine.f
f2py -c main.f -I subroutine.o -m pymod
python -c "import pymod; print(pymod.subv); print(pymod.sub)"
rm subroutine.o
rm pymod.so



