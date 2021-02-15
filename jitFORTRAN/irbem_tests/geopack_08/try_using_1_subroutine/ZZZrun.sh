#!/bin/sh
## run in terminal to execute example
## to view the source files in example, use $ cat *

gfortran -c subroutine.f
gfortran -c main.f
gfortran main.o subroutine.o
./a.out
rm *.o

f2py -c subroutine.f -m pymod
python -c "import pymod; print(pymod.meadv)"
rm pymod.so
rm a.out
