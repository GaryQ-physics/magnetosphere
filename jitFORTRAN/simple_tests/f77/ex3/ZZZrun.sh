#!/bin/sh
## run in terminal to execute example
## to view the source files in example, use $ cat *

gfortran -c sub.f
gfortran -c outer_sub.f
gfortran -c main.f
gfortran main.o outer_sub.o sub.o
./a.out
rm a.out
rm *.o

gfortran -c -fPIC sub.f
f2py -c outer_sub.f -I sub.o -m pymod
python -c "import pymod; print(pymod.subv); print(pymod.sub)"
rm sub.o
rm pymod.so
