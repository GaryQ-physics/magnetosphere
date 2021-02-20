#!/bin/sh
## run in terminal to execute example
## to view the source files in example, use $ cat *

gfortran -c subroutine.f
gfortran -c main.f
gfortran main.o subroutine.o
./a.out
rm *.o
rm a.out

f2py -c subroutine.f -m pymod
python main.py
rm pymod.so

python jit.py

