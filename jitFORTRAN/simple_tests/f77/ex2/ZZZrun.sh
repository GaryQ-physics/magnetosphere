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
