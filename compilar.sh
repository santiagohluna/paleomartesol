#!/usr/bin/bash

echo "Compilando."

if test ! -d out; then 
    mkdir out
fi

if test ! -d bin; then 
    mkdir bin
fi

cd bin

gfortran -g -fbacktrace -fcheck="all" -Wall -o paleomartesol.x ../src/paleo_martesol.f90 ../src/paleo_martesol_subs.f90

echo "¡Listo!"