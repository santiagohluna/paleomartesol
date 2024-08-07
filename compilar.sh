#!/usr/bin/bash

echo "Compilando."

mkdir bin

cd bin

gfortran -g -fbacktrace -fcheck="all" -Wall -o paleomartesol.x ../src/paleo_martesol.f90 ../src/paleo_martesol_subs.f90

echo "¡Listo!"