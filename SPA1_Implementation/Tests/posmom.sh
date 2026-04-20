#!/bin/bash
# Adds all positions and momenta of the initial conditions to a csv file

echo "position,momentum" > posmom.csv

for i in {1..200}; do
    awk 'NR==3 {x=$2} NR==5 {print x "," $1}' R8/SPA0/ic${i}/Geometry.dat
done >> posmom.csv