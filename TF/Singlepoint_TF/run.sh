#!/bin/bash
file=$1

TCPATH=/opt/ispg/terachem/2023.6/arch/intel2022.2.1-cuda11-mpich3.4.3-openmm7.5
source $TCPATH/setenv-terachem.sh
export TCEXE=$TCPATH/bin/terachem

nohup $TCEXE --inputfile=$file.inp > $file.out &
