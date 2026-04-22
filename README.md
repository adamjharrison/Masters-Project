# Masters Project
Code used for Masters Project: "Method Development of Generalised Ab Initio Multiple Spawning"

Supervisor: Professor Basile Curchod

## SPA1_Implementation
Code used to test the saddle point approximation on the SOC matrix elements between TBFs in Generalised Ab Initio Multiple Spawning.

The `SPA1_SOC_model()` function within the `SPA1.f90` file is purely included to show the main code produced for this project, however this will not run by itself as it is implemented into the larger OpenFMS nonadiabatic dynamics codebase. More details on how this function is interfaced within OpenFMS can be seen in the merged commit here: https://github.com/ispg-group/openfms/commit/bd061019bf5d5c313afa324edec250b0ef357f05

### Tests

This folder contains the inputs and `N.dat` outputs produced for testing by comparing the SPA0 and SPA1 implementations of different variants of the 1D model system. It also contains the positions and momentum distribution of the 200 initial conditions in `posmom.csv`.

### Verification

`SPA1_verify.ipynb` is a jupyter notebook that calculates the first order term of the SPA1 using the Wolfram Client library for Python using input parameters of TBFs for the cases discussed in the thesis, used to verify the SPA1 implementation in OpenFMS. Also calculates exact values of the SOC matrix elements between TBFs. Note that this requires the Wolfram Engine to be installed to run as the Wolfram Client interfaces with the Engine to do the calculations.

## 1D_model
`Matrixeval.f90` contains some initial code used to learn Fortran, used to evaluate the matrix elements of the 1D model system shown in Equation 24 of G. Granucci, M. Persico and G. Spighi, *J. Chem. Phys.*, 2012, **137**, 22A501. This code was not used for any testing purposes as the system was already implemented in OpenFMS but was used for the model plots and learning Fortran.

## TF

`Rotation.py` is a simple script for rotating the momentum and the position of a TF geometry by 90 degrees with respect to the x axis.

### CSThresh

Contains inputs and fms.out file, as well as relevant outputs for plots, for GAIMS simulations where `CSThresh` was varied from 0.0 to 0.1.

### OMin_Parent

Contains inputs and `fms.out` file, as well as relevant outputs for plots, for GAIMS simulations where `OMin_Parent` was varied from 0.2 to 0.8. Also contained within this folder are the simulations on the transformed geometry with `OMin_Parent` set to 0.4.


### Singlepoint_TF

Contains a single point energy calcualtion of the thioformaldehyde structure using SA-CASSCF.

### TF_nomom

Contains inputs and `fms.out` file, as well as relevant outputs for plots, for initial GAIMS simulations of thioformaldehyde structure with zero momentum.