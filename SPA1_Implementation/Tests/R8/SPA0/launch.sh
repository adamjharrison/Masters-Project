#!/bin/bash

export FMS_exe=/Users/adamjharrison/Desktop/Masters/openfms/bin/openfms.zero

for i in {1..200}; do
	cd ic${i}/
	pwd
	$FMS_exe > fms.out
	cd ../
done
