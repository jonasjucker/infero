#!/bin/bash -l

. ../dev/daint_env.sh

INFERO_ROOT='../../infero_build/builds/infero'
LIBS='-linferoapi -linferof -lnetcdff -lnetcdf -lmpichf90'
FCFLAGS="-cpp -Wall -pedantic -fbacktrace -O2 -g -ffree-line-length-256 -I${NETCDF_DIR}/include -I${INFERO_ROOT}/module -Wl,-rpath,${INFERO_ROOT}/lib -L${NETCDF_DIR}/lib -I${MPICH_DIR}/include"
LDFLAGS="-L${INFERO_ROOT}/lib -L${NETCDF_DIR}/lib -L${MPICH_DIR}/lib"

ftn -o multiple_fields_model.exe $FCFLAGS $LDFLAGS $LIBS ecrad_ml.f90 
ftn -o performance.exe $FCFLAGS $LDFLAGS $LIBS performance.f90 
