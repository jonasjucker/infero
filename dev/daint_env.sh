#!/bin/bash -l

module load PrgEnv-gnu
module load daint-gpu
module switch gcc gcc/9.3.0
module load cray-netcdf/4.7.4.4
module load CMake
module load cray-mpich
