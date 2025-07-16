#!/bin/bash
module load gcc/13.1.0
module load intel-oneapi-mkl/2023.1.0
module load intel-oneapi-compilers/2023.1.0
module load intel-oneapi-mpi/2021.9.0

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel64
export LD_LIBRARY_PATH
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/packages/spack/spack-rhel8/opt/spack/linux-rhel8-broadwell/gcc-13.1.0/intel-oneapi-dal-2023.1.0-rwo3dn4gikgsiubrqa4gxaxlqsvn66xx/compiler/2023.1.0/linux/compiler/lib/intel64_lin
export LD_LIBRARY_PATH

#/home/jmclaug2/intel/oneapi/advisor/2025.1/lib64/libstdc++.so.6
export LD_LIBRARY_PATH=/home/jmclaug2/intel/oneapi/advisor/2025.1/lib64:$LD_LIBRARY_PATH