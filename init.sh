#!/bin/bash
module load gcc/12.2.0
MKLROOT=/packages/intel/19/linux/mklls
export MKLROOT
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/jmclaug2/petsc/arch-linux-c-opt/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel64
export LD_LIBRARY_PATH

#gcc a.c -I$MKLROOT/include -L$MKLROOT/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm