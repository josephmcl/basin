#!/bin/bash
module load gcc/13.1.0
module load intel-oneapi-mkl/2023.1.0
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKLROOT/lib/intel64
export LD_LIBRARY_PATH
