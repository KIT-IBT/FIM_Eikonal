#!/bin/bash

var=$(hostname)
if [ ${var:0:4} == 'uc2n' ]; then
    # we are on bwUniCluster2
    module load compiler/clang
    module load devel/cuda/11.0
    module load devel/cmake
fi

mkdir SCI-Solver_Eikonal/build
pushd SCI-Solver_Eikonal/build
cmake ../src
make
popd

mkdir tetFIM/build
pushd tetFIM/build
cmake .. # you might want to specify the path to VTK here, e.g.: cmake .. -D VTK_DIR=../../vtk-8.2.0/build
make
popd
