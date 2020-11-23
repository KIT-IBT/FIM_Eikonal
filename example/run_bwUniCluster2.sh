#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=32000
#SBATCH --gres=gpu:1

module load compiler/clang
module load devel/cuda

../tetFIM/build/tetFIM -i box.vtu -p seeds.txt -o result.bin -a 3
