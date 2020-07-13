# FIM_Eikonal
Extension of SCI-Solver_Eikonal for simulation of anisotropic wavefront propagation in the heart.

### Syntax:
```
tetFIM -i inFile.vtu -p seedsFile.txt -o outFile.txt {additional options}
```

### Available options:
```
-i  or -inFile              VTK unstructured grid dataset (.vtu or .vtk)
-p  or -seedsFile           .txt
-o  or -outFile             .bin or .txt
-a  or -anisotropy          default: 1.0
-s  or -speed               default: 1.0
-m  or -maxIterations       default: 500
-b  or -maxBlocks           default: 16
-vb or -maxVertsPerBlock    default: 24
-v  or -verbose
-h  or -help
```

### Hints:
* The seeds file is used to provide IDs of starting points (value zero in the solution).
  It may consist of multiple lines defining different 'seed configurations', each of which is solved individually.
  Within each line, point IDs have to be separated by a space.

* If the output file has the extension '.bin', the output will be binary, otherwise ascii.
  For binary, the first value will be the number of points as int32, followed by the values for each seed configuration as float.
  For ascii, there will be one column for each point in the mesh and one row for each seed configuration.

* Cell data arrays 'Phi' and 'Theta' containing angles in radian are required in the VTK dataset for anisotropy != 1.

* A cell data array 'Speed' may be provided in the VTK dataset to define different speeds for each tetrahedron.

* Try to reduce maxVertsPerBlock if you get 'cudaCheckError() : invalid configuration argument'.
  This parameter determines maxNumTotalTets, which must be smaller than the max number of threads per block for your CUDA device (usually 512 or 1024).
