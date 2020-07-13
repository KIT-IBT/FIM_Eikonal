cmake_minimum_required(VERSION 2.8)
project(tetFIM)

############### VTK ###############

# Select only VTK components (libraries to be linked to) that are needed.
# (Removing the COMPONENTS argument from find_package links all of them.)
# See folder structure in https://vtk.org/doc/nightly/html/files.html
# to find out to which component a VTK class belongs.
set(VTK_COMPONENTS
    vtkCommonCore
    vtkCommonDataModel
    vtkCommonExecutionModel
    vtkCommonMath
    vtkCommonMisc
    vtkCommonSystem
    vtkFiltersCore
    vtkFiltersExtraction
    vtkFiltersFlowPaths
    vtkFiltersGeneral
    vtkFiltersGeometry
    vtkFiltersModeling
    vtkFiltersParallel
    vtkIOCore
    vtkIOGeometry
    vtkIOLegacy
    vtkIOPLY
    vtkIOXML
)
find_package(VTK REQUIRED COMPONENTS ${VTK_COMPONENTS} NO_MODULE)
if(NOT VTK_FOUND)
    message(FATAL_ERROR "VTK could not be found. Try setting VTK_DIR.")
else()
    message(STATUS "Found VTK: " ${VTK_DIR})
endif()
include(${VTK_USE_FILE})

############### CUDA ###############

find_package(CUDA REQUIRED)

####################################

set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../SCI-Solver_Eikonal/src")
set(BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../SCI-Solver_Eikonal/build")

include_directories(
    ${SRC_DIR}
    ${SRC_DIR}/core3d/include
    ${SRC_DIR}/core2d/include
    ${BUILD_DIR}/externals/cusp/src
)

CUDA_ADD_EXECUTABLE(tetFIM src/main.cu)
TARGET_LINK_LIBRARIES(tetFIM ${VTK_LIBRARIES} ${BUILD_DIR}/libEIKONAL_CORE.a ${BUILD_DIR}/core3d/libEikonal-3D_CORE.a ${BUILD_DIR}/externals/metis/lib/libmetis.a)

