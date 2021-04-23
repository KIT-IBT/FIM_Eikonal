//
//  main.cu
//  tetFIM
//
//  Created by Steffen Schuler in July 2020.
//  Copyright © 2020 IBT. All rights reserved.
//

#include <Eikonal.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include "Matrix3.h"

Matrix3<float> calculateTensor(float phi, float theta, float anisotropy)
{
    // angles are defined according to the ISO convention for rotational coordinates, see
    // https://commons.wikimedia.org/w/index.php?title=File:Kugelkoord-def.svg&oldid=378536762
    // the positive z-axis (element 8 of result) represents the initial fiber direction before rotation
    
    float cp = cos(phi);
    float sp = sin(phi);
    Matrix3<float> rotZ(cp, -sp, 0,
                        sp,  cp, 0,
                         0,   0, 1);
    
    float ct = cos(theta);
    float st = sin(theta);
    Matrix3<float> rotY( ct, 0, st,
                          0, 1,  0,
                        -st, 0, ct);
    
    Matrix3<float> rot = rotZ * rotY;
    
    Matrix3<float> result;
    result.SetToIdentityMatrix();
    result(8) = 1.0 / (anisotropy * anisotropy);
    result = rot * result * rot.GetTranspose();
    
    return result;
}

int main(int argc, char *argv[])
{
    std::string inFile;
    std::string seedsFile;
    std::string outFile;
    bool verbose = false;
    float anisotropy = 1.0;
    bool anisotropyProvided = false;
    float speed = 1.0;
    bool speedProvided = false;
    int maxBlocks = 16;
    int maxVertsPerBlock = 24;
    int maxIterations = 2000;
    
    for(int i = 0; i < argc; i++)
    {
        if(argc < 2 || strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "-h") == 0)
        {
            std::cout << std::endl;
            std::cout << "Syntax:" << std::endl << std::endl;
            std::cout << "  tetFIM -i inFile.vtu -p seedsFile.txt -o outFile.txt {additional options}" << std::endl;
            std::cout << std::endl;
            std::cout << "Available options:" << std::endl << std::endl;
            std::cout << "  -i  or -inFile              VTK unstructured grid dataset (.vtu or .vtk)" << std::endl;
            std::cout << "  -p  or -seedsFile           .txt" << std::endl;
            std::cout << "  -o  or -outFile             .bin or .txt" << std::endl;
            std::cout << "  -a  or -anisotropy          default: 1.0" << std::endl;
            std::cout << "  -s  or -speed               default: 1.0" << std::endl;
            std::cout << "  -m  or -maxIterations       default: 2000" << std::endl;
            std::cout << "  -b  or -maxBlocks           default: 16" << std::endl;
            std::cout << "  -vb or -maxVertsPerBlock    default: 24" << std::endl;
            std::cout << "  -v  or -verbose" << std::endl;
            std::cout << "  -h  or -help" << std::endl;
            std::cout << std::endl;
            std::cout << "Hints:" << std::endl << std::endl;
            std::cout << "  The seeds file is used to provide IDs of starting points (value zero in the solution)." << std::endl;
            std::cout << "  It may consist of multiple lines defining different 'seed configurations', each of which is solved individually." << std::endl;
            std::cout << "  Within each line, point IDs have to be separated by a space." << std::endl;
            std::cout << std::endl;
            std::cout << "  If the output file has the extension '.bin', the output will be binary, otherwise ascii." << std::endl;
            std::cout << "  For binary, the first value will be the number of points as int32, followed by the values for each seed configuration as float." << std::endl;
            std::cout << "  For ascii, there will be one column for each point in the mesh and one row for each seed configuration." << std::endl;
            std::cout << std::endl;
            std::cout << "  Cell data arrays 'Phi' and 'Theta' containing angles in radian are required in the VTK dataset for anisotropy != 1." << std::endl;
            std::cout << std::endl;
            std::cout << "  A cell data array 'Speed' may be provided in the VTK dataset to define different speeds for each tetrahedron." << std::endl;
            std::cout << std::endl;
            std::cout << "  Try to reduce maxVertsPerBlock if you get 'cudaCheckError() : invalid configuration argument'." << std::endl;
            std::cout << "  This parameter determines maxNumTotalTets, which must be smaller than the max number of threads per block for your CUDA device (usually 512 or 1024)." << std::endl;
            std::cout << std::endl;
            return(EXIT_SUCCESS);
        }
        else if((strcmp(argv[i], "-inFile") == 0 || strcmp(argv[i], "-i") == 0) && argc > i+1)
            inFile = argv[i+1];
        else if((strcmp(argv[i], "-seedsFile") == 0 || strcmp(argv[i], "-p") == 0) && argc > i+1)
            seedsFile = argv[i+1];
        else if((strcmp(argv[i], "-outFile") == 0 || strcmp(argv[i], "-o") == 0) && argc > i+1)
            outFile = argv[i+1];
        else if((strcmp(argv[i], "-anisotropy") == 0 || strcmp(argv[i], "-a") == 0) && argc > i+1)
        {
            anisotropy = atof(argv[i+1]);
            anisotropyProvided = true;
        }
        else if((strcmp(argv[i], "-speed") == 0 || strcmp(argv[i], "-s") == 0) && argc > i+1)
        {
            speed = atof(argv[i+1]);
            speedProvided = true;
        }
        else if((strcmp(argv[i], "-maxIterations") == 0 || strcmp(argv[i], "-m") == 0) && argc > i+1)
            maxIterations = atoi(argv[i+1]);
        else if((strcmp(argv[i], "-maxBlocks") == 0 || strcmp(argv[i], "-b") == 0) && argc > i+1)
            maxBlocks = atoi(argv[i+1]);
        else if((strcmp(argv[i], "-maxVertsPerBlock") == 0 || strcmp(argv[i], "-vb") == 0) && argc > i+1)
            maxVertsPerBlock = atoi(argv[i+1]);
        else if(strcmp(argv[i], "-verbose") == 0 || strcmp(argv[i], "-v") == 0)
            verbose = true;
    }
    
    if(inFile.empty())
    {
        std::cout << "ERROR:  Missing parameter -inFile or -i" << std::endl;
        return(EXIT_FAILURE);
    }
    if(seedsFile.empty())
    {
        std::cout << "ERROR:  Missing parameter -seedsFile or -p" << std::endl;
        return(EXIT_FAILURE);
    }
    if(outFile.empty())
    {
        std::cout << "ERROR:  Missing parameter -outFile or -o" << std::endl;
        return(EXIT_FAILURE);
    }
    
    std::cout << "INFO:   maxIterations is " << maxIterations << std::endl;
    std::cout << "INFO:   maxBlocks is " << maxBlocks << std::endl;
    std::cout << "INFO:   maxVertsPerBlock is " << maxVertsPerBlock << std::endl;
    
    std::ofstream outFileStream(outFile);
    if(outFileStream.is_open())
        outFileStream.close();
    else
    {
        std::cout << "ERROR:  Could not open output file" << std::endl;
        return(EXIT_FAILURE);
    }
    
    int pos = inFile.find_last_of(".");
    std::string extension = inFile.substr(pos+1, inFile.size()-pos);
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
    
    std::cout << "STATUS: Loading input file " << inFile << "..." << std::endl;
    
    vtkSmartPointer<vtkUnstructuredGrid> uGrid;
    if(extension == "vtk")
    {
        vtkSmartPointer<vtkDataSetReader> vtkReader = vtkSmartPointer<vtkDataSetReader>::New();
        vtkReader->SetFileName(inFile.c_str());
        if(!vtkReader->OpenVTKFile())
        {
            std::cout << "ERROR:  VTK file could not be opened - does it exist?" << std::endl;
            return(EXIT_FAILURE);
        }
        vtkReader->Update();
        if(!vtkReader->IsFileUnstructuredGrid())
        {
            std::cout << "ERROR:  VTK file does not contain an unstructured grid" << std::endl;
            return(EXIT_FAILURE);
        }
        uGrid = vtkReader->GetUnstructuredGridOutput();
    }
    else if(extension == "vtu")
    {
        vtkSmartPointer<vtkXMLUnstructuredGridReader> vtuReader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        if(!vtuReader->CanReadFile(inFile.c_str()))
        {
            std::cout << "ERROR:  VTU file could not be opened - does it exist?" << std::endl;
            return(EXIT_FAILURE);
        }
        vtuReader->SetFileName(inFile.c_str());
        vtuReader->Update();
        uGrid = vtuReader->GetOutput();
    }
    
    // create cell list
    int numCells = (int)uGrid->GetNumberOfCells();
    int* cells = new int[4*numCells];
    for(int i = 0; i < numCells; i++)
    {
        vtkCell* cell = uGrid->GetCell(i);
        if(cell->GetCellType() != VTK_TETRA)
        {
            std::cout << "ERROR:  Only linear tetrahedrons are supported" << std::endl;
            return(EXIT_FAILURE);
        }
        for(int j = 0; j < 4; j++)
            cells[4*i+j] = (int)cell->GetPointId(j);
    }
    
    // create point list
    int numPoints = uGrid->GetNumberOfPoints();
    float* points = new float[3*numPoints];
    for(int i = 0; i < numPoints; i++)
    {
        for(int j = 0; j < 3; j++)
            points[3*i+j] = (float)uGrid->GetPoint(i)[j];
    }
    
    std::cout << "STATUS: Loading seeds file " << seedsFile << "..." << std::endl;
    
    std::vector< std::vector<int> > seeds;
    std::ifstream seedsFileStream(seedsFile);
    if(!seedsFileStream.is_open())
    {
        std::cout << std::endl << "ERROR:  Could not open seeds file" << std::endl;
        return(EXIT_FAILURE);
    }
    std::string line;
    while(std::getline(seedsFileStream, line))
    {
        std::istringstream iss(line);
        std::vector<int> s;
        int i;
        while(iss >> i)
            s.push_back(i);
        if(s.size() > 0)
            seeds.push_back(s);
    }
    seedsFileStream.close();
    
    if(seeds.size() == 0)
    {
        std::cout << "ERROR:  No seed configuration found" << std::endl;
        return(EXIT_FAILURE);
    }
    std::cout << "INFO:   " << seeds.size() << " seed configuration(s) found" << std::endl;
    
    // create slowness values
    std::vector<float> slownessMtx;
    bool isotropicCase = true;
    vtkCellData* cellData = uGrid->GetCellData();
    
    std::vector<float> slownessVals(numCells, 1.0/(speed*speed));
    if(speedProvided)
        std::cout << "INFO:   Using global speed of " << speed << std::endl;
    else if(cellData->HasArray("Speed"))
    {
        vtkDataArray* speedArray = cellData->GetArray("Speed");
        for(int i = 0; i < numCells; i++)
        {
            float s = speedArray->GetComponent(i, 0);
            slownessVals[i] = 1.0/(s*s);
        }
        std::cout << "INFO:   Using local speeds from cell data array 'Speed'" << std::endl;
    }
    else
        std::cout << "INFO:   Using default global speed of " << speed << std::endl;
    
    std::vector<float> anisotropyVals;
    if(anisotropyProvided && anisotropy != 1.0)
    {
        anisotropyVals.resize(numCells);
        std::fill(anisotropyVals.begin(), anisotropyVals.end(), anisotropy);
        isotropicCase = false;
        std::cout << "INFO:   Anisotropic case" << std::endl;
        std::cout << "INFO:   Using global anisotropy ratio of " << anisotropy << std::endl;
    }
    else if(cellData->HasArray("Anisotropy"))
    {
        anisotropyVals.resize(numCells);
        vtkDataArray* anisotropyArray = cellData->GetArray("Anisotropy");
        for(int i = 0; i < numCells; i++)
            anisotropyVals[i] = anisotropyArray->GetComponent(i, 0);
        isotropicCase = false;
        std::cout << "INFO:   Anisotropic case" << std::endl;
        std::cout << "INFO:   Using local anisotropy ratios from cell data array 'Anisotropy'" << std::endl;
    }
    
    if(isotropicCase)
    {
        // isotropic case: only 1 value for each cell in slownessMtx
        slownessMtx = slownessVals;
        std::cout << "INFO:   Isotropic case" << std::endl;
    }
    else
    {
        // anisotropic case: 6 values for each cell in slownessMtx,
        // representing one triangular part of the symmetric matrix
        if(cellData->HasArray("Phi") && cellData->HasArray("Theta"))
        {
            slownessMtx.resize(6*numCells);
            vtkDataArray* phiArray = cellData->GetArray("Phi");
            vtkDataArray* thetaArray = cellData->GetArray("Theta");
            for(int i = 0; i < numCells; i++)
            {
                float phi = phiArray->GetComponent(i, 0);
                float theta = thetaArray->GetComponent(i, 0);
                Matrix3<float> M = calculateTensor(phi, theta, anisotropyVals[i]);
                slownessMtx[6*i+0] = slownessVals[i] * M(0);
                slownessMtx[6*i+1] = slownessVals[i] * M(1);
                slownessMtx[6*i+2] = slownessVals[i] * M(2);
                slownessMtx[6*i+3] = slownessVals[i] * M(4);
                slownessMtx[6*i+4] = slownessVals[i] * M(5);
                slownessMtx[6*i+5] = slownessVals[i] * M(8);
            }
            std::cout << "INFO:   Using cell data arrays 'Phi' and 'Theta' to define orientation of anisotropy" << std::endl;
        }
        else
        {
            std::cout << "ERROR:  Cell data arrays 'Phi' and 'Theta' are required for anisotropic case" << std::endl;
            return(EXIT_FAILURE);
        }
    }
    
    std::cout << "STATUS: Converting mesh..." << std::endl;
    
    TetMesh* tetMesh = new TetMesh();
    tetMesh->init(points, numPoints,
                  NULL, 0, // trilist, numtri
                  cells, numCells,
                  NULL, // attrlist
                  slownessMtx,
                  verbose);
    tetMesh->need_neighbors(verbose);
    tetMesh->need_adjacenttets(verbose);
    tetMesh->need_tet_virtual_tets(verbose);
    
    meshFIM3dEikonal* fim3d = new meshFIM3dEikonal;
    fim3d->SetMesh(tetMesh);
    
    std::cout << "STATUS: Partitioning..." << std::endl;
    
    fim3d->GraphPartition_METIS2(maxBlocks, maxVertsPerBlock, verbose);
    fim3d->m_numBlock = maxBlocks;
    fim3d->PartitionTets(maxBlocks, verbose);
    
    std::cout << "INFO:   maxNumTotalTets is " << fim3d->m_maxNumTotalTets << std::endl;
    
    // create a copy of fim3d, so that it can be restored for multiple runs
    meshFIM3dEikonal fim3d_copy = *fim3d;
    
    pos = outFile.find_last_of(".");
    extension = outFile.substr(pos+1, outFile.size()-pos);
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
    bool binaryOutput = (extension == "bin") ? true : false;
    
    if(binaryOutput)
    {
        std::cout << "INFO:   Using binary output format" << std::endl;
        outFileStream.open(outFile, std::ofstream::out | std::ofstream::app | std::ofstream::binary);
    }
    else
    {
        std::cout << "INFO:   Using ascii output format" << std::endl;
        outFileStream.open(outFile, std::ofstream::out | std::ofstream::app);
    }
    if(!outFileStream.is_open())
    {
        std::cout << "ERROR:  Could not open output file" << std::endl;
        return(EXIT_FAILURE);
    }
    if(binaryOutput)
        outFileStream.write(reinterpret_cast<const char*>(&numPoints), sizeof(int));
    else
        outFileStream.precision(7);
    
    for(int i = 0; i < seeds.size(); i++)
    {
        std::cout << "STATUS: Solving - seed configuration " << i+1 << " of " << seeds.size() << "..." << std::endl;
        
        // restore fim3d
        *fim3d = fim3d_copy;
        
        fim3d->SetSeedPoint(seeds[i]);
        std::vector< std::vector<float> > iterVals;
        iterVals = fim3d->GenerateData(maxIterations, verbose);
        std::vector<float> vals = iterVals.back();
        
        if(binaryOutput)
        {
            for(int i = 0; i < vals.size(); i++)
                outFileStream.write(reinterpret_cast<const char*>(&vals[i]), sizeof(float));
        }
        else
        {
            for(int i = 0; i < vals.size(); i++)
                outFileStream << vals[i] << " ";
            outFileStream << std::endl;
        }
        
        //std::vector< std::vector<float> > exportVals;
        //exportVals.push_back(iterVals.back());
        //fim3d->writeVTK(exportVals);
    }
    
    outFileStream.close();
    
    return(EXIT_SUCCESS);
}
