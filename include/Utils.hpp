#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include "Matrix.hpp"

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>

void writeVtkOutput(Matrix<Vec> U, unsigned int Nx, unsigned int Ny, unsigned int timestep){

   // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    double dx = 1 / (Nx + 2);
    double dy = 1 / (Ny + 2);

    double x = 0;
    double y = 0;

    for (int col = 0; col < Ny + 2; col++) {
        for (int row = 0; row < Nx + 2; row++) {
            points->InsertNextPoint(x, y, 0);
            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(Nx + 2, Ny + 2, 1);
    structuredGrid->SetPoints(points);

    // Velocity Array
    vtkDoubleArray *Velocity = vtkDoubleArray::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Print Velocity from bottom to top
    for (int j = 0; j < Ny + 2; j++) {
        for (int i = 0; i < Nx + 2; i++) {
            vel[0] = U(i,j).comp[0];
            vel[1] = U(i,j).comp[1];
            Velocity->InsertNextTuple(vel);
        }
    }

    // Add Velocity to Structured Grid
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname = "output" + "." + std::to_string(timestep) + ".vtk";

    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}

#endif //__UTILS_HPP__
