#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include "Matrix.hpp"
#include "Definitions.hpp"

#include <fstream>
#include <filesystem>
#include <algorithm>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>

namespace Utils{


void set_file_names(std::string file_name, std::string& case_name, std::string& dict_name) {
    std::string temp_dir;
    bool case_name_flag = true;

    for (int i = file_name.size() - 1; i > -1; --i) {
        if (file_name[i] == '/') {
            case_name_flag = false;
        }
        if (case_name_flag) {
            case_name.push_back(file_name[i]);
        }
    }

    for (int i = file_name.size() - case_name.size() - 1; i > -1; --i) {
        temp_dir.push_back(file_name[i]);
    }

    std::reverse(case_name.begin(), case_name.end());
    std::reverse(temp_dir.begin(), temp_dir.end());

    case_name.erase(case_name.size() - 4);
    dict_name = temp_dir;
    dict_name.append(case_name);
    dict_name.append("_Output");

    // Create output directory
    std::filesystem::path folder(dict_name);
    try {
        std::filesystem::create_directory(folder);
    } catch (const std::exception &e) {
        std::cerr << "Output directory could not be created." << std::endl;
        std::cerr << "Make sure that you have write permissions to the "
                     "corresponding location"
                  << std::endl;
    }
}


void read_inputs(std::string file_name, Inputs& input) {
    // Read input parameters
    const int MAX_LINE_LENGTH = 1024;
    std::ifstream file(file_name);

    if (file.is_open()) {
        std::string var;
        while (!file.eof() && file.good()) {
            file >> var;
            if (var[0] == '#') { /* ignore comment lines*/
                file.ignore(MAX_LINE_LENGTH, '\n');
            } else {
                if (var == "imax") file >> input.x_dim;
                if (var == "jmax") file >> input.y_dim;
                if (var == "timesteps") file >> input.timesteps;
                if (var == "plot_int") file >> input.plot_interval;
                if (var == "Re") file >> input.reynolds;
                if (var == "wall_vel") file >> input.wall_vel;
            }
        }
    }

    file.close();

    set_file_names(file_name, input.case_name, input.dict_name);
}


void writeVtkOutput(Matrix<Vec>& U, size_t Nx, size_t Ny, size_t timestep, size_t my_rank, std::string case_name, std::string dict_name){

   // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    double dx = 1.0 / Nx;
    double dy = 1.0 / Ny;

    double x = 0;
    double y = 0;

    { y += dy; }
    { x += dx; }

    for (size_t col = 0; col < Ny + 1; col++) {
        x = 0;
        { x += dx; }
        for (size_t row = 0; row < Nx + 1; row++) {
            points->InsertNextPoint(x, y, 0);
            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(Nx + 1, Ny + 1, 1);
    structuredGrid->SetPoints(points);

    // Velocity Array
    vtkDoubleArray *Velocity = vtkDoubleArray::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Print Velocity from bottom to top
    for (size_t j = 1; j < Ny + 1; j++) {
        for (size_t i = 1; i < Nx + 1; i++) {
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
    std::string outputname =
    dict_name + '/' + case_name + "_" + std::to_string(my_rank) + "." + std::to_string(timestep) + ".vtk";
    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}


} // namespace Utils

#endif //__UTILS_HPP__
