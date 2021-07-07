#include <cmath>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "cuda.h"

#include <fstream>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>

/******************
 *   6   2   5
 *     \ | /
 *   3 - 0 - 1
 *     / | \
 *   7   4   8
******************/
// D2Q9 LATTICE CONSTANTS
static const int N_DIM = 2;
static const int N_DIRECTIONS = 9;
// Index of lattice velocity in each opposite direction
__constant__ int OPP[N_DIRECTIONS] = {0,   3,   4,   1,   2,   7,   8,   5,   6};
__constant__ int LATTICE_VELOCITIES[N_DIRECTIONS][N_DIM] =  {{ 0, 0},
                                                             { 1, 0},
                                                             { 0, 1},
                                                             {-1, 0},
                                                             { 0,-1},
                                                             { 1, 1},
                                                             {-1, 1},
                                                             {-1,-1},
                                                             { 1,-1}};
__constant__ double LATTICE_WEIGHTS[N_DIRECTIONS] = {4.0/9.0,
                                                     1.0/9.0,
                                                     1.0/9.0,
                                                     1.0/9.0,
                                                     1.0/9.0,
                                                     1.0/36.0,
                                                     1.0/36.0,
                                                     1.0/36.0,
                                                     1.0/36.0};
// TODO
// __constant__ double CS = 1.0/sqrt(3.0);
__constant__ double CS = 0.57735;

typedef struct Node {
    double dir[N_DIRECTIONS];
} Node;

typedef struct Vec {
    double comp[N_DIM];
} Vec;

typedef struct Inputs {
    size_t x_dim;
    size_t y_dim;
    size_t timesteps;
    size_t plot_interval;

    double reynolds;
    double wall_vel;

    std::string case_name;
    std::string dict_name;
} Inputs;



/**********************************************************************************************/






/**********************************************************************************************/



#define gpuErrChk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess)
  {
    fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}


/**********************************************************************************************/






/**********************************************************************************************/

// namespace filesystem = std::filesystem;
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
    dict_name.append("_CUDA_Output");

    int status = mkdir(dict_name.c_str(), 0777);
    if ((status < 0) && (errno != EEXIST)) {exit(EXIT_FAILURE);}
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


void writeVtkOutput(Vec * U, size_t Nx, size_t Ny, size_t timestep, size_t my_rank, std::string case_name, std::string dict_name){

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

    for (size_t col = 1; col < Ny + 1; col++) {
        x = 0;
        { x += dx; }
        for (size_t row = 1; row < Nx + 1; row++) {
            points->InsertNextPoint(x, y, 0);
            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(Nx, Ny, 1);
    structuredGrid->SetPoints(points);

    // Velocity Array
    vtkDoubleArray *Velocity = vtkDoubleArray::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Print Velocity from bottom to top
    for (size_t j = Ny; j > 0; j--) {
        for (size_t i = 1; i < Nx + 1; i++) {
            vel[0] = U[i + j * Nx].comp[0];
            vel[1] = U[i + j * Nx].comp[1];
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


/**********************************************************************************************/






/**********************************************************************************************/


namespace Initialize {

__host__ void initMovingwall(Vec * U, double uMax, size_t Nx){
    /// Initialize Moving Wall
    for(size_t k = 0; k < Nx + 2; ++k){
         U[k + 0 * (Nx + 2)].comp[0] = uMax;
    }
}

// MICROSCOPIC INITIAL CONDITION
// TODO: use CS
__global__ void initDistfunc(Node * fIn, Vec * U, double * RHO, size_t Nx, size_t Ny){

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    double c_u;
    for (size_t q = 0; q < N_DIRECTIONS; ++q) {
        c_u = (U[i + j * (Nx + 2)].comp[0] * LATTICE_VELOCITIES[q][0]
            +  U[i + j * (Nx + 2)].comp[1] * LATTICE_VELOCITIES[q][1]);

        fIn[i + j * (Nx + 2)].dir[q] = RHO[i + j * (Nx + 2)] * LATTICE_WEIGHTS[q]
                              * ( 1 + c_u / (CS * CS) + (c_u * c_u) / (2 * CS * CS * CS * CS)
                              - ( U[i + j * (Nx + 2)].comp[0] * U[i + j * (Nx + 2)].comp[0]
                              +   U[i + j * (Nx + 2)].comp[1] * U[i + j * (Nx + 2)].comp[1] ) / (2 * CS * CS));
    }
}

} // namespace Initialize



/**********************************************************************************************/






/**********************************************************************************************/


namespace Boundary{
/// SETTING BOUNDARIES
// Moving wall
__global__ void setMovingwall(Node * fIn, Vec * U, double * RHO, double uMax, size_t Nx){

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if(j == 0) {
        // Macroscopic Dirichlet boundary conditions
        U[i].comp[0] = uMax;
        U[i].comp[1] = 0.0;
        RHO[i] = (fIn[i].dir[0] + fIn[i].dir[1] + fIn[i].dir[3]
        + 2*(fIn[i].dir[2] + fIn[i].dir[5] + fIn[i].dir[6])) / (1 - U[i].comp[1]);

        // Microscopic  Zou/He boundary conditions
        fIn[i].dir[4] = fIn[i].dir[2] - 2 / 3 * RHO[i] * U[i].comp[1];
        fIn[i].dir[7] = fIn[i].dir[5] + 0.5 * (fIn[i].dir[1] - fIn[i].dir[3])
                    - 0.5 * (RHO[i] * U[i].comp[0]) - 1/6 * (RHO[i] * U[i].comp[1]);
        fIn[i].dir[8] = fIn[i].dir[6] - 0.5 * (fIn[i].dir[1] - fIn[i].dir[3])
                    + 0.5 * (RHO[i] * U[i].comp[0]) - 1/6 * (RHO[i] * U[i].comp[1]);
    }
}


// Set bounce back cells
__global__ void setBounceback(Node * fOut, Node * fIn, size_t Nx, size_t Ny){

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if (j == Ny + 1) {
        for(size_t q = 0; q < N_DIRECTIONS; ++q){
            fOut[i + (Ny + 1) * (Nx + 2)].dir[q] = fIn[i + (Ny + 1) * (Nx + 2)].dir[OPP[q]];
        }
    }

    if (i == 0) {
        for(size_t q = 0; q < N_DIRECTIONS; ++q){
            fOut[0 + j * (Nx + 2)].dir[q] = fIn[0 + j * (Nx + 2)].dir[OPP[q]];
        }
    }

    if (i == Ny + 1) {
        for(size_t q = 0; q < N_DIRECTIONS; ++q){
            fOut[(Nx + 1) + j * (Nx + 2)].dir[q] = fIn[(Nx + 1) + j * (Nx + 2)].dir[OPP[q]];
        }
    }
}

} // namespace Boundary



/**********************************************************************************************/






/**********************************************************************************************/



namespace Processes {

__device__ void computeDensity(const Node * currentNode, double * density){
    *density = 0;
    for (int q = 0; q < N_DIRECTIONS; ++q) {
        *density += currentNode->dir[q];
    }
}

__device__ void computeVelocity(const Node * currentNode, Vec * velocity, double density) {

    for (int d = 0; d < N_DIM; ++d) {
        // Initialize each velocity to zero
        velocity->comp[d] = 0.0;
        for (int q = 0; q < N_DIRECTIONS; ++q) {
           velocity->comp[d] += currentNode->dir[q] * LATTICE_VELOCITIES[q][d];
        }

        velocity->comp[d] = (density == 0) ? 0: velocity->comp[d] / density;
    }
}

__device__ void computefEq(Node * eqField, const Vec * velocity, double density){

	for (int q = 0; q < N_DIRECTIONS; ++q) {

        double c_dot_u = 0.0;
        double u_dot_u = 0.0;
		for (int d = 0; d < N_DIM; ++d) {
            c_dot_u += LATTICE_VELOCITIES[q][d] * velocity->comp[d];
            u_dot_u += velocity->comp[d] * velocity->comp[d];
		}

		eqField->dir[q] =  LATTICE_WEIGHTS[q] * density * (1 + c_dot_u / (CS * CS)
                    + (c_dot_u * c_dot_u)/(2 * CS * CS * CS * CS)
                    - u_dot_u / (2 * CS * CS));
	}
}

__device__ void computePostCollisionDistributions(Node * outputNode, Node * currentNode, const Node * eqField, double omega){

    for (int q = 0; q < N_DIRECTIONS; ++q) {
        outputNode->dir[q] = currentNode->dir[q] - omega *( currentNode->dir[q] - eqField->dir[q] );
    }
}

__global__ void calcQuantities(Node * fIn, Vec * U, double * RHO, size_t Nx, size_t Ny){

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    size_t j = blockIdx.y * blockDim.y + threadIdx.y;
    /// CALCULATE MACROSCOPIC QUANTITIES FOR EACH CELL
    computeDensity( &fIn[i + j * (Nx + 2)], &RHO[i + j * (Nx + 2)]);
    computeVelocity( &fIn[i + j * (Nx + 2)], &U[i + j * (Nx + 2)], RHO[i + j * (Nx + 2)]);
}

__global__ void doCollision(Node * fOut, Node * fIn, Node * fEq, Vec * U, double * RHO,
                                                            double omega, size_t Nx, size_t Ny) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    size_t j = blockIdx.y * blockDim.y + threadIdx.y;
    /// COLLISION STEP
    computefEq( &fEq[i + j * (Nx + 2)], &U[i + j * (Nx + 2)], RHO[i + j * (Nx + 2)] );
    computePostCollisionDistributions( &fOut[i + j * (Nx + 2)], &fIn[i + j * (Nx + 2)], &fEq[i + j * (Nx + 2)], omega );
}

__global__ void doStreaming(Node * fOut, Node * fIn, size_t Nx, size_t Ny){

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    int Cx, Cy;
    for (size_t q = 0; q < N_DIRECTIONS; ++q) {
        Cx = LATTICE_VELOCITIES[q][0];
        Cy = LATTICE_VELOCITIES[q][1];
        fOut[i + j * (Nx + 2)].dir[q]
        = fIn[ (i - Cx + Nx + 2) % (Nx + 2) + ((j - Cy + Ny + 2) % (Ny + 2)) * (Nx + 2) ].dir[q];
    }
}

} // namespace Processes



/**********************************************************************************************/






/**********************************************************************************************/



int main(int argn, char **args){

    if (argn <= 1) {
        std::cout << "ERROR: No input file is provided. Exiting..." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string file_name{args[1]};
    Inputs input;
    Utils::read_inputs(file_name, input);

    /// DIMENSIONS
    size_t Nx      = input.x_dim;        // number of cells in x-direction
    size_t Ny      = input.y_dim;        // number of cells in y-direction
    size_t Nt      = input.timesteps;    // number of time steps
    size_t plot_interval = input.plot_interval;   // How many timesteps for the next plot update

    /// FLOW PARAMETERS
    double Re      = input.reynolds;    // Reynolds number (Main flow parameter)
    double uMax    = input.wall_vel;    // Velocity of lid
    double nu      = uMax * Nx / Re;    // kinematic viscosity (is chosen such that the given Raynolds number is achieved)
    double omega   = 2 / (6 * nu + 1);  // relaxation parameter

    /// Input distribution:
    Node * fIn;
    gpuErrChk(cudaMalloc(&fIn, (Nx + 2) * (Ny + 2) * sizeof(Node)));
    /// Output distributions (post collision):
    Node * fOut;
    gpuErrChk(cudaMalloc(&fOut, (Nx + 2) * (Ny + 2) * sizeof(Node)));
    /// Equilibrium distribution:
    Node * fEq;
    gpuErrChk(cudaMalloc(&fEq, (Nx + 2) * (Ny + 2) * sizeof(Node)));

    /// Macroscopic velocity vector
    Vec * U;
    gpuErrChk(cudaMalloc(&U, (Nx + 2) * (Ny + 2) * sizeof(Vec)));
    /// Density
    double * RHO;
    gpuErrChk(cudaMalloc(&RHO, (Nx + 2) * (Ny + 2) * sizeof(double)));
    // Host side memeory for U
    Vec * Uin;
    Uin = (Vec *)malloc((Nx + 2) * (Ny + 2) * sizeof(Vec));

    // TODO: Initialize all cuda memory to zero

    dim3 gridSize(Nx / 32, Ny / 32);
    dim3 blockSize( 32, 32);
    /// SET INITIAL CONDITIONS
    // Moving wall velocity
    Initialize::initMovingwall(Uin, uMax, Nx);
    // Copy input data to device
    gpuErrChk(cudaMemcpy(U, Uin, (Nx + 2) * (Ny + 2) * sizeof(Vec), cudaMemcpyHostToDevice));
    // Initialize Microscopic quantities
    Initialize::initDistfunc<<<gridSize, blockSize>>>(fIn, U, RHO, Nx, Ny);

/*************************************************************************************
*   SIMULATION
*************************************************************************************/

    for(size_t t_step = 0; t_step < Nt; ++t_step){

        /// CALCULATE MACROSCOPIC QUANTITIES FOR EACH CELL
        Processes::calcQuantities<<<gridSize, blockSize>>>(fIn, U, RHO, Nx, Ny);

        /// SETTING BOUNDARIES
        // Moving wall
        Boundary::setMovingwall<<<gridSize, blockSize>>>(fIn, U, RHO, uMax, Nx);

        /// COLLISION STEP
        Processes::doCollision<<<gridSize, blockSize>>>(fOut, fIn, fEq, U, RHO, omega, Nx, Ny);

        /// SETTING BOUNDARIES
        // Set bounce back cells
        Boundary::setBounceback<<<gridSize, blockSize>>>(fOut, fIn, Nx, Ny);

        // STREAMING STEP
        Processes::doStreaming<<<gridSize, blockSize>>>(fIn, fOut, Nx, Ny);

        if(t_step % plot_interval == 0){
            // TODO: Print message
            gpuErrChk(cudaMemcpy(Uin, U, (Nx + 2) * (Ny + 2) * sizeof(Vec), cudaMemcpyDeviceToHost));
            Utils::writeVtkOutput(Uin, Nx, Ny, t_step, 0, input.case_name, input.dict_name);
        }
    }
}
