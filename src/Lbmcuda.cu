#include <cmath>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <chrono>
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
__constant__ double LATTICE_WEIGHTS[N_DIRECTIONS][2] = {{4.0, 9.0},
                                                        {1.0, 9.0},
                                                        {1.0, 9.0},
                                                        {1.0, 9.0},
                                                        {1.0, 9.0},
                                                        {1.0, 36.0},
                                                        {1.0, 36.0},
                                                        {1.0, 36.0},
                                                        {1.0, 36.0}};

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
    bool output_results;

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
                if (var == "output_results"){
                    std::string state;
                    file >> state;
                    if (state == "on") input.output_results = true;
                    else input.output_results = false;
                }
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

    for (size_t col = 0; col < Ny; col++) {
        x = 0;
        for (size_t row = 0; row < Nx; row++) {
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
    for (size_t j = 0; j < Ny; j++) {
        for (size_t i = 0; i < Nx; i++) {
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

__global__ void initMovingwall(Node * fIn, Node * fOut, Node * fEq, Vec * U, double * RHO, double uMax, size_t Nx, size_t Ny){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if((i < Nx + 2) && (j < Ny + 2)) {
        /// Zero initialize all data
        for(size_t d = 0; d < N_DIM; d++){
            U[i + j * (Nx + 2)].comp[d] = 0.0;
        }
        for(size_t q = 0; q < N_DIRECTIONS; q++){
            fIn[i + j * (Nx + 2)].dir[q] = 0.0;
            fOut[i + j * (Nx + 2)].dir[q] = 0.0;
            fEq[i + j * (Nx + 2)].dir[q] = 0.0;
        }
        /// Initialize density to 1
        RHO[i + j * (Nx + 2)] = 1.0;

        /// Initialize Moving Wall
        if(j == Ny + 1){
             U[i + (Ny + 1) * (Nx + 2)].comp[0] = uMax;
        }
    }
}

// MICROSCOPIC INITIAL CONDITION
__global__ void initDistfunc(Node * fIn, Vec * U, double * RHO, size_t Nx, size_t Ny){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if((i < Nx + 2) && (j < Ny + 2)) {
        double c_u;
        double u_u;
        for (size_t q = 0; q < N_DIRECTIONS; ++q) {
            c_u = (U[i + j * (Nx + 2)].comp[0] * LATTICE_VELOCITIES[q][0]
                +  U[i + j * (Nx + 2)].comp[1] * LATTICE_VELOCITIES[q][1]);
            u_u = (U[i + j * (Nx + 2)].comp[0] * U[i + j * (Nx + 2)].comp[0]
                 + U[i + j * (Nx + 2)].comp[1] * U[i + j * (Nx + 2)].comp[1]);

            double weight = __ddiv_rn(LATTICE_WEIGHTS[q][0], LATTICE_WEIGHTS[q][1]);
            fIn[i + j * (Nx + 2)].dir[q] = RHO[i + j * (Nx + 2)] * weight
                                  * ( 1 + 3 * c_u + 4.5 * c_u * c_u - 1.5 * u_u );
        }
    }
}

} // namespace Initialize



/**********************************************************************************************/






/**********************************************************************************************/


namespace Boundary{
/// SETTING BOUNDARIES
// Moving wall
__global__ void setMovingwall(Node * fIn, Vec * U, double * RHO, double uMax, size_t Nx, size_t Ny){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if((i < Nx + 2) && (j < Ny + 2)) {
        if(j == Ny + 1) {
            // Macroscopic Dirichlet boundary conditions
            U[i + (Ny + 1) * (Nx + 2)].comp[0] = uMax;
            U[i + (Ny + 1) * (Nx + 2)].comp[1] = 0.0;
            RHO[i + (Ny + 1) * (Nx + 2)] = (fIn[i + (Ny + 1) * (Nx + 2)].dir[0] + fIn[i + (Ny + 1) * (Nx + 2)].dir[1] + fIn[i + (Ny + 1) * (Nx + 2)].dir[3]
            + 2*(fIn[i + (Ny + 1) * (Nx + 2)].dir[2] + fIn[i + (Ny + 1) * (Nx + 2)].dir[5] + fIn[i + (Ny + 1) * (Nx + 2)].dir[6])) / (1 - U[i + (Ny + 1) * (Nx + 2)].comp[1]);

            // Microscopic  Zou/He boundary conditions
            fIn[i + (Ny + 1) * (Nx + 2)].dir[4] = fIn[i + (Ny + 1) * (Nx + 2)].dir[2] - 2 / 3 * RHO[i + (Ny + 1) * (Nx + 2)] * U[i + (Ny + 1) * (Nx + 2)].comp[1];
            fIn[i + (Ny + 1) * (Nx + 2)].dir[7] = fIn[i + (Ny + 1) * (Nx + 2)].dir[5] + 0.5 * (fIn[i + (Ny + 1) * (Nx + 2)].dir[1] - fIn[i + (Ny + 1) * (Nx + 2)].dir[3])
                        - 0.5 * (RHO[i + (Ny + 1) * (Nx + 2)] * U[i + (Ny + 1) * (Nx + 2)].comp[0]) - 1/6 * (RHO[i + (Ny + 1) * (Nx + 2)] * U[i + (Ny + 1) * (Nx + 2)].comp[1]);
            fIn[i + (Ny + 1) * (Nx + 2)].dir[8] = fIn[i + (Ny + 1) * (Nx + 2)].dir[6] - 0.5 * (fIn[i + (Ny + 1) * (Nx + 2)].dir[1] - fIn[i + (Ny + 1) * (Nx + 2)].dir[3])
                        + 0.5 * (RHO[i + (Ny + 1) * (Nx + 2)] * U[i + (Ny + 1) * (Nx + 2)].comp[0]) - 1/6 * (RHO[i + (Ny + 1) * (Nx + 2)] * U[i + (Ny + 1) * (Nx + 2)].comp[1]);
        }
    }
}


// Set bounce back cells
__global__ void setBounceback(Node * fOut, Node * fIn, size_t Nx, size_t Ny){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if((i < Nx + 2) && (j < Ny + 2)) {
        if (j == 0) {
            for(size_t q = 0; q < N_DIRECTIONS; ++q){
                fOut[i + 0 * (Nx + 2)].dir[q] = fIn[i + 0 * (Nx + 2)].dir[OPP[q]];
            }
        }

        if ((i == 0) && (j < Ny + 1)) {
            for(size_t q = 0; q < N_DIRECTIONS; ++q){
                fOut[0 + j * (Nx + 2)].dir[q] = fIn[0 + j * (Nx + 2)].dir[OPP[q]];
            }
        }

        if ((i == Nx + 1) && (j < Ny + 1)) {
            for(size_t q = 0; q < N_DIRECTIONS; ++q){
                fOut[(Nx + 1) + j * (Nx + 2)].dir[q] = fIn[(Nx + 1) + j * (Nx + 2)].dir[OPP[q]];
            }
        }
    }
}

} // namespace Boundary



/**********************************************************************************************/






/**********************************************************************************************/



namespace Processes {

__device__ void computeDensity(Node * currentNode, double * density){
    *density = 0;
    for (int q = 0; q < N_DIRECTIONS; ++q) {
        *density = __dadd_rn(*density, currentNode->dir[q]);
    }
}

__device__ void computeVelocity(Node * currentNode, Vec * velocity, double density) {

    for (int d = 0; d < N_DIM; ++d) {
        // Initialize each velocity to zero
        velocity->comp[d] = 0.0;
        for (int q = 0; q < N_DIRECTIONS; ++q) {
            velocity->comp[d] = __dadd_rn(velocity->comp[d], currentNode->dir[q] * LATTICE_VELOCITIES[q][d]);
        }
        velocity->comp[d] = (fabs(density) < 0.01) ? 0: __ddiv_rn(velocity->comp[d], density);
    }
}

__device__ void computefEq(Node * eqField, const Vec * velocity, double density){

	for (size_t q = 0; q < N_DIRECTIONS; ++q) {

        double c_dot_u = 0.0;
        double u_dot_u = 0.0;
		for (size_t d = 0; d < N_DIM; ++d) {
            c_dot_u = __dadd_rn(c_dot_u, __dmul_rn(LATTICE_VELOCITIES[q][d], velocity->comp[d]));
            u_dot_u = __dadd_rn(u_dot_u, __dmul_rn(velocity->comp[d], velocity->comp[d]));
		}
        double c_dot_u_2 = __dmul_rn(c_dot_u, c_dot_u);
        double weight = __ddiv_rn(LATTICE_WEIGHTS[q][0], LATTICE_WEIGHTS[q][1]);
        double weight_density = __dmul_rn(weight, density);
		eqField->dir[q] = __dmul_rn(weight_density, (1 + 3.0 * c_dot_u + 4.5 * c_dot_u_2 - 1.5 * u_dot_u));
	}
}

__device__ void computePostCollisionDistributions(Node * outputNode, Node * currentNode, Node * eqField, double omega){

    for (size_t q = 0; q < N_DIRECTIONS; ++q) {
        // Accuracy hotspot
        outputNode->dir[q] = currentNode->dir[q] - omega * (currentNode->dir[q] - eqField->dir[q]);
        // outputNode->dir[q] = (1 - omega) * currentNode->dir[q] + omega * eqField->dir[q];
        // outputNode->dir[q] = __dmul_rn((1 - omega), currentNode->dir[q]) + __dmul_rn(omega, eqField->dir[q]);
    }
}

__global__ void doCollision(Node * fOut, Node * fIn, Node * fEq, Vec * U, double * RHO,
                                                            double omega, size_t Nx, size_t Ny) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if((i < Nx + 2) && (j < Ny + 2)) {
        /// CALCULATE MACROSCOPIC QUANTITIES FOR EACH CELL
        computeDensity( &fIn[i + j * (Nx + 2)], &RHO[i + j * (Nx + 2)]);
        computeVelocity( &fIn[i + j * (Nx + 2)], &U[i + j * (Nx + 2)], RHO[i + j * (Nx + 2)]);
        /// COLLISION STEP
        computefEq( &fEq[i + j * (Nx + 2)], &U[i + j * (Nx + 2)], RHO[i + j * (Nx + 2)] );
        computePostCollisionDistributions( &fOut[i + j * (Nx + 2)], &fIn[i + j * (Nx + 2)], &fEq[i + j * (Nx + 2)], omega );
    }
}

__global__ void doStreaming(Node * fOut, Node * fIn, size_t Nx, size_t Ny){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if((i < Nx + 2) && (j < Ny + 2)) {
        int Cx, Cy;
        for (size_t q = 0; q < N_DIRECTIONS; ++q) {
            Cx = LATTICE_VELOCITIES[q][0];
            Cy = LATTICE_VELOCITIES[q][1];
            fOut[i + j * (Nx + 2)].dir[q]
            = fIn[ (i - Cx + Nx + 2) % (Nx + 2) + ((j - Cy + Ny + 2) % (Ny + 2)) * (Nx + 2) ].dir[q];
        }
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
    assert(Nx == Ny);
    size_t Nt      = input.timesteps;    // number of time steps

    bool output_results = input.output_results;
    size_t plot_interval = input.plot_interval;   // How many timesteps for the next plot update

    /// FLOW PARAMETERS
    double Re      = input.reynolds;    // Reynolds number (Main flow parameter)
    double uMax    = input.wall_vel;    // Velocity of lid
    double nu      = uMax * (Nx + 2) / Re; // kinematic viscosity (is chosen such that the given Raynolds number is achieved)
    double omega   = 2 / (6 * nu + 1);  // relaxation parameter

    auto start = std::chrono::high_resolution_clock::now();

    /// Input distribution:
    Node * fIn = NULL;
    gpuErrChk(cudaMalloc(&fIn, (Nx + 2) * (Ny + 2) * sizeof(Node)));
    /// Output distributions (post collision):
    Node * fOut = NULL;
    gpuErrChk(cudaMalloc(&fOut, (Nx + 2) * (Ny + 2) * sizeof(Node)));
    /// Equilibrium distribution:
    Node * fEq = NULL;
    gpuErrChk(cudaMalloc(&fEq, (Nx + 2) * (Ny + 2) * sizeof(Node)));

    /// Macroscopic velocity vector
    Vec * U = NULL;
    gpuErrChk(cudaMalloc(&U, (Nx + 2) * (Ny + 2) * sizeof(Vec)));
    /// Density
    double * RHO = NULL;
    gpuErrChk(cudaMalloc(&RHO, (Nx + 2) * (Ny + 2) * sizeof(double)));

    Vec * Uin = NULL;
    Uin = (Vec *)malloc((Nx + 2) * (Ny + 2) * sizeof(Vec));

    /// SET-UP GRID
    size_t gSize_x = floor((Nx + 2)/8) + 1;
    size_t gSize_y = floor((Ny + 2)/16) + 1;
    size_t bSize_x = min(int(Nx + 2), 8);
    size_t bSize_y = min(int(Ny + 2), 16);

    dim3 gridSize(gSize_x, gSize_y);
    dim3 blockSize(bSize_x, bSize_y);
    /// SET INITIAL CONDITIONS
    // Moving wall velocity and zero initialize other memory locations
    Initialize::initMovingwall<<<gridSize, blockSize>>>(fIn, fOut, fEq, U, RHO, uMax, Nx, Ny);
    // gpuErrChk(cudaDeviceSynchronize());
    // Initialize Microscopic quantities
    Initialize::initDistfunc<<<gridSize, blockSize>>>(fIn, U, RHO, Nx, Ny);
    // gpuErrChk(cudaDeviceSynchronize());

/*************************************************************************************
*   SIMULATION
*************************************************************************************/
    std::cout << "Running simulation..." << std::endl;
    for(size_t t_step = 0; t_step < Nt; ++t_step){

        /// SETTING BOUNDARIES - Set Moving wall
        Boundary::setMovingwall<<<gridSize, blockSize>>>(fIn, U, RHO, uMax, Nx, Ny);
        // gpuErrChk(cudaDeviceSynchronize());

        /// COLLISION STEP
        Processes::doCollision<<<gridSize, blockSize>>>(fOut, fIn, fEq, U, RHO, omega, Nx, Ny);
        // gpuErrChk(cudaDeviceSynchronize());

        /// SETTING BOUNDARIES - Set bounce back cells
        Boundary::setBounceback<<<gridSize, blockSize>>>(fOut, fIn, Nx, Ny);
        gpuErrChk(cudaDeviceSynchronize());

        /// STREAMING STEP
        Processes::doStreaming<<<gridSize, blockSize>>>(fIn, fOut, Nx, Ny);
        // gpuErrChk(cudaDeviceSynchronize());

        if(output_results && t_step % plot_interval == 0){
            std::cout << "Writing vtk file at t = " << t_step << std::endl;
            gpuErrChk(cudaMemcpy(Uin, U, (Nx + 2) * (Ny + 2) * sizeof(Vec), cudaMemcpyDeviceToHost));
            // gpuErrChk(cudaDeviceSynchronize());
            Utils::writeVtkOutput(Uin, Nx+2, Ny+2, t_step, 0, input.case_name, input.dict_name);
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << std::endl;
    std::cout << "LBSIM CUDA VERSION" << std::endl;
    std::cout << "Simulating " << input.case_name << std::endl;
    std::cout << "Dimensions = " << Nx + 2 << " x " << Ny + 2 << std::endl;
    std::cout << "No. of iterations = " << Nt << std::endl;
    std::cout << "Execution time = " << duration.count() / 1e6 << "s" << std::endl;
}
