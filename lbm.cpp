#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>

//***********************
//   6   2   5
//     \ | /
//   3 - 0 - 1
//     / | \
//   7   4   8
//**********************
// D2Q9 LATTICE CONSTANTS
static const int N_DIM = 2;
static const int N_DIRECTIONS = 9;
// Index of lattice velocity in each opposite direction
static const int OPP[N_DIRECTIONS] = {0,   3,   4,   1,   2,   7,   8,   5,   6};
static const int LATTICE_VELOCITIES[N_DIRECTIONS][N_DIM] =  {{ 0, 0},
                                                             { 1, 0},
                                                             { 0, 1},
                                                             {-1, 0},
                                                             { 0,-1},
                                                             { 1, 1},
                                                             {-1, 1},
                                                             {-1,-1},
                                                             { 1,-1}};
static const double LATTICE_WEIGHTS[N_DIRECTIONS] = {4.0/9.0,
                                                     1.0/9.0,
                                                     1.0/9.0,
                                                     1.0/9.0,
                                                     1.0/9.0,
                                                     1.0/36.0,
                                                     1.0/36.0,
                                                     1.0/36.0,
                                                     1.0/36.0};
static const double CS = 1.0/sqrt(3.0);

typedef struct Node {
    double dir[N_DIRECTIONS];
} Node;

typedef struct Vec {
    double comp[N_DIM];
} Vec;


/*************************************************************************************/



/*************************************************************************************/


template <typename T> class Matrix {

  public:
    Matrix<T>() = default;

    /// Constructor without an initial value.
    Matrix<T>(size_t i_max, size_t j_max) : _imax(i_max), _jmax(j_max) {
        _container = (T *)malloc(_imax * _jmax * sizeof(T));
    }

    /// Element access and modify using index
    T &operator()(size_t i, size_t j) {
        return _container[j * _imax + i];
    }

    /// Access of the size of the structure
    size_t size() const { return _imax * _jmax; }

    /// get the number of elements in x direction
    size_t imax() const { return _imax; }

    /// get the number of elements in y direction
    size_t jmax() const { return _jmax; }

  private:
    /// Number of elements in x direction
    size_t _imax;
    /// Number of elements in y direction
    size_t _jmax;

    /// Data container
    T * _container;
};

void writeVtkOutput(Matrix<Vec> U, unsigned int Nx, unsigned int Ny, unsigned int timestep);

/*************************************************************************************/



/*************************************************************************************/


void computeDensity(const Node * currentNode, double * density){
    *density = 0;
    for (int q = 0; q < N_DIRECTIONS; ++q) {
        *density += currentNode->dir[q];
    }
}

void computeVelocity(const Node * currentNode, Vec * velocity, double density) {

    for (int d = 0; d < N_DIM; ++d) {
        // Initialize each velocity to zero
        velocity->comp[d] = 0.0;
        for (int q = 0; q < N_DIRECTIONS; ++q) {
           velocity->comp[d] += currentNode->dir[q] * LATTICE_VELOCITIES[q][d];
        }

        velocity->comp[d] = (density == 0) ? 0: velocity->comp[d] / density;
    }
}

void computefEq(Node * eqField, const Vec * velocity, double density){

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


void computePostCollisionDistributions(Node * currentNode, const Node * eqField, double omega){

    for (int q = 0; q < N_DIRECTIONS; ++q) {
        currentNode->dir[q] = currentNode->dir[q] - omega *( currentNode->dir[q] - eqField->dir[q] );
    }
}


/*************************************************************************************/



/*************************************************************************************/


int main(){

    /// DIMENSIONS
    int Nx      = 4;         // number of cells in x-direction
    int Ny      = 4;         // number of cells in y-direction
    int Nt      = 1;        // number of time steps
    int plot_interval = 100;   // How many timesteps for the next plot update

    /// FLOW PARAMETERS
    double Re      = 300;               // Reynolds number (Main flow parameter)0000
    double uMax    = 0.08;              // Velocity of lid
    double nu      = uMax * Nx / Re;    // kinematic viscosity (is chosen such that the given Raynolds number is achieved)
    double omega   = 2 / (6 * nu + 1);  // relaxation parameter

    /// Input distribution:
    Matrix<Node> fIn  = Matrix<Node>(Nx + 2, Ny + 2);
    /// Output distributions (post collision):
    Matrix<Node> fOut = Matrix<Node>(Nx + 2, Ny + 2);
    /// Equilibrium distribution:
    Matrix<Node> fEq  = Matrix<Node>(Nx + 2, Ny + 2);

    /// Macroscopic velocity vector
    Matrix<Vec> U = Matrix<Vec>(Nx + 2, Ny + 2);
    /// Density
    Matrix<double> RHO = Matrix<double>(Nx + 2, Ny + 2);

    /// Zero initialize all
    for (size_t i = 0; i < Nx + 2; ++i) {
        for (size_t j = 0; j < Ny + 2; ++j) {
            for(size_t q = 0; q < N_DIRECTIONS; ++q){
                fIn(i, j).dir[q] = 0.0;
                fOut(i, j).dir[q] = 0.0;
                fEq(i, j).dir[q] = 0.0;
            }
            for(size_t d = 0; d < N_DIM; ++d){
                U(i, j).comp[d] = 0.0;
            }
            RHO(i, j) = 1.0;
        }
    }
    /// Initialize Moving Wall
    for(size_t i = 0; i < Nx + 2; ++i){
         U(i, 0).comp[0] = uMax;
    }

    // MICROSCOPIC INITIAL CONDITION
    // TODO: Separate function, use C_S
    double c_u;
    for (size_t i = 0; i < Nx + 2; ++i) {
        for (size_t j = 0; j < Ny + 2; ++j) {
            for (size_t q = 0; q < N_DIRECTIONS; ++q) {
                c_u = 3 * (U(i,j).comp[0] * LATTICE_VELOCITIES[q][0]
                        +  U(i,j).comp[1] * LATTICE_VELOCITIES[q][1]);

                fIn(i,j).dir[q] = RHO(i, j) * LATTICE_WEIGHTS[q]
                                   * ( 1 + c_u + 0.5 * (c_u * c_u)
                                      - 1.5 * ( U(i,j).comp[0] * U(i,j).comp[0] + U(i,j).comp[1] * U(i,j).comp[1] ));
            }
        }
    }


/*************************************************************************************/





/*************************************************************************************/


    for(size_t t_step = 0; t_step < Nt; ++t_step){


        /// CALCULATE MACROSCOPIC QUANTITIES FOR EACH CELL
        for (size_t i = 0; i < Nx + 2; ++i) {
            for (size_t j = 0; j < Ny + 2; ++j) {
                computeDensity( &fIn(i, j), &RHO(i, j));
                computeVelocity( &fIn(i, j), &U(i, j), RHO(i, j));
            }
        }


        /// SETTING BOUNDARIES
        // Moving wall
        for (size_t i = 0; i < Nx + 2; ++i) {
            // Macroscopic Dirichlet boundary conditions
            U(i, 0).comp[0] = uMax;
            U(i, 0).comp[1] = 0.0;
            RHO(i, 0) =   (fIn(i, 0).dir[0] + fIn(i, 0).dir[1] + fIn(i, 0).dir[3]
                      + 2*(fIn(i, 0).dir[2] + fIn(i, 0).dir[5] + fIn(i, 0).dir[6])) / (1 - U(i, 0).comp[1]);

            // Microscopic  Zou/He boundary conditions
            fIn(i, 0).dir[4] = fIn(i, 0).dir[2] - 2 / 3 * RHO(i, 0) * U(i, 0).comp[1];
            fIn(i, 0).dir[7] = fIn(i, 0).dir[5] + 0.5 * (fIn(i, 0).dir[1] - fIn(i, 0).dir[3])
                        - 0.5 * (RHO(i, 0) * U(i, 0).comp[0]) - 1/6 * (RHO(i, 0) * U(i, 0).comp[1]);
            fIn(i, 0).dir[8] = fIn(i, 0).dir[6] + 0.5 * (fIn(i, 0).dir[1] - fIn(i, 0).dir[3])
                        + 0.5 * (RHO(i, 0) * U(i, 0).comp[0]) - 1/6 * (RHO(i, 0) * U(i, 0).comp[1]);
        }


        /// COLLISION STEP
        for (size_t i = 0; i < Nx + 2; ++i) {
            for (size_t j = 0; j < Ny + 2; ++j) {
                computefEq( &fEq(i, j), &U(i, j), RHO(i, j));
                computePostCollisionDistributions( &fOut(i, j), &fEq(i,j), omega);
            }
        }


        /// SETTING BOUNDARIES
        // Set bounce back cells
        for(size_t i = 0; i < Nx + 2; ++i){
            for(size_t q = 0; q < N_DIRECTIONS; ++q){
                fOut(i, Ny + 1).dir[q] = fIn(i, Ny + 1).dir[OPP[q]];
            }
        }
        for(size_t j = 1; j < Ny + 2; ++j){
            for(size_t q = 0; q < N_DIRECTIONS; ++q){
                fOut(0, j).dir[q] = fIn(0, j).dir[OPP[q]];
                fOut(Nx + 1, j).dir[q] = fIn(Nx + 1, j).dir[OPP[q]];
            }
        }


        // STREAMING STEP
        int Cx, Cy;
        for (size_t j = 0; j < Ny + 2; ++j) {
            for (size_t i = 0; i < Nx + 2; ++i) {
                for (size_t q = 0; q < N_DIRECTIONS; ++q) {
                    Cx = LATTICE_VELOCITIES[q][0];
                    Cy = LATTICE_VELOCITIES[q][1];
                    fIn((i + Cx + Nx) % Nx, (j + Cy + Ny) % Ny).dir[q]
                                                = fOut(i, j).dir[q];
                }
            }
        }

        writeVtkOutput(U, Nx, Ny, t_step);
    }
}

/*************************************************************************************/





/*************************************************************************************/


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
