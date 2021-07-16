# LBSIM
LBSIM is a CUDA high-performance Computational Fluid Dynamics simulation code based on the Lattice Boltzmann Method. It's a code for simulating low-speed incompressible flows. It is currently in its first development steps, capable of running the two-dimensional Lid-Driven Cavity benchmark with different grid sizes and flow properties. The code can both run sequentially on the CPU and in parallel on the GPU.

In the LBM part of the code, we have calculated the equilibrium in each cell using the macroscopic quantities and then performed the collision according to the BGK-operator. We directly used the Zou-He scheme for the boundary condition in moving no-slip walls and pure bounce back for fixed walls. 

We are creating two programs, and both are D2Q9 Lattice Boltzmann codes. The first one runs in series and is written in C++. The second one is parallelized with a GPU. Functions in the GPU are the initialization of distribution function, setting the boundary conditions, collision, and streaming. The code will be more noticeably faster when it is run for a high-resolution problem.

## Initialization:
You can change different parameters such as Reynolds number, number of the cells, time steps, and wall velocities to get different results using the cases/LidDrivenCavity/LidDrivenCavity.dat file.

## Requirments:
Ensure you have all requirements before you start building the code since it is written in C++ and CUDA. This code is designed to run on **Linux**. We strongly recommend using Linux for compilation, computation, and post-processing. Therefore, The following prerequisites are demanded:

1. A supported version of Linux with a GCC compiler and tool-chain
2. C++ and its editor, preferably [Visual Studio Code](https://code.visualstudio.com/docs/setup/linux)
3. C++ build system generator, [Cmake](https://cmake.org/install/)
4. You need a CUDA-capable GPU, which can be checked here ([Link](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#verify-you-have-cuda-enabled-system))
5. To use your GPU, you’ll need [NVIDIA CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit)
6. The Visualization Toolkit, VTK 7 or higher


## How to use it:

### Install
```
git clone git@github.com:bamaratunga/lbsim.git
cd lbsim
```

### Compile
```
mkdir build && cd build
cmake ..
make install
```

### Run
To run the code sequentially on the CPU, you can use the following command: 

`lbsim_cpp ../cases/LidDrivenCavity/LidDrivenCavity.dat`

To run the code in parallel on GPU, you can use the following command: 

`lbsim_cuda ../cases/LidDrivenCavity/LidDrivenCavity.dat`

## Contributors
Binu Amaratunga, Erfan Mashayekh, Hessel Juliust




