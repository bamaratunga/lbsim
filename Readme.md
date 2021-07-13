# LBSIM
LBSIM is a CUDA high-performance Computational Fluid Dynamics simulation code based on the Lattice Boltzmann Method. LBSIM is a code for simulating low-speed incompressible flows. It is currently in its first development steps, capable of running the two-dimensional Lid-Driven Cavity benchmark with different grid sizes and flow properties. The code can both run sequentially on the CPU and in parallel on the GPU, and it can be switched using different commands. 


In the LBM part of the code, we have calculated the equilibrium in each cell using the macroscopic quantities and then performed the collision according to the BGK-operator. We directly used the Zou-He scheme for the boundary condition in moving no-slip walls and pure bounce back for fixed walls. 


We are creating two programs, and both are D2Q9 Lattice Boltzmann codes. The first one runs in serial and is written in C++. We implemented some namespaces that are not suggested in the LBM worksheet, which suppose to be in C. The second one is parallelized with a GPU. Here, CUDA API creates a kernel, which means a function will be executed parallelly in the GPU threads.  Functions in the GPU are the initialization of distribution function, setting the boundary conditions, collision, and streaming. The code will be more noticeably faster when it is run for a high-resolution problem.

## Improvements
The purpose of developing LBSIM was to answer whether we can create a fluid dynamics solver that is faster than the codes that solve Navier Stokes equations in their computations. To solve Navier Stokes equations, we often use implicit methods and store huge matrices that may lead to more computational time and energy consumption. To overcome such a problem, we decided to test the Lattice Boltzmann equation as the solver to reduce such costs since this method is explicit and does not store all variables for its calculations. Besides, as GPU programming is helpful for a large amount of data to process, we tried to implement it with CUDA and check whether we can get a boost in the computations. The performance results of the code are compared with the CFD Lab worksheet code for the Lid-Driven Cavity with size 50x50. The time it takes to reach the steady-state condition is computed for the methods in the following table.

| Code             | # of Threads | Time (s)  | Speedup |
| ---------------- | ------------ | --------- | ------- |
| CFD Lab          | 1            | 23.98     | -       |
| CFD Lab + MPI    | 2            | 14.1      | 1.7     |
| Sequential LBSIM | 1            | 8.59      | 2.79    |
| LBSIM on GPU *   | Na           | 3.87      | 6.19    |

* NVIDIA GTX 1060: The number of CUDA cores is 1152

We also calculated the speedup for different sizes of our GPU implementations:

| Size        |	Sequential only time (s) |	GPU time (s) | Speedup          |
| ----------- | ------------------------- | ------------------ | ---------------- |
| 50x50	      | 5.03	                 | 2.39	             | 2.10460251046025 |
| 100x100     | 23.8808	                 | 9.49319            | 2.51557168875794 |
| 300x300     | 291.828	                 | 104.871	     | 2.78273307205996 |
| 500x500     | 838.177	                 | 249.946	     | 3.35343234138574 |
| 800x800     | 2442.92	                 | 711.123	     | 3.43529881609792 |
| 1000x1000   | Na	                 | 1236.38	     | Na               |

In bigger problems, the parts of the code that we parallelize would increase, resulting in greater speed-ups at the end.

## Initialization:
There is a file in the repository called LidDrivenCavity.dat in which you can change different parameters such as Reynolds number, number of the cells, time steps, and wall velocities to get different results.

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
git clone https://gitlab.lrz.de/ge75rud/lbsim.git
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






