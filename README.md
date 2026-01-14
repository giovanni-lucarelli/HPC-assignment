# HPC Exam Project: Hybrid MPI/OpenMP 5-Point Stencil

This project presents the design, implementation, and performance evaluation
of a hybrid MPI/OpenMP parallel solver for the two-dimensional heat diffusion
equation based on the 5-point stencil method.

The main objectives of the project are:

1. To implement a hybrid parallel solution using MPI for distributed-memory
   domain decomposition and OpenMP for shared-memory parallelism, starting
   from the provided code templates.

2. To instrument the code in order to separately measure computation and
   communication times.

3. To perform a comprehensive scalability study, including:
   - OpenMP thread scaling
   - Strong scaling
   - Weak scaling

The program simulates the time evolution of the two-dimensional heat diffusion
equation on a rectangular grid using an explicit finite-difference scheme:

$$
\partial_t u = \alpha \left( \partial_x^2 u + \partial_y^2 u \right)
$$

which leads to the following discrete update rule:

$$
u^{(t+1)}_{i,j} = (1 - 4\alpha)\,u^{(t)}_{i,j}
+ \alpha \sum_{\langle i,j \rangle} u^{(t)}_{i,j}
$$

where $ i \in \{1, \dots, N_x-1\} $ and $ j \in \{1, \dots, N_y-1\} $
denote the internal grid points.

The implementation supports both periodic and non-periodic boundary
conditions, allows the injection of an arbitrary number of heat sources
configured at runtime, and provides optional output of the grid state
and energy statistics during the simulation.

> **Note:** This README is a summary of the project structure, usage instructions,
> and implementation details. For more details, analysis, and results,
> please refer to the [slides](slide/main.pdf) of the presentation and 
> the [Jupyter notebook](plot/results.ipynb).


## Repository Structure

```
HPC-assignment/
├── code-leonardo/                      # HPC cluster code and scripts
│   ├── Makefile                        # Build configuration
│   ├── data/                           # Raw experimental data (CSV files)
│   ├── include/
│   │   └── stencil_template_parallel.h # Parallel implementation header
│   ├── script/
│   │   ├── compile.sbatch                  # compilation job script
│   │   ├── run_scaling.sbatch              # Strong/weak scaling job script
│   │   ├── run_threads.sbatch              # Thread scaling job script
│   │   ├── run_vanilla.sbatch              # Baseline execution script
│   │   ├── submit_strong.sh                # Strong scaling submission
│   │   ├── submit_threads.sh               # Thread scaling submission
│   │   └── submit_weak.sh                  # Weak scaling submission
│   └── src/
│       └── stencil_template_parallel.c # Hybrid MPI/OpenMP implementation
├── code-local/                         # Local development code
├── data/                               # Experimental data
│   ├── strong_*.csv                    # Strong scaling results
│   ├── threads_*.csv                   # Thread scaling results
│   └── weak.csv                        # Weak scaling results
├── doc/                                # Documentation
│   ├── guide_serial.md                 # Serial implementation guide
│   └── node_info.md                    # HPC node specifications
├── plot/                               # Jupyter notebooks for visualization
│   ├── evolution-parallel.ipynb        # Parallel performance analysis
│   ├── evolution-serial.ipynb          # Serial performance analysis
│   └── results.ipynb                   # Main results visualization
└── slide/                              # Presentation

```

The first version of the code has been developed locally, and is the one saved in `code-local/`. The version runned on the cluster, specifically on Leonardo CINECA supercomputer, has been developed in SSH protocol and is the one in `code-leonardo/`. In the folder are also stored the job files and the submission bash scripts used for the scaling study.

In order to reproduce the analysis on the cluster

1. compile the program using the job file `compile.sbatch` which internally loads the dependencies and calls the makefile

2. submit the scaling experiments (threads, strong, weak) using the bash scripts `submit_*.sh`

The program will save the results in `data/results.csv`.

   > **Note:** During the scalability studies, grid printing and output generation are disabled by default via a preprocessor directive to ensure accurate and uncontaminated timing measurements. Output can be enabled at compile time by setting `ENABLE_OUTPUT=1` in the `compile.sbatch` file.

## Command Line Interface (CLI)

The program accepts command-line options to configure the simulation.
All parameters are optional and have default values.

### Usage

```bash
mpirun -np <Ntasks> ./stencil [options]
```
### Options

* `-x <int>` : grid size in the x-direction (default: 10000)
* `-y <int>` : grid size in the y-direction (default: 10000)
* `-e <int>` : total number of heat sources (default: 4)
* `-E <double>` : energy per heat source (default: 1.0)
* `-n <int>` : number of iterations (default: 1000)
* `-p <0|1>` : enable periodic boundary conditions
* `-o <0|1>` : enable output of energy statistics
* `-v <int>` : verbosity level
* `-h` : print help message and exit

## Code Implementation and Details



This section briefly describes the role of the main functions in the program.

* `initialize(...)` handles the global initialization of the program:

    * parsing of command-line arguments
    * domain decomposition among MPI processes
    * identification of neighboring processes (halo exchange)
    * allocation of data structures and communication buffers
    * initialization of heat sources

* `update_internal(...)` and `update_boundary(...)` execute the computational kernel of the **5-point stencil**, respectively for the internal grid points and the boundary grid points, i.e., those that depend on halo region data from neighboring MPI processes. 

* `output_energy_stat(...)` computes the total energy of the system or the avg energy per grid point, performing global reduction in MPI.

* The functions `dump(...)`, `gather_global_plane(...)`, and `dump_global(...)` are used to produce visual output from the simulation. They support both local and global data dumping by collecting the distributed subdomains managed by MPI processes and reconstructing the global solution when needed. These routines are mainly intended for debugging, and visualization.

## Parallelization Strategy

The program employs a hybrid parallelization strategy combining MPI for distributed memory parallelism and OpenMP for shared memory parallelism within each MPI process.

### MPI Parallelization

After decomposing the global domain into subdomains, each MPI process is responsible for updating its local subdomain. Specifically, each process exchanges halo regions with its neighboring processes to update boundary values. This is done using non-blocking MPI communication to overlap computation and communication. The general flow is as follows:

1. Each MPI process initializes its local subdomain and heat sources.
2. For each iteration:
   - Exchange halo regions with neighboring processes using non-blocking MPI calls.
   - Update internal grid points using the 5-point stencil method.
   - Wait for halo exchanges and update boundary grid points.

### OpenMP parallelization

OpenMP is used inside each MPI process to parallelize the stencil update
over the local subdomain, and in general to parallelize loops that operate on the local data structures. The computational loop over grid points is shared among threads, improving cache locality and fully exploiting multi-core architectures.