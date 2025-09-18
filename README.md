# Parallel-scientific-computing

## Hybrid_MPI_OpenMP
Unsteady heat conduction problems in 2D and 3D domains will be solved and parallelized using a
hybrid parallel programming paradigm combining Message Passing Interface (MPI) and Open
Multi-Processing (OpenMP) and then parallel programs for MPI and OpenMP solely will be
coded. This hybrid program approach aims to use MPI for within and across the nodes of a cluster
communication and OpenMP for communication within each node.Different levels of parallelism
can be achieved by combining MPI and OpenMP parallelization.To evaluate the performance of the
hybrid code, parameters like efficiency, execution time, speed up will be evaluated and plotted to
Compare the hybrid MPI-OpenMP approach against solely MPI and OpenMP implementations.
## MPI
This repository contains MPI-based parallel solvers for a variety of scientific computing problems. Implementations include Jacobi and Gauss–Seidel (red–black) methods for Poisson’s equation, as well as first-order upwind and third-order QUICK schemes for convection problems. Both serial and MPI versions are provided to compare accuracy, convergence, and scalability. For detailed source codes and performance results, please see the MPI/ directory, which contains the codes and the full report (report.pdf).
