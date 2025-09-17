# Parallel Unsteady Heat Conduction Solver using MPI and OpenMP

## Overview

This project aims to solve unsteady heat conduction problems in 2D and 3D domains using parallel computing techniques. The solution approach employs a **hybrid parallel programming paradigm** that combines:

- **Message Passing Interface (MPI)** for inter-node communication within a cluster
- **Open Multi-Processing (OpenMP)** for intra-node parallelism

Additionally, implementations using only MPI and only OpenMP are also provided to compare against the hybrid method.

The objective is to compute the temperature distribution within the considered domains by solving the governing equations using appropriate discretization methods and iterative solvers. Performance metrics such as **efficiency**, **execution time**, and **speed-up** are evaluated and plotted.

---

## Governing Equation

The generalized heat conduction equation in Cartesian coordinates is given by:

$$
\frac{\partial}{\partial x}\left(k \frac{\partial T}{\partial x}\right) + \frac{\partial}{\partial y}\left(k \frac{\partial T}{\partial y}\right) + \frac{\partial}{\partial z}\left(k \frac{\partial T}{\partial z}\right) + \dot{Q_g} = \rho C_p \frac{\partial T}{\partial t}
$$

Where:


- \( T \) is the temperature
- \( k \) is thermal conductivity
- \( \rho \) is density
- \( C_p \) is specific heat capacity
- \( \dot{Q_g} \) is the volumetric heat generation rate

Assuming constant properties, the equation simplifies to:

$$
\alpha \left(\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} + \frac{\partial^2 T}{\partial z^2}\right) = \frac{\partial T}{\partial t}
$$

Where \( \alpha = \frac{k}{\rho C_p} \) is the thermal diffusivity.

---

## Discretization

### Spatial Discretization

Using the **Central Difference Scheme**, the second derivatives in each spatial direction are approximated as follows:

For the \( x \)-direction:

$$
\frac{\partial^2 T}{\partial x^2} \approx \frac{T_{i-1,j,k}^{n} - 2T_{i,j,k}^{n} + T_{i+1,j,k}^{n}}{\Delta x^2}
$$

For the \( y \)-direction:

$$
\frac{\partial^2 T}{\partial y^2} \approx \frac{T_{i,j-1,k}^{n} - 2T_{i,j,k}^{n} + T_{i,j+1,k}^{n}}{\Delta y^2}
$$

For the \( z \)-direction:

$$
\frac{\partial^2 T}{\partial z^2} \approx \frac{T_{i,j,k-1}^{n} - 2T_{i,j,k}^{n} + T_{i,j,k+1}^{n}}{\Delta z^2}
$$

---

### Temporal Discretization

Using **Forward Time Marching**, the time derivative is approximated as:

$$
\frac{\partial T}{\partial t} \approx \frac{T_{i,j,k}^{n+1} - T_{i,j,k}^{n}}{\Delta t}
$$

---

### Final Discretized Equation

Defining the coefficients:

$$
k_1 = \frac{\alpha \Delta t}{\Delta x^2}, \quad k_2 = \frac{\alpha \Delta t}{\Delta y^2}, \quad k_3 = \frac{\alpha \Delta t}{\Delta z^2}
$$

The final update equation for the 3D domain becomes:

$$
T_{i,j,k}^{n+1} = T_{i,j,k}^{n} + k_1 \left(T_{i-1,j,k}^{n} - 2T_{i,j,k}^{n} + T_{i+1,j,k}^{n}\right) + k_2 \left(T_{i,j-1,k}^{n} - 2T_{i,j,k}^{n} + T_{i,j+1,k}^{n}\right) + k_3 \left(T_{i,j,k-1}^{n} - 2T_{i,j,k}^{n} + T_{i,j,k+1}^{n}\right)
$$

For the 2D domain, the update equation simplifies to:

$$
T_{i,j}^{n+1} = T_{i,j}^{n} + k_1 \left(T_{i-1,j}^{n} - 2T_{i,j}^{n} + T_{i+1,j}^{n}\right) + k_2 \left(T_{i,j-1}^{n} - 2T_{i,j}^{n} + T_{i,j+1}^{n}\right)
$$

---

## Implementation Details

### Hybrid Approach
- MPI is used for communication between nodes in the cluster.
- OpenMP is used for parallelizing computations within each node.

### Serial and Parallel Implementations
- A serial solver is implemented as a baseline.
- MPI-only and OpenMP-only solvers are implemented for comparison.
- The hybrid MPI-OpenMP approach is evaluated for scalability and efficiency.

### Performance Metrics
- Efficiency
- Execution time
- Speed-up

Plots and comparative analysis are provided to highlight the advantages and limitations of each approach.

---

