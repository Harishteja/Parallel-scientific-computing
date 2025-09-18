# 1D Convective Travelling Wave — Upwind and QUICK Schemes

## Overview

This project implements a solver for the one-dimensional travelling wave equation:

$$
\frac{\partial u}{\partial t} + c \frac{\partial u}{\partial x} = 0,\qquad c=1.0,
$$

on the domain \(0 \le x \le L\) with \(L=2.0\).  
The initial condition is:

$$
u(x,0) =
\begin{cases}
\sin(4\pi x), & 0 \le x \le 0.5, \\
0, & 0.5 < x \le 2.0,
\end{cases}
$$

with homogeneous Dirichlet boundary conditions:

$$
u(0,t)=0,\qquad u(L,t)=0.
$$

The exact solution is the shifted profile:

$$
u(x,t) = u_0(x-ct).
$$

---

## Numerical Methods

- **Time integration:** Explicit Euler
  
$$
\frac{\partial u}{\partial t} \approx \frac{u_i^{n+1}-u_i^n}{\Delta t}
$$

- **Spatial discretizations:**
  - **First-order Upwind:**
    
$$
\frac{\partial u}{\partial x}\Big|_i \approx \frac{u_i^n - u_{i-1}^n}{\Delta x}.
$$
  - **Third-order QUICK:**
    
$$
  \frac{\partial u}{\partial x}\Big|_i \approx \frac{1}{\Delta x}\left(\tfrac{3}{8}u_i^n - \tfrac{7}{8}u_{i-1}^n + \tfrac{1}{8}u_{i-2}^n + \tfrac{3}{8}u_{i+1}^n\right).
$$

- **Discretization parameters:**
  
  - Grid spacing: \($$\Delta x = 0.002\$$)
  - Time step: \($$\Delta t = 0.0001\$$)
  - CFL number: \(0.05\)

- Near-boundary points where QUICK cannot be applied use the upwind approximation.

---

## Implementations

1. **Serial Codes**  
   - Implemented for both upwind and QUICK schemes.
   - Analytical and numerical solutions compared at \(t = 0,\;0.5,\;1.0\).

2. **MPI Parallel Codes**  
   - Domain decomposition in \(x\) with halo/ghost cells.  
   - Runs tested with \(p=2\) and \(p=4\) processes.  
   - Proper communication of ghost cells ensures correct stencil support for QUICK.

3. **Results & Analysis**  
   - Upwind scheme shows **numerical diffusion** (smeared wave profile).  
   - QUICK scheme provides better accuracy but introduces **oscillations** near discontinuities due to dispersion.  
   - Performance comparisons of serial vs. MPI are included in the report.  
   - All plots show good agreement of numerical solutions with the analytical solution.

# Poisson Solver — Jacobi & Gauss–Seidel (Red–Black) — Serial & MPI

## Overview

This repository contains serial and MPI implementations to solve the Poisson equation

$$
\nabla^2 \phi = -q,\qquad q(x,y) = x^2 + y^2
$$

on a rectangular domain with the boundary conditions shown below. The solvers implemented are:

- Jacobi iterative method (serial & MPI)
- Gauss–Seidel with red–black coloring (serial & MPI)
  
   <img width="600" height="412" alt="Screenshot from 2025-09-18 10-10-58" src="https://github.com/user-attachments/assets/4e05a9d5-e99b-464a-89b2-2a8a4b2f9776" />

Double-precision arithmetic is used throughout. A detailed report with results, convergence analysis, and plots is provided in `report.pdf`.

---

## Governing equation & discretization

The Poisson equation is discretized on a uniform Cartesian mesh with $\Delta x=\Delta y=h$. Using standard 5-point finite differences and Jacobi/Gauss–Seidel iterations, the update formula can be written as

$$
\phi_{i,j}^{(k+1)} = \frac{1}{4}\Big(\phi_{i+1,j}^{(k)} + \phi_{i-1,j}^{(k+1/k)} + \phi_{i,j+1}^{(k)} + \phi_{i,j-1}^{(k+1/k)}\Big) + \frac{h^2}{4}\,q_{i,j},
$$

where the notation $\phi^{(k+1/k)}$ indicates that Jacobi uses the previous iterate for all neighbors, while Gauss–Seidel uses already-updated neighbors when available.

On the boundary at \(x=1\) (right boundary in the domain image) the first derivative is approximated with the second-order one-sided formula:

$$
\phi_i = \frac{4\phi_{i-1} - \phi_{i-2}}{3}.
$$

An initial guess of \(\phi(x,y)=0\) is used for all runs.

---

## Results (summary)

- For \(h=0.1\) (serial Jacobi), the solver converged to the \(10^{-4}\) tolerance in **N\_iter\_serial** iterations (see `report.pdf` for numeric values and plots).
- For \(h=0.01\) (MPI Jacobi), convergence histories and iteration counts are reported for \(p=2,4,8\); parallel solutions match serial solutions to within numerical precision (plots in `report.pdf`).
- Gauss–Seidel (red–black) consistently required **fewer iterations** than Jacobi for the same tolerance; comparative plots and tables are provided in the report.
- Strong-scaling tests for \(h=0.005\) were run with \(p=2,4,8,16\). Speed-up curves and discussion are included in `report.pdf`. Observations:
  - GS red–black generally shows faster convergence (lower iterations) but more synchronization overhead in parallel.
  - Jacobi is easier to parallelize (less synchronization) but needs more iterations.

Please refer to `report.pdf` for full plots, numerical values, and analysis.



