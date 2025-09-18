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


