# Numerical Derivative Calculation Using Padé Scheme

## Overview

This repository contains implementations to compute the derivative of the function:

$$
f(x) = \sin(5x), \quad 0 \leq x \leq 3
$$

using high-order accurate Padé schemes. The derivative is computed using:

- **Fourth-order accurate Padé scheme** for interior points
- **Third-order accurate one-sided Padé scheme** near boundaries

### Governing Equations

For interior points (\(i = 1, 2, ..., N-1\)):

$$
f'_{i+1} + 4f'_i + f'_{i-1} = \frac{3}{h} (f_{i+1} - f_{i-1})
$$

For the left boundary:

$$
f'_0 + 2f'_1 = \frac{1}{h}\left(-\frac{5}{2}f_0 + 2f_1 + \frac{1}{2}f_2\right)
$$

For the right boundary:

$$
f'_N + 2f'_{N-1} = \frac{1}{h}\left(\frac{5}{2}f_N - 2f_{N-1} - \frac{1}{2}f_{N-2}\right)
$$

where \(h\) is the grid spacing and \(N\) is the total number of grid points in the \(x\)-direction.

---

## Part (a) – Serial Implementation

- A **serial program** is developed to solve the tridiagonal system using **LU decomposition**.
- The numerical derivative is compared with the analytical solution.
- Results are plotted for \(N = 25\).

---

## Part (b) – OpenMP Implementation

- An **OpenMP program** is implemented using the **recursive-doubling algorithm** to solve the derivative calculation efficiently in parallel.
- Numerical and analytical solutions are plotted for \(N = 100\) using **2 threads**.
- Performance is evaluated by plotting the time taken for \(N = 1000\) using **2, 4, and 8 threads**.

---


## Code and Report

- The source codes for each part are provided in the directory.
- The results, including plots and analysis, are compiled in `report.pdf`.

---
# Poisson Equation

## Overview

This repository contains implementations to solve the following Poisson equation:

$$
\nabla^2 \phi = -\rho; \quad \rho = 2(2 - x^2 - y^2); \quad (\pm 1, y) = 0; \quad (x, \pm 1) = 0
$$

using the **Gauss-Seidel iterative method** with double precision arithmetic.

The exact solution is given by:

$$
\phi(x, y) = (x^2 - 1)(y^2 - 1)
$$

The discretized form using the Gauss-Seidel method is:

$$
\phi_{i,j}^{(k+1)} = \frac{1}{4}\left(\phi_{i+1,j}^{(k)} + \phi_{i-1,j}^{(k+1)} + \phi_{i,j+1}^{(k)} + \phi_{i,j-1}^{(k+1)} + h^2 \rho_{i,j}\right)
$$

where $\Delta x = \Delta y = h$.

---

## Part (a) – Serial Gauss-Seidel Implementation

- A **serial program** is developed to solve the discretized equations using an initial guess of $\phi(x, y) = 0$ everywhere.
- The system is solved for $h = 0.1$ (21 points along each axis).
- The numerical solution is iterated until it is within 1% of the exact solution.
- The numerical and analytical solutions are plotted for $y = 0.5$.

---

## Part (b) – OpenMP Implementation

- Two parallel Gauss-Seidel methods are implemented using OpenMP:
  1. **Diagonal approach**
  2. **Red-black coloring approach**

---

## Part (c) – Verification and Performance Study

- The parallel programs are verified against the serial program for $h = 0.1$ by plotting both solutions together.
- Upon successful verification, the program is run for $h = 0.1$, $0.01$, and $0.005$ using 8 threads.
- The time taken by both serial and parallel solvers is plotted as a function of grid size ($h$).
- The performance improvement and the better parallel method are analyzed.

---

## Part (d) – Thread Scaling Study

- For a grid size of $h = 0.005$, the time taken by each parallel solver is plotted as a function of the number of threads.
- The tests are conducted for $N_{\text{threads}} = 2, 4, 8, 16$.
- The performance and scalability of both parallel methods are compared.

---

## Code and Report

- The source codes for each part are provided in the directory.
- The results, including plots and analysis, are compiled in `report.pdf`.

---


## Contact

For any questions or suggestions, feel free to open an issue or contact the repository maintainer.

## License

This project is open-source and available under the MIT License.

---

## Contact

For any questions or suggestions, please feel free to open an issue or contact the repository maintainer.

