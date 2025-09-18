# OpenACC — Padé Scheme Derivative Solver

## Overview

This project computes the derivative of the function  

$$
u(x) = \sin(5x), \quad 0 \leq x \leq 3
$$  

using a **fourth-order accurate Padé scheme** for interior points and a **third-order one-sided scheme** near the boundaries.  
The resulting tri-diagonal system is solved using the **Thomas algorithm (LU decomposition)**.

Starting from a **serial implementation**, the solver is extended to **OpenACC**.  
Optimization focuses on reducing host–device data transfer and improving execution time.  
Since the Thomas algorithm itself is inherently sequential, only the parallelizable parts (RHS assembly, result processing) are accelerated.

---

## Governing equations

For an interior point \(i\):  

$$
u'_{i+1} + 4u'_i + u'_{i-1} = \frac{3}{h}\big(u_{i+1} - u_{i-1}\big),
$$  

where \(h\) is the grid spacing.  

Boundary treatment (third-order accurate one-sided Padé scheme):  

$$
u'_0 + 2u'_1 = \frac{1}{h}\left(-\tfrac{5}{2}u_0 + 2u_1 + \tfrac{1}{2}u_2\right),
$$  

$$
u'_N + 2u'_{N-1} = \frac{1}{h}\left(\tfrac{5}{2}u_N - 2u_{N-1} - \tfrac{1}{2}u_{N-2}\right).
$$  

This forms a **tri-diagonal system** solved using the Thomas algorithm.

---

## Implementations

- **Serial version**  
  - Solves the tridiagonal system with the Thomas algorithm.  
  - Plots numerical vs analytical derivative for \(N = 25, 100\).  

- **OpenACC version**  
  - Adds parallelism with OpenACC pragmas for RHS assembly and vectorizable loops.  
  - Reduces unnecessary host–device transfers by keeping arrays resident on GPU.  
  - Thomas algorithm remains serial, but surrounding computations are accelerated.  

---

## Results

- **Accuracy:**  
  Numerical solutions match analytical derivative for \(N = 100\).  

- **Performance:**  
  For \(N = 1000\), the runtime is measured for different numbers of gangs:  
  - 1000 gangs  
  - 100 gangs  
  - 10 gangs  

Plots and detailed performance discussion are in [`report.pdf`](./report.pdf).


