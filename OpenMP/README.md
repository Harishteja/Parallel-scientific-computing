# Numerical Derivative Calculation Using Padé Scheme

## Overview

This repository contains implementations to compute the derivative of the function:

\[
f(x) = \sin(5x), \quad 0 \leq x \leq 3
\]

using high-order accurate Padé schemes. The derivative is computed using:

- **Fourth-order accurate Padé scheme** for interior points
- **Third-order accurate one-sided Padé scheme** near boundaries

### Governing Equations

For interior points (\(i = 1, 2, ..., N-1\)):

\[
f'_{i+1} + 4f'_i + f'_{i-1} = \frac{3}{h} (f_{i+1} - f_{i-1})
\]

For the left boundary:

\[
f'_0 + 2f'_1 = \frac{1}{h}\left(-\frac{5}{2}f_0 + 2f_1 + \frac{1}{2}f_2\right)
\]

For the right boundary:

\[
f'_N + 2f'_{N-1} = \frac{1}{h}\left(\frac{5}{2}f_N - 2f_{N-1} - \frac{1}{2}f_{N-2}\right)
\]

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

## How to Run

1. Clone the repository:
    ```bash
    git clone https://github.com/your-username/parallel-scientific-computing.git
    cd parallel-scientific-computing
    ```

2. Compile the serial version:
    ```bash
    gcc serial_derivative.c -o serial -lm
    ```

3. Compile the OpenMP version:
    ```bash
    gcc openmp_derivative.c -fopenmp -o openmp -lm
    ```

4. Run the serial version:
    ```bash
    ./serial
    ```

5. Run the OpenMP version:
    ```bash
    export OMP_NUM_THREADS=2
    ./openmp
    ```

---

## Results

- Analytical vs numerical derivative plots are available in the `plots/` folder.
- Performance comparison for different thread counts is presented in `performance/`.

---

## Requirements

- GCC with OpenMP support
- Python with `matplotlib` for plotting (optional but recommended)

---

## License

This project is open-source and available under the MIT License.

---

## Contact

For any questions or suggestions, please feel free to open an issue or contact the repository maintainer.

