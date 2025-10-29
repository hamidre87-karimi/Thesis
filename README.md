# spectral Petrov–Galerkin method for two-sided fractional reaction–diffusion equations

MATLAB implementation of the **Petrov-Galerkin spectral method** using **Jacobi bases** and **Gauss-Jacobi quadrature**.

Used in **Chapter 4** of my thesis (October 2025).

---

## Files
| File | Description |
|------|-------------|
| `solveEquations.m` | Compute σ and σ* from θ and α |
| `Error_Convergence_Example1.m` | E1 error convergence for f(x) = sin(x) |
| `Error_Convergence_Example2.m` | E1 error for f(x) = \|sin(x)\| |
| `Error_Convergence_Example3.m` | E2 error for singular source |
| `Solution_Example1.m` | Plot u_N(x) for N=128 |
| `Solution_Example2.m` | Plot u_N(x) for N=128 |
| `Solution_Example3.m` | Plot u_N(x) for N=128 |
| `demo_sigma_computation.m` | Generate σ/σ* table |
| `jags.m` | Jacobi-Gauss quadrature points and weights |
| `japoly.m` | Evaluate Jacobi polynomials and derivatives |

## Jacobi Functions (from Shen et al., 2011)
---

## Requirements
- MATLAB R2023a or later
- Optimization Toolbox (`fsolve`)

---

**Author**: Hamidreza Karimi  
**Date**: October 2025  
