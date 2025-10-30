<div align="center">

# Spectral Petrov–Galerkin Method for Two-Sided Fractional Reaction–Diffusion Equations

**MATLAB Implementation** using **Jacobi Poly-Fractonomials** and **Gauss-Jacobi Quadrature**  
*Chapter 4 – M.Sc. Thesis (October 2025)*

</div>

---

## Author & Contact

| Role | Name | Email |
|------|------|-------|
| **Student** | Hamidreza Karimi | [ha.karimi@sci.ui.ac.ir](mailto:ha.karimi@sci.ui.ac.ir) |
| **Personal** | — | [hamidre87@gmail.com](mailto:hamidre87@gmail.com) |
| **Supervisor** | Dr. Hassan Khosravian Arab | [h.khosravian@sci.ui.ac.ir](mailto:h.khosravian@sci.ui.ac.ir) |

> **We can be reached via email.**

---

## Thesis Title
> **Error estimates of a spectral Petrov–Galerkin method for two-sided fractional reaction–diffusion equations**

**Institution**: Faculty of Mathematical Sciences and Statistics, University of Isfahan, Iran

---

## Repository Contents

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

---

## Method Overview
- **Basis**: Jacobi poly-fractonomials \( P_n^{\alpha,\beta}(x) \)
- **Quadrature**: Gauss-Jacobi (exact for polynomials up to degree 2N+1)
- **Petrov–Galerkin**: Test functions shifted by fractional order
- **Stiffness & Mass Matrices**: Computed via quadrature (no explicit assembly)

> Reference: Shen, J., Tang, T., & Wang, L. L. (2011). *Spectral Methods: Algorithms, Analysis and Applications*.

---

## Requirements
- **MATLAB R2023a or later**
- **Optimization Toolbox** (`fsolve` for σ computation)

---

<div align="center">

**Author**: Hamidreza Karimi  
**Date**: October 2025  
**Location**: Isfahan, Iran  

</div>
