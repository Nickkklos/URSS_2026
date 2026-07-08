# Finite Difference Method for the Ornstein--Uhlenbeck Fokker--Planck Equation

This folder contains the derivation of a finite-difference scheme for the
one-dimensional Ornstein--Uhlenbeck Fokker--Planck equation.

## Model equation

The equation considered is

$$
\partial_t \rho
=\kappa \partial_x(x\rho)
+
D\partial_{xx}\rho,
\qquad
x\in[x_L,x_R],
\qquad
t\geq 0.
$$

The associated probability flux is

$$J(x,t)=-\kappa x\rho(x,t)-
D\partial_x\rho(x,t).$$

Zero-flux boundary conditions are imposed at both endpoints:

$$J(x_L,t)=J(x_R,t)=0.$$

## Numerical method

The numerical scheme uses:

- a uniform spatial grid;
- central differences for the spatial derivatives;
- the forward Euler method for time integration;
- ghost points to impose the zero-flux boundary conditions.

After eliminating the ghost values, the fully discrete method can be written as

$$\boldsymbol{\rho}^{n+1}=\left(I+A+B\right)\boldsymbol{\rho}^{n},$$

where $A$ represents the drift contribution and $B$ represents the diffusion
contribution.

## Files

- [`FD_OU.tex`](FD_OU.tex): LaTeX source containing the derivation.
- [`FD_OU.pdf`](FD_OU.pdf): Compiled version of the derivation.
- [`week_1.ipynb`](../week_1.ipynb): Python implementation of the numerical scheme.

## Purpose

The Ornstein--Uhlenbeck equation is used as an initial benchmark because its
probability density is known analytically. This allows the finite-difference
solution to be checked before applying the same workflow to the double-well
Fokker--Planck equation.

The finite-difference density will later be compared with density estimates
obtained from Euler--Maruyama simulations and Monte Carlo sampling.

