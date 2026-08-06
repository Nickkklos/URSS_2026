# URSS 2026: Data-Driven Discovery of Fokker–Planck Equations

This repository contains the code and numerical experiments for my 2026
Undergraduate Research Support Scheme (URSS) project at the University of Warwick,
supervised by Dr Radu Cimpeanu.

## Project Overview

The aim of this project is to investigate whether the governing Fokker–Planck
equation can be identified directly from observed probability density data.

A general one-dimensional Fokker–Planck equation has the form

$$ \frac{\partial \rho}{\partial t}=
-\frac{\partial}{\partial x}\bigl(b(x)\rho\bigr)
+
\frac{1}{2}
\frac{\partial^2}{\partial x^2}
\bigl(\sigma^2(x)\rho\bigr).$$


where:

- $\rho(x,t)$ is the probability density;
- $b(x)$ is the drift coefficient;
- $\sigma(x)$ is the diffusion coefficient.

The main objective is to recover the unknown drift and diffusion terms from
numerically generated or observed density data.

In this project we focus on the ***Ornstein-Uhlenbeck*** process, which is a stochastic differential equation of the following form

$$dX_t = -\kappa X_t \, dt+ \sqrt{2D} \,dW_t$$

and the associated Fokker-Planck equation is

$$\rho_t = \kappa \partial_x(x\rho)+ D \partial_{xx}\rho.$$

The main tasks are:

1. Solve the corresponding Fokker–Planck equation using finite differences.
2. Simulate the associated stochastic differential equation using the
   Euler–Maruyama method.
3. Estimate the probability density from an ensemble of stochastic trajectories.
4. Compare the finite-difference solution with the empirical density.
5. Use the resulting data to test equation-identification methods.
6. How to improve our methods if the result is not satisfied.




## Numerical Methods

The project involves:

- finite-difference methods for the Fokker–Planck equation;
- Euler–Maruyama simulation of stochastic differential equations;
- probability-density estimation from trajectory data;
- sparse regression and data-driven equation discovery;
- error analysis and weak formulation.

