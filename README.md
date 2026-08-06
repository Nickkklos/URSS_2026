# URSS 2026: Data-Driven Discovery of Fokker–Planck Equations

This repository contains the code, mathematical derivations, and numerical experiments for my 2026 Undergraduate Research Support Scheme (URSS) project at the University of Warwick, supervised by Dr Radu Cimpeanu.

The project investigates whether the governing Fokker–Planck equation can be recovered directly from probability-density data using sparse equation-identification methods.

## Research Question

A general one-dimensional Fokker–Planck equation has the form

$$
\frac{\partial \rho}{\partial t}
================================

-\frac{\partial}{\partial x}\bigl(b(x)\rho\bigr)
+
\frac{1}{2}
\frac{\partial^2}{\partial x^2}\bigl(\sigma^2(x)\rho\bigr),
$$

where

* $\rho(x,t)$ is the probability density;
* $b(x)$ is the drift coefficient;
* $\sigma(x)$ is the diffusion coefficient.

The main question is whether the drift and diffusion terms can be identified from sampled values of $\rho(x,t)$, and how the reliability of the identification depends on the data source, candidate library, and regression method.

## Benchmark Problem: Ornstein–Uhlenbeck Process

The current identification experiments focus on the Ornstein–Uhlenbeck process

$$
dX_t=-\kappa X_t,dt+\sqrt{2D},dW_t,
$$

whose probability density satisfies

$$
\rho_t
======

\kappa,\partial_x(x\rho)
+
D,\partial_{xx}\rho.
$$

Because the exact density is known analytically, this model provides a controlled benchmark for comparing the recovered coefficients with the true values $\kappa$ and $D$.

On a truncated computational domain, the finite-difference solver uses zero-flux boundary conditions.

## Methodology

Probability-density data are generated in three ways:

1. **Exact density:** analytical Ornstein–Uhlenbeck density;
2. **Finite-difference density:** numerical solution of the Fokker–Planck equation;
3. **Euler–Maruyama density:** histogram estimates constructed from simulated stochastic trajectories.

Two equation-identification approaches are then studied:

* **Pointwise SINDy**, in which the required derivatives are approximated numerically before sparse regression;
* **Weak SINDy**, in which the equation is integrated against local test functions, transferring derivatives from the noisy density onto smooth test functions.

The recovered equations are assessed using their selected terms, coefficient errors, regression residuals, and agreement with the exact time derivative.

## Repository Structure

| File or folder                                     | Description                                                                                                                              |
| -------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `ou_fd.py`                                         | Reusable finite-difference solver for the Ornstein–Uhlenbeck Fokker–Planck equation                                                      |
| `ou_em.py`                                         | Euler–Maruyama simulator and histogram-density constructor                                                                               |
| `week_1.ipynb`                                     | Mathematical background, finite-difference and Euler–Maruyama methods, convergence tests, and preliminary double-well experiments        |
| `week_1_external_solvers.ipynb`                    | Refactored version of the Week 1 experiments using `ou_fd.py` and `ou_em.py`                                                             |
| `week_5.ipynb`                                     | Combined exploratory workflow for exact, finite-difference, and Euler–Maruyama density identification                                    |
| `week_6_Exact.ipynb`                               | Pointwise SINDy identification using the exact Ornstein–Uhlenbeck density                                                                |
| `week_6_FD.ipynb`                                  | Pointwise SINDy identification using finite-difference density data and different candidate libraries                                    |
| `week_7_em.ipynb`                                  | Weak-SINDy identification from Euler–Maruyama histogram data, including comparisons over candidate-library size, $K$, and subdomain size |
| `FD_OU/`                                           | Derivation of the finite-difference scheme and zero-flux boundary treatment                                                              |
| `Stochastic_process_and_FPE_background_knowledge/` | Background notes on stochastic processes and Fokker–Planck equations                                                                     |

## Running the Notebooks

The experiments were developed using Python 3.11.

Clone the repository and install the main dependencies:

```bash
git clone https://github.com/Nickkklos/URSS_2026.git
cd URSS_2026

python -m pip install numpy scipy pandas matplotlib pysindy jupyterlab
jupyter lab
```

A suggested reading order is:

1. `week_1_external_solvers.ipynb`;
2. `week_6_Exact.ipynb`;
3. `week_6_FD.ipynb`;
4. `week_7_em.ipynb`.

## Preliminary Findings

The experiments show that Weak-SINDy performance depends jointly on

* the number $K$ of sampled weak equations;
* the size of the space–time subdomains;
* the number of terms in the candidate library.

Small subdomains may not sufficiently average the fluctuations in Euler–Maruyama histogram data. Very large subdomains suppress noise more strongly, but may also make different candidate terms difficult to distinguish. Consequently, a small weak residual does not necessarily imply that the correct PDE has been identified.

In the current single-seed parameter sweep, a subdomain full width equal to $25%$ of each complete space–time axis is the most robust common choice for both candidate libraries. Increasing $K$ is not automatically beneficial and cannot fully compensate for an unsuitable subdomain size.

These numerical ranges are specific to the current data resolution, regression threshold, and random seed, and will be tested further using repeated samples.

## Current Scope

The Ornstein–Uhlenbeck equation is used as the main identification benchmark. Preliminary numerical experiments for a double-well model are also included, with extension of the full identification workflow to more complex drift fields forming a natural next step.
