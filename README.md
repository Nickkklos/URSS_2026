# URSS 2026: Data-Driven Discovery of Fokker-Planck Equations

This repository contains the mathematical derivations, numerical solvers, and equation-discovery experiments from my 2026 Undergraduate Research Support Scheme (URSS) project at the University of Warwick, supervised by Dr Radu Cimpeanu.

The project asks whether a Fokker-Planck equation can be recovered directly from probability-density data by sparse regression, and how reliable that recovery remains when the density is estimated from stochastic trajectories and affected by realistic observation errors.

## Research question

For the one-dimensional stochastic differential equation

$$
dX_t=b(X_t)\,dt+\sigma(X_t)\,dW_t,
$$

the probability density $\rho(x,t)$ satisfies

$$
\frac{\partial \rho}{\partial t}
=-\frac{\partial}{\partial x}\bigl(b(x)\rho\bigr)
+\frac{1}{2}\frac{\partial^2}{\partial x^2}
\bigl(\sigma^2(x)\rho\bigr).
$$

The central question is:

> Can sparse equation identification recover the active drift and diffusion terms from sampled values of $\rho(x,t)$, and how do data generation, numerical differentiation, candidate-library design, weak-form parameters, and camera errors affect the result?

## Benchmark: Ornstein-Uhlenbeck process

The main benchmark is

$$
dX_t=-\kappa X_t\,dt+\sqrt{2D}\,dW_t,
$$

whose density satisfies the conservative Fokker-Planck equation

$$
\rho_t=\kappa\,\partial_x(x\rho)+D\rho_{xx}.
$$

The experiments use

$$
\kappa=D=1,
\qquad
X_0\sim\mathcal N(-2,0.5^2),
\qquad
(x,t)\in[-5,5]\times[0,3].
$$

Because the OU density and its time derivative are available analytically, the model separates errors caused by data generation, numerical differentiation, library design, and sparse regression. The finite-difference solver uses a zero-flux boundary condition on the truncated spatial domain.

## Data and identification methods

Three density sources are compared:

1. **Exact density** - analytical OU density;
2. **Finite-difference density** - numerical solution of the Fokker-Planck equation;
3. **Euler-Maruyama density** - histograms constructed from 60,000 simulated trajectories.

Two sparse-identification approaches are studied:

- **Pointwise SINDy:** numerical derivatives are computed before sparse regression;
- **Weak SINDy:** the PDE is integrated against compactly supported test functions, and integration by parts transfers derivatives from the noisy density to the smooth test functions.

The recovered model is assessed using support recovery, identified coefficients, relative coefficient error, weak regression residual, and agreement with the exact time derivative when that reference is available.

## Main findings

- Pointwise SINDy accurately recovers the OU equation from smooth exact and finite-difference densities.
- Pointwise differentiation strongly amplifies bin-to-bin fluctuations in Euler-Maruyama histograms. Weak integration restores the correct two-term support for a suitably chosen weak domain and candidate library.
- The number of weak equations $K$ and the support width of the test functions must be selected jointly. A small residual alone does not guarantee the correct PDE.
- Algebraically dependent candidate terms can make individual coefficients non-identifiable even when the combined expression is correct. Compact conservative libraries are preferable when the structure is known.
- In the Week 9 test with $\Delta t_{\rm cam}=0.04$, $\widetilde\sigma=0.5$, $K=500$, and 25% requested full widths, Weak SINDy selects the correct OU operator support for clean, motion-blurred, localisation-only, and combined observations in both the 5- and 15-candidate libraries.
- Camera errors do not leave the observed density unchanged. In the linear-Gaussian OU setting they preserve the two-term OU structure but alter the effective diffusion coefficient. The Week 9 coefficients should therefore be compared with a camera-aware reference as well as with the underlying physical value $D=1$.

These results demonstrate robustness within the tested regime; they are not yet a universal camera-error robustness claim. Multi-seed and multi-noise-level sweeps remain important follow-up tests.

## Repository structure

| File or folder | Purpose |
| --- | --- |
| `ou_fd.py` | Reusable finite-difference solver for the OU Fokker-Planck equation |
| `ou_em.py` | Reusable Euler-Maruyama simulator and histogram-density constructor |
| `week_1_data.ipynb`  | Original self-contained background, solvers, and convergence experiments |
| `week_5.ipynb` | Exploratory exact, finite-difference, Euler-Maruyama, and weak-form workflow |
| `week_6_Exact.ipynb` | Pointwise SINDy on the exact OU density |
| `week_6_FD.ipynb` | Pointwise SINDy on finite-difference density with several candidate libraries |
| `week_7_em.ipynb` | Weak-SINDy identification from Euler-Maruyama histograms; library, $K$, and subdomain-width comparisons |
| `week_9_OU_camera_errors.ipynb` | Trajectory-level motion blur and localisation noise followed by Weak-SINDy identification |
| `FD_OU/` | Derivation of the finite-difference scheme and zero-flux boundary treatment |
| `Stochastic_process_and_FPE_background_knowledge/` | Background notes on stochastic processes and Fokker-Planck equations |
| `poster/` | Editable LaTeX poster source and figure assets |

## Suggested reading order

1. `week_1_data.ipynb`
2. `week_5.ipynb`
3. `week_6_Exact.ipynb`
4. `week_6_FD.ipynb`
5. `week_7_em.ipynb`
6. `week_9_OU_camera_errors.ipynb`

## Running the notebooks

The notebooks were developed with Python 3.11. A minimal environment requires NumPy, SciPy, pandas, Matplotlib, PySINDy, JupyterLab, and IPython.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy scipy pandas matplotlib pysindy jupyterlab ipython
jupyter lab
```

Run notebooks from the repository root so that the local modules `ou_fd.py` and `ou_em.py` can be imported. For a reproducibility check, restart the kernel and use **Run All** rather than relying on variables left from an earlier session.

The stochastic and weak-subdomain experiments use explicit random seeds. A single seed gives a reproducible realisation, not an uncertainty estimate; repeated-seed statistics should be used for stronger claims.



## References

1. S. L. Brunton, J. L. Proctor, and J. N. Kutz, "Discovering governing equations from data by sparse identification of nonlinear dynamical systems," *Proceedings of the National Academy of Sciences*, 113(15), 3932-3937, 2016. [doi:10.1073/pnas.1517384113](https://doi.org/10.1073/pnas.1517384113)
2. S. H. Rudy, S. L. Brunton, J. L. Proctor, and J. N. Kutz, "Data-driven discovery of partial differential equations," *Science Advances*, 3(4), e1602614, 2017. [doi:10.1126/sciadv.1602614](https://doi.org/10.1126/sciadv.1602614)
3. D. A. Messenger and D. M. Bortz, "Weak SINDy for partial differential equations," *Journal of Computational Physics*, 443, 110525, 2021. [doi:10.1016/j.jcp.2021.110525](https://doi.org/10.1016/j.jcp.2021.110525)
4. A. J. Berglund, "Statistics of camera-based single-particle tracking," *Physical Review E*, 82, 011917, 2010. [doi:10.1103/PhysRevE.82.011917](https://doi.org/10.1103/PhysRevE.82.011917)

## Acknowledgements

This project was supported by the University of Warwick Undergraduate Research Support Scheme (URSS) and supervised by Dr Radu Cimpeanu at the Warwick Mathematics Institute.
