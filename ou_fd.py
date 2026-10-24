"""
Finite-difference solver for the Ornstein-Uhlenbeck Fokker-Planck equation.

The density rho(x, t) of dX_t = -k X_t dt + sqrt(2D) dW_t satisfies

    rho_t = k (x rho)_x + D rho_xx,    x in [x_L, x_R],  t in [0, T],

and the truncated domain uses the zero-flux boundary condition

    k x rho + D rho_x = 0    at x = x_L and x = x_R.

The derivation of the scheme is in FD_OU/FD_OU.pdf.
"""

import numpy as np
from scipy.sparse import lil_matrix, eye


def solve_ou_fd(
    Nx,
    Nt,
    *,
    k,
    D,
    x_L,
    x_R,
    m0,
    s0,
    T,
    save_times=None,
    track_history=False,
):
    """
    Solve the OU Fokker-Planck equation with central differences in space
    and forward Euler in time.

    Spatial derivatives use second-order central differences on a uniform
    grid. The zero-flux condition is imposed with one ghost point at each
    boundary; eliminating the ghost values gives the update

        rho^{n+1} = M rho^n,    M = I + A + (D dt / dx**2) B,

    where A holds the drift contribution and B the second-difference
    stencil, both with boundary rows modified by the zero-flux condition.

    The scheme is explicit. A frozen-coefficient von Neumann analysis
    requires D * dt / dx**2 <= 1/2 and (k * max(|x_L|, |x_R|))**2 * dt <= 2 * D.
    These conditions are not checked here.

    Parameters
    ----------
    Nx : int
        Number of spatial intervals. The grid has Nx + 1 points, including
        both boundary points.
    Nt : int
        Number of forward Euler time steps on [0, T].
    k, D : float
        Drift rate and diffusion coefficient of the OU process
        (k is kappa in the README). D must be positive.
    x_L, x_R : float
        Endpoints of the truncated spatial domain.
    m0, s0 : float
        Mean and standard deviation of the Gaussian initial density.
    T : float
        Final time.
    save_times : sequence of float or None
        Times at which to store snapshots of rho. Every time must lie in
        [0, T] and on the time grid, i.e. be a multiple of dt = T / Nt.
    track_history : bool
        If True, record the total mass (trapezoidal rule) and the minimum
        of rho at every time step.

    Returns
    -------
    dict
        "x"                 spatial grid, shape (Nx + 1,)
        "rho_initial"       initial density on the grid
        "rho_final"         density at time T
        "dx", "dt"          grid spacing and time step
        "M"                 sparse CSR update matrix, rho^{n+1} = M rho^n
        "snapshots"         dict mapping each saved grid time to rho
        "times"             time grid, shape (Nt + 1,), or None
        "mass_history"      total mass at each time step, or None
        "minimum_history"   minimum of rho at each time step, or None

        The last three entries are None unless track_history is True.
    """

    # ------------------------------------------------------------
    # Basic checks
    # ------------------------------------------------------------

    if Nx < 2:
        raise ValueError("Nx must be at least 2.")

    if Nt < 1:
        raise ValueError("Nt must be a positive integer.")

    # ------------------------------------------------------------
    # Spatial and temporal grids
    # ------------------------------------------------------------

    x = np.linspace(x_L, x_R, Nx + 1)
    dx = x[1] - x[0]

    dt = T / Nt

    # Initial Gaussian density
    rho_initial = (
        np.exp(-0.5 * ((x - m0) / s0)**2)
        / (np.sqrt(2.0 * np.pi) * s0)
    )

    # rho is updated during time stepping
    rho = rho_initial.copy()

    # ------------------------------------------------------------
    # Drift matrix A: forward Euler contribution of k * d(x rho)/dx
    # ------------------------------------------------------------

    A = lil_matrix((Nx + 1, Nx + 1))
    B = lil_matrix((Nx + 1, Nx + 1))

    for j in range(Nx + 1):

        # Left boundary (ghost value eliminated by the zero-flux condition)
        if j == 0:
            A[0, 0] = -k**2 * dt / D * (x_L - dx) * x_L
            A[0, 1] = k * dt

        # Right boundary
        elif j == Nx:
            A[j, j] = -k**2 * dt / D * x_R * (x_R + dx)
            A[j, j - 1] = k * dt

        # Interior points
        else:
            A[j, j - 1] = -k * dt / (2.0 * dx) * x[j - 1]
            A[j, j + 1] = k * dt / (2.0 * dx) * x[j + 1]

    # ------------------------------------------------------------
    # Diffusion matrix B: second-difference stencil
    # ------------------------------------------------------------

    for j in range(Nx + 1):

        # Left boundary
        if j == 0:
            B[0, 0] = 2.0 * k * dx / D * x_L - 2.0
            B[0, 1] = 2.0

        # Right boundary
        elif j == Nx:
            B[j, j] = -2.0 * k * dx / D * x_R - 2.0
            B[j, j - 1] = 2.0

        # Interior points
        else:
            B[j, j - 1] = 1.0
            B[j, j] = -2.0
            B[j, j + 1] = 1.0

    # Forward Euler update matrix: rho^{n+1} = M rho^n
    M = eye(Nx + 1, format="csr") + A.tocsr() + D * dt / dx**2 * B.tocsr()

    # ------------------------------------------------------------
    # Time steps at which snapshots are saved
    # ------------------------------------------------------------

    if save_times is None:
        save_times = []

    save_steps = {}

    for requested_time in save_times:

        # The requested time must lie in [0, T]
        if requested_time < 0.0 or requested_time > T:
            raise ValueError(
                f"Save time {requested_time} is outside [0, T]."
            )

        # Nearest temporal grid point
        step = int(np.rint(requested_time / dt))
        actual_time = step * dt

        # The requested time must lie on the time grid
        if not np.isclose(
            requested_time,
            actual_time,
            rtol=0.0,
            atol=1e-12 * max(1.0, T),
        ):
            raise ValueError(
                f"Save time {requested_time} is not on the temporal grid. "
                f"The nearest grid time is {actual_time}."
            )

        # Two requested times must not map to the same step
        if step in save_steps:
            raise ValueError(
                f"Multiple save times correspond to step {step}."
            )

        # Store only after all checks have passed
        save_steps[step] = actual_time

    snapshots = {}

    # The time-stepping loop starts from step 1, so the initial
    # condition is saved separately
    if 0 in save_steps:
        snapshots[save_steps[0]] = rho.copy()

    # ------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------

    if track_history:

        times = np.arange(Nt + 1) * dt

        mass_history = np.empty(Nt + 1)
        minimum_history = np.empty(Nt + 1)

        mass_history[0] = np.trapezoid(rho, x)
        minimum_history[0] = np.min(rho)

    else:

        times = None
        mass_history = None
        minimum_history = None

    # ------------------------------------------------------------
    # Forward Euler time stepping
    # ------------------------------------------------------------

    for step in range(1, Nt + 1):

        rho = M @ rho

        # Diagnostics at every time step
        if track_history:
            mass_history[step] = np.trapezoid(rho, x)
            minimum_history[step] = np.min(rho)

        # Snapshots are saved independently of track_history
        if step in save_steps:
            snapshots[save_steps[step]] = rho.copy()

    return {
        "x": x,
        "rho_initial": rho_initial,
        "rho_final": rho,
        "dx": dx,
        "dt": dt,
        "M": M,
        "snapshots": snapshots,
        "times": times,
        "mass_history": mass_history,
        "minimum_history": minimum_history
    }
