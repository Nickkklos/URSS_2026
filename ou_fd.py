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


    # Basic checks
    if Nx < 2:
        raise ValueError("Nx must be at least 2.")

    if Nt < 1:
        raise ValueError("Nt must be a positive integer.")


    # Spatial grid
    x = np.linspace(x_L, x_R, Nx + 1)
    dx = x[1] - x[0]

   
    # Temporal grid
    dt = T / Nt

    
    # Initial Gaussian density
    rho_initial = (
        np.exp(-0.5 * ((x - m0) / s0)**2)
        / (np.sqrt(2.0 * np.pi) * s0)
    )

    # rho will be updated during time stepping
    rho = rho_initial.copy()

    # Sparse matrices
    A = lil_matrix((Nx + 1, Nx + 1))
    B = lil_matrix((Nx + 1, Nx + 1))


    # Forward-Euler contribution from
    # k * d(x rho)/dx

    for j in range(Nx + 1):

        # Left boundary
        if j == 0:
            A[0, 0] = ( -k**2 * dt / D* (x_L - dx)* x_L)
            A[0, 1] = k * dt

        # Right boundary
        elif j == Nx:
            A[j, j] = (-k**2 * dt / D* x_R* (x_R + dx))
            A[j, j - 1] = k * dt

        # Interior points
        else:
            A[j, j - 1] = ( -k * dt / (2.0 * dx) * x[j - 1])
            A[j, j + 1] = ( k * dt / (2.0 * dx)* x[j + 1])


    # Discretisation of the second derivative
    for j in range(Nx + 1):

        # Left boundary
        if j == 0:
            B[0, 0] = (2.0 * k * dx / D * x_L - 2.0)
            B[0, 1] = 2.0

        # Right boundary
        elif j == Nx:
            B[j, j] = (-2.0 * k * dx / D * x_R - 2.0 )
            B[j, j - 1] = 2.0

        # Interior points
        else:
            B[j, j - 1] = 1.0
            B[j, j] = -2.0
            B[j, j + 1] = 1.0


    # Assemble the Forward Euler update matrix
    # rho^{n+1} = M rho^n
    M = (eye(Nx + 1, format="csr")+ A.tocsr() + D * dt / dx**2 * B.tocsr())


    
    # Determine the time steps at which snapshots are saved
    if save_times is None:
        save_times = []

    save_steps = {}

    for requested_time in save_times:

        # Check whether the requested time lies in [0, T]
        if requested_time < 0.0 or requested_time > T:
            raise ValueError(
                f"Save time {requested_time} is outside [0, T]."
            )

        # Find the corresponding temporal grid point
        step = int(np.rint(requested_time / dt))
        actual_time = step * dt

        # Require the requested time to lie on the time grid
        if not np.isclose(requested_time, actual_time, rtol=0.0, atol=1e-12 * max(1.0, T)):
            raise ValueError(
                f"Save time {requested_time} is not on the temporal grid. "
                f"The nearest grid time is {actual_time}."
            )

        # Avoid assigning two requested times to the same step
        if step in save_steps:
            raise ValueError(
                f"Multiple save times correspond to step {step}."
            )

        # Store only after all checks have passed
        save_steps[step] = actual_time

    snapshots = {}

    # Save the initial condition separately because the
    # time-stepping loop starts from step 1
    if 0 in save_steps:
        snapshots[save_steps[0]] = rho.copy()

    # Diagnostics
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

    # Forward Euler time stepping
    for step in range(1, Nt + 1):

        rho = M @ rho

        # Save diagnostics at every time step
        if track_history:
            mass_history[step] = np.trapezoid(rho, x)
            minimum_history[step] = np.min(rho)

        # Save selected snapshots independently of track_history
        if step in save_steps:
            snapshots[save_steps[step]] = rho.copy()


    # Return results   
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