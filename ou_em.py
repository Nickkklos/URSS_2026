import numpy as np


def solve_ou_em(
    num_paths,
    N_em,
    Nx,
    *,
    k,
    D,
    x_L,
    x_R,
    m0,
    s0,
    T,
    save_times=None,
    seed=319,
    store_particles=False,
):
    """
    Simulate the Ornstein-Uhlenbeck SDE using Euler-Maruyama,

        dX_t = -k X_t dt + sqrt(2D) dW_t,

    and construct histogram density snapshots.

    The main output conventions are

        result["snapshots"][time]
            Histogram density at the requested time.

        result["particle_snapshots"][time]
            Particle ensemble at the requested time, provided only when
            store_particles=True.

    Parameters
    ----------
    num_paths : int
        Number of Monte Carlo paths.
    N_em : int
        Number of Euler-Maruyama time steps on [0, T].
    Nx : int
        Number of spatial intervals. The density grid has Nx + 1 points.
    k, D : float
        Ornstein-Uhlenbeck drift and diffusion parameters.
    x_L, x_R : float
        Histogram domain.
    m0, s0 : float
        Mean and standard deviation of the Gaussian initial distribution.
    T : float
        Final time.
    save_times : array-like or None
        Requested snapshot times. Every time must lie on the EM time grid.
    seed : int or None
        Random seed.
    store_particles : bool
        If True, save particle ensembles at requested times. This is useful
        for moment estimation and Monte Carlo convergence tests, but should
        normally be False when many times are saved for SINDy identification.

    Returns
    -------
    dict
        Dictionary containing the spatial grid, histogram densities,
        optional particle snapshots and diagnostics.
    """

    # ============================================================
    # Basic checks
    # ============================================================

    if not isinstance(num_paths, (int, np.integer)) or num_paths < 1:
        raise ValueError("num_paths must be a positive integer.")

    if not isinstance(N_em, (int, np.integer)) or N_em < 1:
        raise ValueError("N_em must be a positive integer.")

    if not isinstance(Nx, (int, np.integer)) or Nx < 2:
        raise ValueError("Nx must be an integer at least 2.")

    if not isinstance(store_particles, (bool, np.bool_)):
        raise TypeError("store_particles must be True or False.")

    if x_L >= x_R:
        raise ValueError("x_L must be smaller than x_R.")

    if s0 <= 0.0:
        raise ValueError("s0 must be positive.")

    if D < 0.0:
        raise ValueError("D must be non-negative.")

    if T <= 0.0:
        raise ValueError("T must be positive.")

    # ============================================================
    # Time grid
    # ============================================================

    dt_em = T / N_em

    # ============================================================
    # Spatial density grid and histogram bins
    # ============================================================

    x = np.linspace(
        x_L,
        x_R,
        Nx + 1,
    )

    dx = x[1] - x[0]

    # One histogram bin is centred at each x-grid point.
    bin_edges = np.linspace(
        x_L - 0.5 * dx,
        x_R + 0.5 * dx,
        Nx + 2,
    )

    bin_widths = np.diff(bin_edges)

    bin_centres = 0.5 * (
        bin_edges[:-1] + bin_edges[1:]
    )

    # ============================================================
    # Determine the requested EM save steps
    # ============================================================

    if save_times is None:
        save_times = np.empty(0, dtype=float)
    else:
        save_times = np.asarray(
            save_times,
            dtype=float,
        ).reshape(-1)

    save_steps = {}

    for requested_time in save_times:

        if not np.isfinite(requested_time):
            raise ValueError("All save_times must be finite.")

        if requested_time < 0.0 or requested_time > T:
            raise ValueError(
                f"Save time {requested_time} is outside [0, T]."
            )

        step = int(
            np.rint(requested_time / dt_em)
        )

        actual_time = step * dt_em

        if not np.isclose(
            requested_time,
            actual_time,
            rtol=0.0,
            atol=1e-12 * max(1.0, T),
        ):
            raise ValueError(
                f"Save time {requested_time} is not on the EM time grid. "
                f"The nearest grid time is {actual_time}."
            )

        if step in save_steps:
            raise ValueError(
                f"Multiple save times correspond to EM step {step}."
            )

        # Use the actual grid time as the dictionary key.
        save_steps[step] = float(actual_time)

    # ============================================================
    # Initial particle ensemble
    # ============================================================

    rng = np.random.default_rng(seed)

    X_current = rng.normal(
        loc=m0,
        scale=s0,
        size=num_paths,
    )

    # ============================================================
    # Convert particles into a histogram density
    # ============================================================

    def histogram_density(particles):
        counts, _ = np.histogram(
            particles,
            bins=bin_edges,
        )

        density = counts / (
            num_paths * bin_widths
        )

        mass_inside = (
            counts.sum() / num_paths
        )

        return density, mass_inside

    # Density snapshots are always stored because they are the data
    # required for plotting and Fokker-Planck identification.
    density_snapshots = {}

    # Particle snapshots are stored only when explicitly requested.
    particle_snapshots = {} if store_particles else None

    inside_mass = {}

    # ============================================================
    # Initial density and optional initial particle snapshot
    # ============================================================

    rho_initial, initial_inside_mass = histogram_density(
        X_current
    )

    if 0 in save_steps:
        time = save_steps[0]

        density_snapshots[time] = rho_initial.copy()
        inside_mass[time] = initial_inside_mass

        if store_particles:
            particle_snapshots[time] = X_current.copy()

    # ============================================================
    # Euler-Maruyama time stepping
    # ============================================================

    noise_scale = np.sqrt(
        2.0 * D * dt_em
    )

    rho_final = None
    final_inside_mass = None
    X_final = None

    for step in range(1, N_em + 1):

        X_current = (
            X_current
            - k * X_current * dt_em
            + noise_scale
            * rng.standard_normal(num_paths)
        )

        # Construct a histogram only at requested times and at T.
        if step in save_steps or step == N_em:

            rho_current, mass_current = histogram_density(
                X_current
            )

            if step in save_steps:
                time = save_steps[step]

                density_snapshots[time] = rho_current.copy()
                inside_mass[time] = mass_current

                if store_particles:
                    particle_snapshots[time] = X_current.copy()

            if step == N_em:
                rho_final = rho_current.copy()
                final_inside_mass = mass_current
                X_final = X_current.copy() if store_particles else None

    snapshot_times = np.array(
        sorted(density_snapshots.keys()),
        dtype=float,
    )

    return {
        "x": x,
        "rho_initial": rho_initial,
        "rho_final": rho_final,
        "dx": dx,
        "dt_em": dt_em,
        "bin_edges": bin_edges,
        "bin_widths": bin_widths,
        "bin_centres": bin_centres,

        # Density data:
        # snapshots[time] has shape (Nx + 1,).
        "snapshots": density_snapshots,
        
        "times": snapshot_times,

        # Particle data:
        # particle_snapshots[time] has shape (num_paths,) when enabled.
        "particle_snapshots": particle_snapshots,
        "particles_final": X_final,

        # Histogram-domain diagnostics.
        "inside_mass": inside_mass,
        "initial_inside_mass": initial_inside_mass,
        "final_inside_mass": final_inside_mass,
    }
