import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-GUI backend for headless environments
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import math
from compute import compute_breakup_omega, remove_near_duplicates_sorted
from derivative import find_slope_savgol

def plot_all_masses1(rho, J, radius, spin, mass, pressure, B, alpha, mass_tolerance=1e-4):
    """
    Plot log10(central density) vs. log10(angular momentum) for a set of white dwarf models.

    Parameters:
        rho (list): Central densities.
        J (list): Angular momenta.
        radius (list): Radii.
        spin (list): Spin frequencies (rad/s).
        mass (list): Masses (M_sun).
        pressure (list): Central pressures.
        B (float): Magnetic field strength.
        alpha (float): Magnetic axis inclination.
        mass_tolerance (float): Allowed tolerance when selecting models by mass.

    Returns:
        tuple: (mass_slope, mass_dj_dt, mass_pc, rho_mass, J_mass)
    """
    # Define mass targets to analyze
    mass_targets = np.linspace(1.3, 1.45, 1)
    cmap = matplotlib.colormaps['viridis']
    norm = Normalize(vmin=mass_targets.min(), vmax=mass_targets.max())

    # Containers for grouped results
    J_mass = {}
    rho_mass = {}
    R_mass = {}
    omega_mass = {}
    mass_dj_dt = {}
    mass_pc = {}

    c = 2.99792458e10  # Speed of light in cm/s
    sin2_alpha = np.sin(alpha)**2

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    for mass_target in mass_targets:
        # Temporary containers for each mass_target
        group_rho = []
        group_J = []
        group_R = []
        group_omega = []
        group_pc = []

        # Filter data points matching the target mass
        for i in range(len(mass)):
            if abs(mass[i] - mass_target) <= mass_tolerance:
                if rho[i] > 0 and J[i] > 0:
                    group_rho.append(rho[i])
                    group_J.append(J[i])
                    group_R.append(radius[i])
                    group_omega.append(spin[i])
                    group_pc.append(pressure[i])

        if group_rho:
            # Convert to arrays and sort by central density
            log_rho = np.log10(np.array(group_rho))
            log_J = np.log10(np.array(group_J))
            sort_idx = np.argsort(log_rho)

            # Sort arrays
            group_rho_s = np.array(group_rho)[sort_idx]
            group_J_s = np.array(group_J)[sort_idx]
            group_R_s = np.array(group_R)[sort_idx]
            group_omega_s = np.array(group_omega)[sort_idx]
            group_pc_s = np.array(group_pc)[sort_idx]

            # Remove near-duplicate points for smooth curves
            group_rho_s, group_J_s, group_R_s, group_omega_s, group_pc_s = remove_near_duplicates_sorted(
                group_rho_s, group_J_s, group_R_s, group_omega_s, group_pc_s
            )

            # Plot log10(J) vs. log10(rho)
            ax.plot(np.log10(group_J_s), np.log10(group_rho_s),
                    color=cmap(norm(mass_target)), linewidth=1)

            # Store processed data for this mass_target
            J_mass[mass_target] = group_J_s
            rho_mass[mass_target] = group_rho_s
            mass_dj_dt[mass_target] = - (2/3) * B**2 * group_R_s**6 * group_omega_s**3 * sin2_alpha / c**3
            mass_pc[mass_target] = group_pc_s

    # Compute slope d(rho)/dJ using Savitzky-Golay smoothing
    mass_slope = find_slope_savgol(rho_mass, J_mass)

    # Add colorbar and labels
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Mass (M☉)", rotation=270, labelpad=15)

    ax.set_xlabel("log₁₀(J) [g·cm²/s]")
    ax.set_ylabel("log₁₀(ρ$_{c}$) [g/cm$^{3}$]")
    ax.set_title("log(ρ_c) vs. log(J) for Various Masses")
    ax.grid(True)

    fig.tight_layout()
    fig.savefig("Rho_vs_J_all_masses1e-4.png", dpi=300)

    return mass_slope, mass_dj_dt, mass_pc, rho_mass, J_mass
