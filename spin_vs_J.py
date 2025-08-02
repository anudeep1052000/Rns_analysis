import matplotlib
matplotlib.use('Agg')  # Non-GUI backend
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import numpy as np
from matplotlib.colors import Normalize
import math
from matplotlib.colors import Normalize
from compute  import compute_breakup_omega

mass_tolerance = 1e-4
C = 2.99792458e10       # speed of light in cm/s
G = 6.67430e-8          # gravitational constant in cm^3 / g / s^2
seconds_in_year = 365.25 * 24 * 3600 

# Plotting function
def plot_all_masses(rho, J, spin,radius_values, mass):
    mass_targets = np.linspace(1.3, 1.45, 50)
    #cmap = get_cmap('viridis')
    cmap = matplotlib.colormaps['viridis']
    norm = Normalize(vmin=mass_targets.min(), vmax=mass_targets.max())
    fig, ax = plt.subplots(figsize=(8, 6))
    plt.figure(figsize=(8,6))
    omega_limit = []
    omega_max = []
    
    for mass_target in mass_targets:
        group_spin = []
        group_J = []
        group_r = []

        for i in range(len(mass)):
            if abs(mass[i] - mass_target) <= mass_tolerance:
                if spin[i] > 0 and J[i] > 0:
                    group_spin.append(spin[i])
                    group_J.append(J[i])
                    group_r.append(radius_values[i])
        
        if group_spin:
            log_spin = np.log10(np.array(group_spin))
            log_J = np.log10(np.array(group_J))
            log_r = np.log10(np.array(group_r))
            sort_idx = np.argsort(log_J)
            log_J_s = log_J[sort_idx]
            log_spin_s = log_spin[sort_idx]
            log_radius_s = log_r[sort_idx]
            omega_limit.append(compute_breakup_omega(mass_target,10**log_radius_s[-1]))
            omega_max.append(10**log_spin_s[-1])
            


            plt.plot(log_J_s, log_spin_s, color=cmap(norm(mass_target)), linewidth=1)
        
    #print(omega_limit)
    #print(omega_max)
    #print(mass_slope)
    #sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    #sm.set_array([])
    #cbar = plt.colorbar(sm, pad=0.02)
    #cbar.set_label("Mass (M☉)", rotation=270, labelpad=15)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Mass (M☉)", rotation=270, labelpad=15)
    plt.xlabel("log₁₀(J) [g·cm²/s]")
    plt.ylabel("log₁₀(Spin) [rad/s]")
    plt.title("log(Spin) vs log(J) for Various Masses")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("spin_vs_J_all_masses1e-4.png", dpi=300)