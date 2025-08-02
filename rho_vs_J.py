import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-GUI backend
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize
import math
from compute  import compute_breakup_omega
from derivative import find_slope_savgol



def plot_all_masses1(rho, J,radius, spin, mass,pressure,B,alpha,mass_tolerance=1e-4):
    mass_targets = np.linspace(1.3,1.45,1)
    cmap = matplotlib.colormaps['viridis']
    norm = Normalize(vmin=mass_targets.min(), vmax=mass_targets.max())
    J_mass = {}
    rho_mass = {}
    R_mass = {}
    omega_mass = {}
    mass_dj_dt = {}
    mass_pc = {}
    c = 2.99792458e10  # speed of light in cm/
    sin2_alpha = np.sin(alpha)**2
    fig, ax = plt.subplots(figsize=(8, 6))

    for mass_target in mass_targets:
        group_rho = []
        group_J = []
        group_R = []
        group_omega = []
        group_pc  = []
 
        for i in range(len(mass)):
            if abs(mass[i] - mass_target) <= mass_tolerance:
                if rho[i] > 0 and J[i] > 0:
                    group_rho.append(rho[i])
                    group_J.append(J[i])
                    group_R.append(radius[i])
                    group_omega.append(spin[i])
                    group_pc.append(pressure[i])
        

        if group_rho:
            log_rho = np.log10(np.array(group_rho))
            log_J = np.log10(np.array(group_J))
            #log_omega = np.log10(np.array(group_omega))
            #log_R = np.log10(np.array(group_R))
            sort_idx = np.argsort(log_rho)
            log_J_s = log_J[sort_idx]
            log_rho_s = log_rho[sort_idx]
            group_J_s = np.array(group_J)[sort_idx]
            group_omega_s = np.array(group_omega)[sort_idx]
            group_R_s = np.array(group_R)[sort_idx]
            group_rho_s = np.array(group_rho)[sort_idx]
            group_pc_s = np.array(group_pc)[sort_idx]
            ax.plot(log_J_s, log_rho_s, color=cmap(norm(mass_target)), linewidth=1)
        J_mass[mass_target] = group_J_s
        rho_mass[mass_target] = group_rho_s
        mass_dj_dt[mass_target] = - (2/3) * B**2 * group_R_s**6 * group_omega_s**3 * sin2_alpha / c**3
        mass_pc[mass_target] = group_pc_s
    
    
    mass_slope = find_slope_savgol(rho_mass,J_mass)
    #print(mass_slope)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Mass (M☉)", rotation=270, labelpad=15)
    ax.set_xlabel("log₁₀(J) [g·cm²/s]")
    ax.set_ylabel("log₁₀(ρ$_{c}$) [g/cm$^{3}$]")  # Adjust labels if needed
    ax.set_title("log(ρ_c) vs log(J) for Various Masses")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig("Rho_vs_J_all_masses1e-4.png", dpi=300)
    #plt.show()
    return mass_slope,mass_dj_dt,mass_pc,rho_mass,J_mass
