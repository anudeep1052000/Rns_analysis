import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-GUI backend for saving figures
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize
from scipy.signal import savgol_filter
M_sun = 1.989e33
mass_tolerance = 1e-5

# Data holders
mass_values = []
spin_values = []
momentum_values = []
rho_values = []
spin_mass = {}
J_mass = {}


def find_slope_savgol(spin_mass, J_mass, window_length=5, polyorder=2):
    mass_slope = {}
    for i in spin_mass.keys():
        a = np.array(spin_mass[i])
        b = np.array(J_mass[i])
        # Sort by independent variable
        sorted_idx = np.argsort(b)
        a_sorted = a[sorted_idx]
        b_sorted = b[sorted_idx]
        # Compute derivative using Savitzky-Golay
        if len(a_sorted) >= window_length:
            slope = savgol_filter(a_sorted, window_length, polyorder, deriv=1, delta=np.mean(np.diff(b_sorted)))
            mass_slope[i] = slope

    return mass_slope

# Detect turning points (maxima/minima)
def inflection_point(mass_slope):
    mass_flat = []
    mass_spinup = []
    mass_spindown = []

    for mass_val, slopes in mass_slope.items():
        threshold = 1e-4
        turning = None

        for i in range(1, len(slopes)):
            prev_s = slopes[i - 1]
            curr_s = slopes[i]
            if prev_s > threshold and curr_s < -threshold:
                turning = 'maximum'
                break
            elif prev_s < -threshold and curr_s > threshold:
                turning = 'minimum'
                break

        if turning == 'maximum':
            mass_spindown.append(mass_val)
        elif turning == 'minimum':
            mass_spinup.append(mass_val)
        else:
            if all(abs(s) < 1e-5 for s in slopes):
                mass_flat.append(mass_val)

    return mass_flat, mass_spinup, mass_spindown

# Compute slopes d(log_spin)/d(log_J) for each mass track
def find_slope(spin_mass, J_mass):
    mass_slope = {}
    for mass_val in spin_mass:
        log_spin = spin_mass[mass_val]
        log_J = J_mass[mass_val]
        slope = []
        for j in range(1, len(log_spin)):
            dy = log_spin[j] - log_spin[j - 1]
            dx = log_J[j] - log_J[j - 1]
            slope.append(dy / dx if dx != 0 else 0.0)
        mass_slope[mass_val] = slope
    return mass_slope

# Plot log(Spin) vs log(J) for all mass bins
def plot_all_masses(rho, J, spin, mass):
    mass_targets = np.linspace(0.5, 1.49, 50)
    cmap = matplotlib.colormaps['viridis']
    norm = Normalize(vmin=mass_targets.min(), vmax=mass_targets.max())
    fig, ax = plt.subplots(figsize=(8, 6))

    for mass_target in mass_targets:
        group_spin = []
        group_J = []

        for i in range(len(mass)):
            if abs(mass[i] - mass_target) <= mass_tolerance:
                if spin[i] > 0 and J[i] > 0:
                    group_spin.append(spin[i])
                    group_J.append(J[i])

        if group_spin:
            log_spin = np.log10(np.array(group_spin))
            log_J = np.log10(np.array(group_J))
            sort_idx = np.argsort(log_J)
            log_spin_s = log_spin[sort_idx]
            log_J_s = log_J[sort_idx]
            spin_mass[mass_target] = log_spin_s
            J_mass[mass_target] = log_J_s
            ax.plot(log_J_s, log_spin_s, color=cmap(norm(mass_target)), linewidth=1)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Mass (M☉)", rotation=270, labelpad=15)
    ax.set_xlabel("log₁₀(J) [g·cm²/s]")
    ax.set_ylabel("log₁₀(Spin) [rad/s]")
    ax.set_title("log(Spin) vs log(J) for Various Masses")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig("spin_vs_J_all_masses1.png", dpi=300)

# Read all .txt files
txt_files = [f for f in os.listdir('.') if f.endswith('.txt')]

for filename in txt_files:
    with open(filename, 'r') as file:
        lines = file.readlines()
        rho_str = filename.replace("rho_", "").replace(".txt", "")
        try:
            rho_c1 = float(rho_str)
        except ValueError:
            continue

    start_idx = 0
    for i, line in enumerate(lines):
        if line.startswith('ratio'):
            start_idx = i + 1
            break

    for line in lines[start_idx:]:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 9:
            continue
        try:
            mass = float(parts[2])         # in M_sun
            spin = float(parts[5])         # rad/s
            J_over_M2 = float(parts[8])    # J/M^2
            mass_cgs = mass * M_sun
            J = J_over_M2 * (mass_cgs)**2
            if J > 1e69:
                print(f"Large J: {J:.3e} at mass = {mass:.4f}")
            mass_values.append(mass)
            spin_values.append(spin)
            momentum_values.append(J)
            rho_values.append(rho_c1)
        except ValueError:
            continue

# Generate plots and slope analysis
plot_all_masses(rho_values, momentum_values, spin_values, mass_values)
#mass_slope = find_slope(spin_mass, J_mass)
mass_slope = find_slope_savgol(spin_mass,J_mass)
mass_flat, mass_spinup, mass_spindown = inflection_point(mass_slope)

# Results
print("Spin-down (max slope) masses:", mass_spindown)
print("Spin-up (min slope) masses:", mass_spinup)
print("Flat slope masses:", mass_flat)
print(f"Counts → Spin-down: {len(mass_spindown)}, Spin-up: {len(mass_spinup)}, Flat: {len(mass_flat)}")

# Optional: quick plot of classification
def plot_classification(mass_spinup, mass_spindown):
    plt.figure(figsize=(7, 2))
    plt.plot(mass_spinup, [1]*len(mass_spinup), 'g^', label='Spin-Up (Minimum)')
    plt.plot(mass_spindown, [1]*len(mass_spindown), 'rv', label='Spin-Down (Maximum)')
    plt.xlabel("Mass (M☉)")
    plt.yticks([])
    plt.title("Masses with Spin-Up and Spin-Down Turning Points")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("turning_point_classification.png", dpi=300)

plot_classification(mass_spinup, mass_spindown)
