import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-GUI backend
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize

M_sun = 1.989e33
mass_tolerance = 1e-5
mass_values = []       # All mass values from all files
spin_values = []       # All spin values from all files
momentum_values = []   # All J values from all files
rho_values = []        # rho_c for each data point
spin_mass = {}
J_mass = {}

def find_turning_point(slopes):
    for i in range(1, len(slopes)):
        if slopes[i - 1] > 0 and slopes[i] < 0:
            return 'maximum'  # Positive to negative slope → max
        elif slopes[i - 1] < 0 and slopes[i] > 0:
            return 'minimum'  # Negative to positive slope → min

    return None  # No turning point found

def find_turning_point1(slope):
    for i in range(1, len(slope)):
        if slope[i - 1] > 0 and slope[i] < 0:
            print(slope[i-1],slope[i])
            return 'maximum'  # Positive to negative slope → max
        elif slope[i - 1] < 0 and slope[i] > 0:
            print(slope[i-1],slope[i])
            return 'minimum'  # Negative to positive slope → min

    return None

def find_slope(spin_mass,J_mass):
    mass_slope = {}
    for i in spin_mass.keys():
        slope = []
        a = spin_mass[i]
        b = J_mass[i]
        for j in range(1,len(a)):
            dy = a[j]-a[j-1]
            dx = b[j]-b[j-1]
            slope.append(dy/dx)
        mass_slope[i] = slope
    #print(mass_slope)
    return mass_slope
# Plotting function
def inflection_point(mass_slope):
    mass1 = []
    mass_spindown = []
    mass_spinup = []
    for i in mass_slope.keys():
        #print(i)
        c = mass_slope[i]
        #print(c)
        for j, s in enumerate(c):
            #print(j,s)
            if abs(s)==0:
                print("yes")
                mass1.append(i)

        if(find_turning_point(c)=="maximum"):
            mass_spindown.append(i)
        elif(find_turning_point(c)=="minimum"):
            mass_spinup.append(i)

    return mass1,mass_spinup,mass_spindown
def plot_all_masses(rho, J, spin, mass):
    mass_targets = np.linspace(1.2, 1.49, 50)
    #cmap = get_cmap('viridis')
    cmap = matplotlib.colormaps['viridis']
    norm = Normalize(vmin=mass_targets.min(), vmax=mass_targets.max())
    fig, ax = plt.subplots(figsize=(8, 6))
    plt.figure(figsize=(8,6))

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
            log_J_s = log_J[sort_idx]
            log_spin_s = log_spin[sort_idx]
            spin_mass[mass_target] = log_spin_s #added sorted spin and J into the mass dictionary
            J_mass[mass_target] = log_J_s

            plt.plot(log_J_s, log_spin_s, color=cmap(norm(mass_target)), linewidth=1)

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
    plt.savefig("spin_vs_J_all_masses1.png", dpi=300)

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

# Call the new plotting function
plot_all_masses(rho_values, momentum_values, spin_values, mass_values)
#print(spin_mass)
#print(J_mass)
a = find_slope(spin_mass,J_mass)
mass1,mass_spinup,mass_spindown = inflection_point(a)
#print(mass_spindown)
#print(mass_spinup)
print(mass1)
print(len(mass_spindown))
print(len(mass_spinup))
#print(f"Processed {len(txt_files)} files.")
#print(f"Collected {len(mass_values)} data points.")
print(len(mass1))
def plot_classification(mass_spinup, mass_spindown):
    plt.figure(figsize=(7, 2))
    plt.plot(mass_spinup, [1]*len(mass_spinup), 'g^', label='Spin-Up (Minimum)')
    plt.plot(mass_spindown, [1]*len(mass_spindown), 'rv', label='Spin-Down (Maximum)')
    plt.xlabel("Mass (M☉)")
    plt.yticks([])
    plt.title("Masses with Spin-Up and Spin-Down Turning Points")
    plt.legend()
    plt.grid(True)
    #plt.xlim(1.2, 1.5)
    plt.tight_layout()
    plt.savefig("turning_point_classification2.png", dpi=300)

plot_classification(mass_spinup, mass_spindown)
