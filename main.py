import os
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-GUI backend for headless systems
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

# Local module imports
from integrator import integrator1, stepwise_coupled_integrator,classical_integration_case
from derivative import find_slope_savgol
from compute import drho_dt, epsilon_grav, epsilon_nuc_nu, evolution, guess_timelist
from mpi4py import MPI
from spin_vs_J import plot_all_masses
from rho_vs_J import plot_all_masses1
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-evolution_type", help="Type of evolution to run")
args = parser.parse_args()
evolution_type = args.evolution_type
print(args.evolution_type)
# ====================
# Physical constants
# ====================
M_sun = 1.989e33             # Solar mass in g
C = 2.99792458e10           # Speed of light in cm/s
G = 6.67430e-8              # Gravitational constant in CGS
seconds_in_year = 365.25 * 24 * 3600
mass_tolerance = 1e-4

# Scaling constants (consistent with C code)
KAPPA = 1.0e-15 * C**2 / G
KSCALE = KAPPA * G / C**4

# ====================
# Data containers
# ====================
mass_values = []       # Stellar mass in M_sun
spin_values = []       # Spin (rad/s)
momentum_values = []   # Angular momentum J (CGS)
rho_values = []        # Central density for each data point (g/cm^3)
radius_values = []     # Equatorial radius (cm)
pressure_c_values = [] # Central pressure (scaled)

# ====================
# Helper function: initial guess for central temperature
# ====================
def guess_central_temperature(rho_c):
    """Return an estimated central temperature for a WD given central density."""
    if rho_c < 1e7:
        return 1e7  # Cold WD regime (~10 MK)
    elif 1e7 <= rho_c < 1e9:
        # Interpolate between 10 MK and 200 MK
        return 1e7 + (rho_c - 1e7) * (2e8 - 1e7) / (1e9 - 1e7)
    else:
        # High-density regime, cap at ~300 MK
        return 3e8

# ====================
# Load data from ./results/ directory
# ====================
print("Current working directory:", os.getcwd())
txt_files = [f for f in os.listdir('./results/') if f.endswith('.txt')]

for filename in txt_files:
    full_path = os.path.join('./results', filename)
    with open(full_path, 'r') as file:
        lines = file.readlines()
        rho_str = filename.replace("rho_", "").replace(".txt", "")
        try:
            rho_c1 = float(rho_str)
        except ValueError:
            continue

    # Locate the start of data (after 'ratio' header)
    start_idx = 0
    for i, line in enumerate(lines):
        if line.startswith('ratio'):
            start_idx = i + 1
            break

    count = 0
    Pressure_central = None

    # Parse each data line
    for line in lines[start_idx:]:
        if line.startswith('DEBUG'):
            # Extract central pressure from debug lines
            if count == 1:
                Pressure_central = line.split("=")[-1]
            else:
                count += 1
            continue

        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 9:
            continue

        try:
            mass = float(parts[2])          # in M_sun
            spin = float(parts[5])          # rad/s
            J_over_M2 = float(parts[8])     # J/M^2
            I = float(parts[7])             # Moment of inertia (not used here)
            radius = float(parts[4])        # in km (convert later)

            # Convert mass to CGS and compute J
            mass_cgs = mass * M_sun
            J = (J_over_M2 * (mass_cgs)**2 * G) / C

            if J > 1e69:
                print(f"Warning: Large J {J:.3e} at mass {mass:.4f}")

            # Store results
            mass_values.append(mass)
            spin_values.append(spin)
            momentum_values.append(J)
            rho_values.append(rho_c1)
            radius_values.append(radius * 1e5)  # km -> cm
            pressure_c_values.append(Pressure_central)
        except ValueError:
            continue

# Scale central pressure
pressure_c_values = [float(x.strip()) / KSCALE for x in pressure_c_values]

# ====================
# Plot spin vs J and related diagnostics
# ====================
print("-------------------------------------")
plot_all_masses(rho_values, momentum_values, spin_values, radius_values, mass_values)
print("Spin vs J plot generated successfully.")

mass_slope, mass_dj_dt, mass_pc, rho_mass, J_mass = plot_all_masses1(
    rho_values, momentum_values, radius_values, spin_values, mass_values, pressure_c_values, B=1e10, alpha=math.pi/2
)

# Compute d(rho)/dt for each mass
mass_drho_dt = drho_dt(mass_slope, mass_dj_dt, 5, 6)
print("dJ/dt computed.")

# Generate time lists for evolution
mass_time_list = guess_timelist(mass_drho_dt, rho_mass)
mass_time_list_years = {m: np.array(t) / seconds_in_year for m, t in mass_time_list.items()}

# Compute gravitational energy release\mass_E_gravity = epsilon_grav(mass_drho_dt, mass_pc, rho_mass)


time_array = np.array(mass_time_list_years[1.3])
#evolution_type = input("Enter which evolution to consider")

# ====================
# Stepwise thermal + structural evolution
# ====================
for m in mass_drho_dt.keys():
    drho_dt_i = mass_drho_dt[m]
    pc = mass_pc[m]
    rho = rho_mass[m]
    time_list = mass_time_list[m]

    print(f"Evolving WD mass {m}")
    if(evolution_type=="classical"):
        #classical ideal gas law case and without E_nuclear and E_nuetrino 
        Temp_list,time_list,rho_list  = classical_integration_case(time_list=time_list,drho_dt=drho_dt_i,pressure_cenrtral=pc,rho=rho,T_initial=1e8,Xc12_init=0.5,Xo16_init=0.5,t_max= 8e11,mass=i)
        print("Temperature",Temp_list)
        print("Time",time_list)
        print("Rho",rho_list)
    else:
        Temp_list, time_list, t_heating, runaway_event = stepwise_coupled_integrator(
        time_list=time_list,
        drho_dt=drho_dt_i,
        pressure_cenrtral=pc,
        rho=rho,
        T_initial=1e8,
        Xc12_init=0.5,
        Xo16_init=0.5,
        t_max=3e12,
        mass=m
    )
        print("Temperature profile:", Temp_list)
        print("Time profile:", time_list)
        print("Heating timescale:", t_heating)
        print("Runaway event:", runaway_event)

    

    


    
print(f"Processed {len(txt_files)} files.")
print(f"Collected {len(mass_values)} data points.")
