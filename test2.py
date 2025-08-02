import os
import numpy as np
import csv
import matplotlib
matplotlib.use('Agg')  # Non-GUI backend
import matplotlib.pyplot as plt
import math

M_sun = 1.989e33
mass_tolerance = 1e-5
G = 6.67430e-8
C = 2.99792458e10
mass_values = []     # store all mass values from all files
spin_values = []     # store all spin values from all files
momentum_values = [] # store all J values from all files
rho_c = []
radius_values = []

def compute_breakup_omega(mass,Radius):
    M = mass*M_sun
    omega = math.sqrt(G * M / (Radius*10**5)**3)
    return omega

def plot(rho_c,J,spin,mass,radius_values,mass_value):
    group_rho =[]
    group_spin = []
    group_J = []
    group_radius = []
    for  i in range(0,len(mass)):
        if (abs(mass[i]-mass_value)<=mass_tolerance):
            if spin[i]>0 and J[i]>0:
                group_rho.append(rho_c[i])
                group_spin.append(spin[i])
                group_J.append(J[i])
                group_radius.append(radius_values[i])
    log_spin = np.log10(group_spin)
    log_rho = np.log10(group_rho)
    log_radius = np.log10(group_radius)

    log_J = np.log10(group_J)
    for i in range(0,len(log_J)):
        if log_J[i]>=69:
            print("central_density"+str(group_rho[i]))
    #print(log_J)
    #print(log_spin)[sort_indices]
    
    sort_indices = np.argsort(log_J)
    print(np.array(group_rho)[sort_indices])
    log_J_s = log_J[sort_indices]
    log_spin_s = log_spin[sort_indices]
    log_rho_s = log_rho[sort_indices]
    log_radius_s = log_radius[sort_indices]
    Omega  = compute_breakup_omega(mass_value,10**log_radius_s[-1])
    print(Omega)
    print(10**log_spin_s[-1])
    # Plot spin vs log J with dual y-axes (log_spin on left, log_rho on right)
    fig, ax1 = plt.subplots()

# Left y-axis: log(spin)
    ax1.plot(log_J_s, log_spin_s, 'o-', color='blue', label='log Spin')
    ax1.set_xlabel("Log J (g cm²/s)")
    ax1.set_ylabel("log Spin (rad/s)", color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')

# Right y-axis: log(rho)
    ax2 = ax1.twinx()
    ax2.plot(log_J_s, log_rho_s, 'o--', color='orange', label='log ρ_c')
    ax2.set_ylabel("log ρ_c (g/cm³)", color='orange')
    ax2.tick_params(axis='y', labelcolor='orange')

# Final formatting
    fig.tight_layout()
    plt.savefig("spin_and_rho_vs_logJ.png", dpi=300)
#plt.show()

#log_J_s  = log_J.sort()
#log_spin_s = log_spin.sort()
    """   # Plot spin vs rho_c
    plt.figure()
    plt.plot(log_J_s, log_spin_s, 'o-', label=f"Spin vs ρ_c for M ~ {mass_value:.5f}")
    #plt.plot(log_J_s, log_rho_s, 'o-', label=f"Spin vs ρ_c for M ~ {mass_value:.5f}")
    plt.xlabel("Log J (g cm^2)")
    plt.ylabel("log_Spin (rad/s)")
    #Xplt.grid(True)
    #plt.legend()
    plt.tight_layout()
    plt.savefig("spin_vs_rhoc.png", dpi=300)
    """
   

# List all .txt files in current directory
txt_files = [f for f in os.listdir('.') if f.endswith('.txt')]

for filename in txt_files:
    #print(f"Processing file: {filename}")
    with open(filename, 'r') as file:
        lines = file.readlines()
        rho_str = filename.replace("rho_", "").replace(".txt", "")
        rho_c1 = float(rho_str)
        #print(rho_c1) 
        
    # Find the line where the table starts (looking for 'ratio' in the header)
    start_idx = 0
    for i, line in enumerate(lines):
        if line.startswith('ratio'):
            start_idx = i + 1  # Skip the header line
            break

    for line in lines[start_idx:]:
        if not line.strip():
            continue
        parts = line.split()
        # We need at least 9 columns to safely extract parts[8]
        if len(parts) < 9:
            continue

        try:
            mass = float(parts[2])         # 3rd column, solar masses
            spin = float(parts[5])         # 6th column, rad/s
            J_over_M2 = float(parts[8])     # 9th column, J/M^2
            
            radius = float(parts[4])
            #print(J_over_M2)
            mass_cgs = mass * M_sun
            J = (G*J_over_M2 * (mass_cgs)**2)/C
            #print(J)
            if J > 1e69:
                print(f"Large angular momentum detected: J = {J:.3e} at mass = {mass}")
            mass_values.append(mass)
            spin_values.append(spin)
            momentum_values.append(J)
            rho_c.append(rho_c1)
            radius_values.append(radius)
            #with open('results.csv', 'a', newline='') as csvfile:
                #writer = csv.writer(csvfile)

                # Write header only if file is empty
                #if csvfile.tell() == 0:
                    #writer.writerow(['rho_c', 'mass (g)', 'spin (rad/s)', 'angular_momentum (g·cm²/s)'])
                #writer.writerow([rho_c, mass_values, spin_values, momentum_values])
        except ValueError:
            #print(f"Skipping line due to ValueError: {line.strip()}")
            continue

plot(rho_c,momentum_values,spin_values,mass_values,radius_values,1.42)
print(f"Processed {len(txt_files)} files.")
print(f"Collected {len(mass_values)} data points.")
