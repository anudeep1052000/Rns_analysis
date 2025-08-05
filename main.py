import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-GUI backend
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize
import math
from matplotlib.colors import Normalize
from scipy.signal import savgol_filter
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

from integrator import integrator1,stepwise_coupled_integrator
from derivative import find_slope_savgol
from compute import drho_dt,epsilon_grav,epsilon_nuc_nu,evolution,guess_timelist
from mpi4py import MPI
from spin_vs_J import plot_all_masses
from rho_vs_J import plot_all_masses1
import os
import subprocess
"""
# === MPI setup ===
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
# === Create output dir (only once) ===
if rank == 0:
    os.makedirs("energy", exist_ok=True)
comm.Barrier()
"""



M_sun = 1.989e33
mass_tolerance = 1e-4

C = 2.99792458e10       # speed of light in cm/s
G = 6.67430e-8          # gravitational constant in cm^3 / g / s^2
seconds_in_year = 365.25 * 24 * 3600 
# Compute KAPPA and KSCALE as in the C code
KAPPA = 1.0e-15 * C**2 / G
KSCALE = KAPPA * G / C**4
mass_values = []       # All mass values from all files
spin_values = []       # All spin values from all files
momentum_values = []   # All J values from all files
rho_values = [] # rho_c for each data point
radius_values = []
pressure_c_values = []


    
def guess_central_temperature(rho_c):

    if rho_c < 1e7:
        # Cold WD regime
        return 1e7  # ~10 million K
    elif 1e7 <= rho_c < 1e9:
        # Warm, possibly accreting WD
        # Linear interpolation from 1e7 to 2e8 K
        return 1e7 + (rho_c - 1e7) * (2e8 - 1e7) / (1e9 - 1e7)
    else:
        # High-density regime near ignition
        # Cap at ~3e8 K
        return 3e8
    


print("Current working directory:", os.getcwd())
#print("Files inside ./results/:", os.listdir('./results'))
txt_files = [f for f in os.listdir('./results/') if f.endswith('.txt')]
#print(txt_files)
for filename in txt_files:
    full_path = os.path.join('./results', filename)
    with open(full_path, 'r') as file:
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
    #print(start_idx)
    count = 0

    for line in lines[start_idx:]:
        
        if line.startswith('DEBUG'):
            if(count ==1):
                Pressure_central = line.split("=")[-1]
                #print("pressure "+str(Pressure_central))
            else:
                count = count+1
        #print(not line.strip())
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 9:
            continue
        try:
            mass = float(parts[2])         # in M_sun
            spin = float(parts[5])         # rad/s
            J_over_M2 = float(parts[8])# J/M^2
            I = float(parts[7])
            radius = float(parts[4])
            mass_cgs = mass * M_sun
            J = (J_over_M2 * (mass_cgs)**2*G)/C
            if J > 1e69:
                print(f"Large J: {J:.3e} at mass = {mass:.4f}")
            mass_values.append(mass)
            spin_values.append(spin)
            momentum_values.append(J)
            rho_values.append(rho_c1)
            radius_values.append(radius*10**5)
            pressure_c_values.append(Pressure_central)
        except ValueError:
            continue





# Call the new plotting function
pressure_c_values = [float(x.strip())/KSCALE for x in pressure_c_values]
#print(pressure_c_values)

print("-------------------------------------")
plot_all_masses(rho_values, momentum_values, spin_values,radius_values,mass_values)
print("plotted successfully")
mass_slope,mass_dj_dt,mass_pc,rho_mass,J_mass = plot_all_masses1(rho_values,momentum_values,radius_values,spin_values,mass_values,pressure_c_values,B=1e10,alpha=math.pi/2)
#print("rho_mass",rho_mass[1.3])
#print("J_mass",J_mass[1.3])
#print(mass_pc)
#print("mass_dhro_dj",mass_slope[1.3])
#print("mass_dj_dt",mass_dj_dt[1.3])

mass_drho_dt = drho_dt(mass_slope,mass_dj_dt,5,6)
print(mass_dj_dt)
#print("mass_dhro_dt",mass_drho_dt[1.3])
#print(rho_mass)
mass_time_list = guess_timelist(mass_drho_dt,rho_mass)
mass_time_list_years = {}
for i in mass_time_list.keys():
    a = np.array(mass_time_list[i])
    a_years = a / seconds_in_year
    mass_time_list_years[i] = a_years
#print("time_list",mass_time_list_years[1.3])


print("-------------------------------")
mass_E_gravity = epsilon_grav(mass_drho_dt,mass_pc,rho_mass)
#mass_dtemp_dt = evolution(mass_E_gravity=mass_E_gravity,mass_rho = rho_mass)
#mass_dtemp_dt = {1.3: [862363.6804808223, 862363.6804808223, 862363.6804808223, 862363.6804808223, 862363.6804808223, 865363.3699941944, 868362.4946031364, 871361.1244857997, 900466.0212702699, 926701.0500753771, 935813.2074136983, 947985.695560131, 961207.5581251568, 973437.4355650673, 986714.6407904312, 1015429.1546918534, 1032927.0384756505, 1049449.4239607672, 1067051.5447818756, 1105530.982094518, 1127478.07973327, 1195869.3270120274, 1219174.2157707433, 1271355.0950564267, 1301341.0684796653, 1330366.7406822105, 1361666.0737976297, 1393093.065852007, 1424636.2410679436, 1457394.688659501, 1495773.7900523737, 1529906.7939009375, 1569712.2062435, 1610799.6788300008, 1787094.8931798362, 1833636.2381294498, 1881516.314393551, 1936500.8705666685, 1987125.0335865538, 2046081.761418171, 2106489.3015035503, 2224655.174136148, 2290626.4016951118, 2356919.7046334143, 2506021.998955617, 2579331.4931737576, 2741585.6825141353, 2830666.6307895356, 2914073.8060532394, 3008991.590128385, 3106900.6919727996, 3207829.2835174347, 3311792.6011559945, 3420090.363630483, 3531500.645629793, 3646064.343624919]}
mass_dtemp_dt = {1.3: [1156120.441065494, 1156120.441065494, 1156120.441065494, 1156120.441065494, 1156120.441065494, 1160135.6903123786, 1164150.921680996, 1168172.0618444202, 1207147.9894203842, 1242273.4748937234, 1254473.07997387, 1270768.5464290695, 1288466.5588089856, 1304840.2716789148, 1322615.8329503264, 1361044.185030056, 1384462.5457662607, 1406565.7623123901, 1430115.837666846, 1481600.344066835, 1510954.4077705943, 1602429.200387344, 1633590.5956927827, 1703354.5876514534, 1743433.2675535588, 1782229.6543743555, 1824051.3285955382, 1866034.4131205904, 1908176.0955355666, 1951934.1866128005, 2003183.5100913853, 2048758.75994134, 2101897.8447056194, 2156742.1747222426, 2391955.3621076285, 2454017.336731757, 2517860.546780623, 2591142.3580331584, 2658612.522249495, 2737164.4485106957, 2817627.076208128, 2974962.96874952, 3062776.8213589555, 3150992.299096158, 3349314.619181404, 3446785.356153355, 3662403.6698225136, 3780728.4729465833, 3891467.479582928, 4017464.029664282, 4147380.542255716, 4281241.324538733, 4419088.061052556, 4562628.5973866945, 4710232.74874877, 4861937.4445397165]}
#print(mass_dtemp_dt)
time_array = np.array(mass_time_list_years[1.3])
dtemp_array = np.array(mass_dtemp_dt[1.3])
for i in mass_drho_dt.keys():
    drho_dt_i = mass_drho_dt[i]
    pc =  mass_pc[i]
    rho = rho_mass[i]
    time_list = mass_time_list[i]
    #print(drho_dt_i)
    #print(pc)
    #print(rho)
    Temp_list,time_list,t_heating,runaway_event  = stepwise_coupled_integrator(time_list=time_list,drho_dt=drho_dt_i,pressure_cenrtral=pc,rho=rho,T_initial=1e8,Xc12_init=0.5,Xo16_init=0.5,t_max= 3e8)
    print(Temp_list)
    print(time_list)
    print(t_heating)
    print(runaway_event)
    
#check_lengths_equal([mass_slope,mass_dj_dt,mass_drho_dt])
print(f"Processed {len(txt_files)} files.")
print(f"Collected {len(mass_values)} data points.")
#print(len(pressure_c_values))
#print(len(rho_values))

