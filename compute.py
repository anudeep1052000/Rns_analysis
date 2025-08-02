import math
from call_self_heat import run_self_heat
M_sun = 1.989e33

mass_tolerance = 1e-4

C = 2.99792458e10       # speed of light in cm/s
G = 6.67430e-8          # gravitational constant in cm^3 / g / s^2
seconds_in_year = 365.25 * 24 * 3600 


def compute_breakup_omega(mass,Radius):
    M = mass*M_sun
    omega = math.sqrt(G * M / (Radius*10**5)**3)
    return omega

def drho_dt(mass_slope,mass_dj_dt,e_nuc,e_nuetrino):
    mass_drho_dt = {}
    for i in mass_slope.keys():
        slope_drho_dt = []
        a = len(mass_slope)
        b = len(mass_dj_dt)
        #print(a)
        #print(b)
        for j in  range(0,len(mass_slope[i])):
            drho = mass_slope[i][j] * mass_dj_dt[i][j]
            slope_drho_dt.append(drho)
        mass_drho_dt[i] = slope_drho_dt
    #print(mass_drho_dt)
    return mass_drho_dt
    #dT_c_dt = (1 / C_V) * [ e_nuc- e_nuetrino + (P_c / rho_c^2) * (drho_c/dt) ]

def guess_timelist(mass_dhro_dt,rho_mass):
      mass_time_list = {}
      for i in mass_dhro_dt.keys():
          time_list = [0.0]
          for j in range(1,len(mass_dhro_dt[i])):
              if j==0:
                  time_list[j] = 0
              else:
                  if abs(mass_dhro_dt[i][j])<1e-30:
                      dt =1e3
                  else:
                      dt = (rho_mass[i][j]-rho_mass[i][j-1])/mass_dhro_dt[i][j]
                  if(dt<0):
                      dt = abs(dt)
                      
                  elif(dt==0):
                      dt =1e3
                  time_list_temp = time_list[j-1]+dt
                  if math.isclose(time_list_temp, time_list[-1], rel_tol=1e-12, abs_tol=1e-6):
                      time_list_temp+= 1e30
                      
                  time_list.append(time_list_temp)
          mass_time_list[i] = time_list
      return mass_time_list



def epsilon_grav(mass_drho_dt,mass_pc,rho_mass):
    mass_E_gravity ={}
    for i in mass_drho_dt.keys():
        E_grvaity = (np.array(mass_drho_dt[i])*np.array(mass_pc[i]))/np.array(rho_mass[i])**2
        mass_E_gravity[i] = E_grvaity
    return mass_E_gravity



def epsilon_nuc_nu(rho,mass):
    #T = guess_central_temperature(rho)
    eps_nuc, eps_nu, cv = run_self_heat(rho,1e9)
    print(eps_nuc,eps_nu,cv)
    print("Mass = "+str(mass)+"rho ="+str(rho)+"successfully done")
    return eps_nuc,eps_nu,cv


def evolution(mass_E_gravity,mass_rho):
    mass_Temp_dt = {}
    for i in mass_E_gravity.keys():
        dTemp_dt = []
        rho_temp = mass_rho[i]
        print(rho_temp)
        for j  in range(0,len(rho_temp)):
            eps_nuc,eps_nu,cv = epsilon_nuc_nu(rho_temp[j],i)
            eps_gravity_temp  = mass_E_gravity[i][j]
            dTemp_dt1 = (eps_nuc-eps_nu+eps_gravity_temp)/cv
            dTemp_dt.append(dTemp_dt1)
        mass_Temp_dt[i] = dTemp_dt
    return mass_Temp_dt

def check_lengths_equal(dicts):
    """
    Check if all the dictionaries have lists of the same length for each mass key.

    Parameters:
        dicts : list of dictionaries to compare (e.g., [mass_slope, mass_dj_dt, mass_drho_dt])

    Returns:
        True if all keys exist in all dicts and have same length of lists, else False with details.
    """
    keys = dicts[0].keys()
    for key in keys:
        lengths = [len(d[key]) for d in dicts if key in d]
        if len(set(lengths)) != 1:
            print(f"❌ Mismatch at mass={key}: lengths = {lengths}")
            return False
    print("✅ All dictionaries have lists of equal length for each mass key.")
    return True


def evolution1(mass_E_gravity, mass_rho):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # 1) Prepare tasks as list of (mass, rho, idx) for indexing
    tasks = []
    for mass in mass_rho:
        print(mass)
        rhos = mass_rho[mass]
        for idx, rho in enumerate(rhos):
            tasks.append((mass, rho, idx))

    # 2) Scatter tasks evenly to each rank
    my_tasks = tasks[rank::size]

    # 3) Each rank computes results for their tasks
    my_results = []
    for mass, rho, idx in my_tasks:
        eps_nuc, eps_nu, cv = epsilon_nuc_nu(rho, mass)
        eps_gravity_temp = mass_E_gravity[mass][idx]
        dTemp_dt = (eps_nuc - eps_nu + eps_gravity_temp) / cv
        my_results.append((mass, idx, dTemp_dt))

    # 4) Gather all results back to rank 0
    all_results = comm.gather(my_results, root=0)

    # 5) Rank 0 reassembles the dictionary
    mass_Temp_dt = None
    if rank == 0:
        mass_Temp_dt = {}
        # Flatten list of lists
        flat_results = [item for sublist in all_results for item in sublist]
        # Group by mass and sort by index
        for mass, idx, dTemp_dt in flat_results:
            if mass not in mass_Temp_dt:
                mass_Temp_dt[mass] = []
        # Initialize empty lists with the right size
        for mass in mass_Temp_dt:
            mass_Temp_dt[mass] = [0] * len(mass_rho[mass])
        # Fill in the values at correct indices
        for mass, idx, dTemp_dt in flat_results:
            mass_Temp_dt[mass][idx] = dTemp_dt

    return mass_Temp_dt
