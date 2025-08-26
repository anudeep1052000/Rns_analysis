import matplotlib
matplotlib.use('Agg')  # Non-GUI backend
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from call_self_heat import run_self_heat,run_self_heat_classical
import math
import numpy as np 

K =1.38071e-16
mp = 1.672e-24
G = 6.67430e-8
Gamma =5/3
cv_global = K/(2*(mp*(Gamma-1)))
print(cv_global)


def runaway_event(t, T):
    dTdt_interp = interp1d(time_list, dTdt_list, kind='cubic', fill_value="extrapolate")
    t_dynamic_interp = interp1d(t_list, t_dynamic_list, kind='linear', fill_value="extrapolate")
    dT = dTdt_interp(t)
    t_dyn = t_dynamic_interp(t)
    t_heating = T[0] / abs(dT) if dT != 0 else 1e99
    return t_heating - t_dyn

def integrator(dtemp_list,t_list,t_dynamic):
    time_list = np.array(t_list)
    dTdt_list = np.array(dtemp_list)
    T0 = 1e9
    t_final = 3.15576e10
    

    # Interpolate dT/dt as a function of time
    dTdt_interp = interp1d(time_list, dTdt_list, kind='cubic', fill_value="extrapolate")
    t_dynamic_interp = interp1d(t_list, t_dynamic_list, kind='linear', fill_value="extrapolate")

    # ODE function
    def dTdt_ode(t, T):
        return np.array([dTdt_interp(t)])

    def runaway_event(t, T):
        dT = dTdt_interp(t)
        t_dyn = t_dynamic_interp(t)   # Interpolated t_dynamic at current time
        t_heating = T[0] / abs(dT) if dT != 0 else 1e99
        return t_heating - t_dyn
    runaway_event.terminal = True
    runaway_event.direction = -1
    
    # Integrate using stiff solver
    sol = solve_ivp(
        dTdt_ode,
        [time_list[0],3e4],
        [T0],
        method='BDF',
        events=runaway_event,
        rtol=1e-4,
        atol=1e-6
    )

    # Plot results
    plt.figure(figsize=(10,5))
    plt.plot(sol.t, sol.y[0], label='Temperature [K]')
    plt.xlabel("Time [years]")
    plt.ylabel("T [K]")
    
    plt.title("Temperature Evolution")
    plt.grid(True)
    plt.legend()
    plt.savefig("Temperature_evolution.png")
    #plt.show()

    if sol.t_events[0].size > 0:
        t_ignition = sol.t_events[0][0]
        print(f"🔥 Runaway detected at t = {t_ignition:.2e} years")
        return sol.t, sol.y[0], t_ignition
    else:
        print("✅ No runaway ignition detected")
        return sol.t, sol.y[0], None


def integrator1(dtemp_dt,time_list,t_dynamic):
    dTdt_interp = interp1d(time_list, dtemp_dt, kind='cubic', fill_value="extrapolate")
    t_dynamic_interp = interp1d(time_list, t_dynamic, kind='linear', fill_value="extrapolate")
    T0 = 1e8
    dt = 1e4
    Temp_list = [T0]
    dt_max = 1e10
    dt_min = 1e8
    times = [0]
    t_init = 0
    t_final = 3e12
    t = t_init
    t_dynamic_l = []
    t_heating_l = []
    runaway_detect = False
    while t<t_final:
        
        T_before = Temp_list[-1]
        Dtemp = dTdt_interp(t)
        #print(Dtemp)
        t_dyn = t_dynamic_interp(t)   # Interpolated t_dynamic at current time
        t_dynamic_l.append(t_dyn)
        T_after = T_before+Dtemp*dt
        Temp_list.append(T_after)
        t = t+dt
        times.append(t)
        if((T_after-T_before)/T_before<=0.01):
            dt = dt*2
        else:
            dt = dt/2
        #print(T_after)
        t_heating = T_after / abs(Dtemp) if Dtemp != 0 else 1e99
        t_heating_l.append(t_heating)
        if (t_heating< t_dyn):
            print("runaway occured at time ",times[-2], T_before)
            runaway_detect = True
            break
    """plt.figure(figsize=(10,5))
    plt.plot(times[0:len(times)-1],Temp_list[0:len(Temp_list)-1], label='Temperature [K]')
    plt.xlabel("Time [years]")
    plt.ylabel("T [K]")
    
    plt.title("Temperature Evolution")
    plt.grid(True)
    plt.legend()
    plt.savefig("Temperature_evolution1.png")
    
    plt.figure(figsize=(10, 5))
# Assume time_l, t_dynamic_l, t_heating_l are lists or arrays of the same length
    plt.plot(times[0:len(times)-1], t_dynamic_l, label='Dynamic Timescale [yr]', color='tab:blue')
    plt.plot(times[0:len(times)-1], t_heating_l, label='Heating Timescale [yr]', color='tab:orange')

    plt.xlabel("Time [years]")
    plt.ylabel("Timescale [years]")
    plt.title("Heating vs Dynamic Timescales Over Time")
    plt.yscale('log')  # Optional: log scale is useful for comparing orders of magnitude
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig("T_heating_vs_T_dynamic.png")
"""


# Ensure equal lengths
    #n = min(len(times), len(Temp_list), len(t_heating_l))

    fig, ax1 = plt.subplots(figsize=(10, 5))

     # Plot Temperature on primary y-axis
    ax1.plot(times[0:len(times)-1], Temp_list[0:len(Temp_list)-1], color='crimson', label='Temperature [K]')
    ax1.set_xlabel("Time [seconds]")
    ax1.set_ylabel("Temperature [K]", color='crimson')
    ax1.tick_params(axis='y', labelcolor='crimson')

# Twin y-axis for Heating Timescale
    ax2 = ax1.twinx()
    ax2.plot(times[0:len(times)-1], t_heating_l, color='tab:orange', label='Heating Timescale [yr]')
    ax2.set_ylabel("Heating Timescale [seconds]", color='tab:orange')
    ax2.set_yscale('log')
    ax2.tick_params(axis='y', labelcolor='tab:orange')

# Grid and title
    fig.suptitle("Temperature and Heating Timescale Over Time")
    ax1.grid(True, which='both', linestyle='--', alpha=0.6)

# Legends
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper left')

    plt.tight_layout()
    plt.savefig("Temp_and_T_heating_combined.png")
    plt.close()
    return Temp_list,times,runaway_detect
            


def stepwise_coupled_integrator(time_list,drho_dt,pressure_cenrtral,rho,mass,T_initial=1e8 , Xc12_init=0.5,Xo16_init =0.5,t_max=3e11):
    # Initial conditions
    print("Started this mass",mass)
    T = T_initial  # K
    t = 0.0
    dt = 1e4
    dt_min = 1e2
    dt_max = 1e8
    #t_dynamic_interp = interp1d(rho, t_dynamic_list, kind='linear', fill_value="extrapolate")
    pressure_cenrtral_interp = interp1d(time_list,pressure_cenrtral, kind='linear', fill_value="extrapolate")
    rho_interp = interp1d(time_list,rho, kind='linear', fill_value="extrapolate")
    drho_dt_interp = interp1d(time_list,drho_dt, kind='linear', fill_value="extrapolate")
    X_p_init =0
    rho_init = rho_interp(0)
    pc_init =  pressure_cenrtral_interp(0)
    X_he4_init =0
    X_ne20_init =0
    X_na23_init =0
    X_mg24_init =0
    t_dyn_init = 1/math.sqrt(G*rho_init)
    epsilon_gravity_init = pc_init*drho_dt_interp(0)/rho_init**2
    T_initial,cv_temp,p_temp = run_self_heat_classical(rho=rho_init,dt=5)
    print(T_initial)
    # Storage
    T_list = [float(T_initial)]
    t_list = [t]
    rho_list = [rho[0]]
    drho_dt = [drho_dt[0]]
    runaway_detected = False
    t_heat_l = [0]
    T = float(T_initial)
    with open("stellar_evolution_output"+str(mass)+".dat", "w") as f:
        f.write("# time(s) T(K) rho(g/cm^3) p_central(dyne/cm^2) epsilon_nuc epsilon_nu epsilon_gravity cv  t_heat t_dyn X_c12 X_o16 X_p X_he4 X_ne20 X_na23 X_mg24\n")
        f.flush()
        f.write(f"{t:.6e} {float(T_initial):.6e} {rho_init:.6e} {pc_init:.6e} "
                f"{0:.6e} {0:.6e} {epsilon_gravity_init:.6e} "
                f"{float(cv_temp):.6e} {0:.6e} {t_dyn_init:.6e} "
                f"{Xc12_init:.6e} {Xo16_init:.6e} "
                f"{X_p_init:.6e} {X_he4_init:.6e} "
                f"{X_ne20_init:.6e} {X_na23_init:.6e} {X_mg24_init:.6e}\n")
        f.flush()
        print(f"T={float(T_initial):.3e} K, t={t:.3e} s, rho={rho_init:.3e}, p_central={pc_init:.3e}, "
              f"E_nuc={0:.3e}, E_nu={0:.3e}, E_grav={epsilon_gravity_init:.3e}, "
              f"t_heat={0:.3e}, t_dyn={t_dyn_init:.3e}, X_C12={Xc12_init:.3e}, X_O16={Xo16_init:.3e}")
        while t <= t_max:
            rho_temp = rho_interp(t)
        # Call self_heat function (user-defined) to get epsilon_nuc, epsilon_nu, updated X_i
            epsilon_nuc, epsilon_nu,cv,x_p,x_he4,x_c12,x_o16,x_ne20,x_na23,x_mg24,dt_error = run_self_heat(dt,rho_temp,T,X_p_init,X_he4_init, Xc12_init,Xo16_init,X_ne20_init,X_na23_init,X_mg24_init)
            Xc12_init = x_c12
            Xo16_init = x_o16
            X_p_init = x_p
            X_he4_init = x_he4
            X_ne20_init = x_ne20
            X_na23_init = x_na23
            X_mg24_init = x_mg24
        
            drho_dt_temp = drho_dt_interp(t+dt+dt_error)
            p_central_temp = pressure_cenrtral_interp(t+dt+dt_error)
            if(dt_error<0):
                print("Nuclear network didnt complete full run updating the time_step ")
          
        

             # Gravitational heating 
            epsilon_gravity =(drho_dt_temp*p_central_temp)/rho_temp**2 
        
    

             # Compute dT/dt
            dTdt = (epsilon_nuc - epsilon_nu + epsilon_gravity) / cv
            
            # Advance time and temperature 
            T_next = T + dTdt * (dt+dt_error)
            t += dt+dt_error

            # Store
            T_list.append(T_next)
            t_list.append(t)
            # Compute dynamic timescale
            #t_dyn = t_dynamic_interp(rho_temp)
            t_dyn = 1/math.sqrt(G*rho_temp)

            # Compute heating timescale
            t_heat = abs(T / dTdt) if dTdt != 0 else 1e99
            t_heat_l.append(t_heat)
            # Runaway condition
            if t_heat < t_dyn:
                print(f"🔥 Runaway at t = {t:.2e} s, T = {T:.2e} K")
                runaway_detected = True
                break
            if(Xc12_init ==0):
                print("carbon Burning is Completed")

        

            # Adaptive time stepping
            delta_T = abs(dTdt * (dt+dt_error))
            if delta_T/T < 0.01:
                #dt = min(dt * 1e2, dt_max)
                dt =dt*2
            elif delta_T/T > 1:
                #dt = max(dt / 1e2, dt_min)
                dt =dt/2
        

        # Update T for next step
            T = T_next
        
            f.write(f"{t:.6e} {T:.6e} {rho_temp:.8e} {p_central_temp:.6e} "
                f"{epsilon_nuc:.6e} {epsilon_nu:.6e} {epsilon_gravity:.6e} "
                f"{cv: .6e} {t_heat:.6e} {t_dyn:.6e} "
                f"{Xc12_init:.6e} {Xo16_init:.6e} "
                f"{X_p_init:.6e} {X_he4_init:.6e} "
                f"{X_ne20_init:.6e} {X_na23_init:.6e} {X_mg24_init:.6e}\n")
            f.flush()
            print(f"T={T_next:.3e} K, t={t:.3e} s, rho={rho_temp:.3e}, p_central={p_central_temp:.3e}, "
              f"E_nuc={epsilon_nuc:.3e}, E_nu={epsilon_nu:.3e}, E_grav={epsilon_gravity:.3e}, "
              f"t_heat={t_heat:.3e}, t_dyn={t_dyn:.3e}, X_C12={Xc12_init:.3e}, X_O16={Xo16_init:.3e}")
    
   
    #print()
    fig, ax1 = plt.subplots(figsize=(10, 5))

     # Plot Temperature on primary y-axis\
         
    ax1.plot(t_list, T_list, color='crimson', label='Temperature [K]')
    ax1.set_xlabel("Time [seconds]")
    ax1.set_ylabel("Temperature [K]", color='crimson')
    ax1.tick_params(axis='y', labelcolor='crimson')

# Twin y-axis for Heating Timescale
    ax2 = ax1.twinx()
    ax2.plot(t_list, t_heat_l, color='tab:orange', label='Heating Timescale [yr]')
    ax2.set_ylabel("Heating Timescale [seconds]", color='tab:orange')
    ax2.set_yscale('log')
    ax2.tick_params(axis='y', labelcolor='tab:orange')

# Grid and title
    fig.suptitle("Temperature and Heating Timescale Over Time")
    ax1.grid(True, which='both', linestyle='--', alpha=0.6)

# Legends
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper left')

    plt.tight_layout()
    plt.savefig("Temp_and_T_heating_combined_updated"+str(mass)+".png")
    plt.close()

    return t_list, T_list,t_heat_l,runaway_detected

def classical_integration_case(time_list,drho_dt,pressure_cenrtral,rho,mass,T_initial=1e8 , Xc12_init=0.5,Xo16_init =0.5,t_max=3e11):
    
    print("Started this mass",mass)
    T = T_initial  # K
    t = 0.0
    dt = 1e4
    dt_min = 1e2
    dt_max = 1e8
    #t_dynamic_interp = interp1d(rho, t_dynamic_list, kind='linear', fill_value="extrapolate")
    pressure_cenrtral_interp = interp1d(time_list,pressure_cenrtral, kind='linear', fill_value="extrapolate")
    rho_interp = interp1d(time_list,rho, kind='linear', fill_value="extrapolate")
    drho_dt_interp = interp1d(time_list,drho_dt, kind='linear', fill_value="extrapolate")
    X_p_init =0
    rho_init = rho_interp(0)
    pc_init =  pressure_cenrtral_interp(t)
    X_he4_init =0
    X_ne20_init =0
    X_na23_init =0
    X_mg24_init =0
    t_dyn_init = 1/math.sqrt(G*rho_init)
    epsilon_gravity_init = pc_init*drho_dt_interp(0)/rho_init**2
    # Storage
    T_list = []
    t_list = []
    rho_list =[]
    count=0
    runaway_detected = False
    T0 = float(T_initial)       # initial temperature (K)
    rho0 = rho[0]  # reference density
    gamma = 5/3    # adiabatic index for ideal, non-relativistic gas
    T_adiabatic = []

    t_heat_l = []
    
    while t<=t_max:
        rho_temp = rho_interp(t) #interpolating the Density with time
        drho_dt_temp = drho_dt_interp(t) #interpolating Drho/dt with time_list
        pc_temp = pressure_cenrtral_interp(t) #interpolating pressure with time_list
        #print(t)
       #print(rho_temp)
        
        
        if count ==0:
            #finding initial temperature where the Epsilon_nuc and Epsilon_nuetrino are zero
            T_initial,cv,pcentral_temp = run_self_heat_classical(dt,rho=rho_temp,T=T,xp=0,xhe4=0,xc12=0.5,xo16=0.5,xne20=0,xna23=0,xmg24=0)
            T=float(T_initial)
        #print("central pressure from nuclear code ",pcentral_temp)
        #print("central pressure from rns ",pc_temp)
        #print("cv",cv)
        #print(type(cv))
        #print("rho_temp",rho_temp)
        #print("drho_dt",drho_dt_temp)
        epsilon_gravity =float(drho_dt_temp*float(pc_temp))/rho_temp**2 #calculating the Epsilon_gravity term
        #print("epsilon_gravity",epsilon_gravity)
        
        dT_dt = epsilon_gravity/float(cv_global) #cv_global = cv for ideal gas law case
        #print("dT_dt",dT_dt)
        T_next = T+dT_dt*dt #updating the temperature form dT_dt and dt
        #print(type(T_next))
        t_adiabatic = float(T_initial) * (rho_temp/ rho0)**(gamma - 1) #calculating the adiabatic temperature
        
        T_list.append(T)
        T_adiabatic.append(t_adiabatic)
        t_list.append(t)
        #print("T_adiabatic",t_adiabatic)
        #print("Temperature_classical",T)
        #print(type(t_adiabatic))
   
        rho_list.append(float(rho_temp))
        # Adaptive time stepping
        delta_T = abs(dT_dt * (dt))
        if delta_T/T < 0.01:
            #dt = min(dt * 1e2, dt_max) 
            dt =2*dt #increasing the timestep
        elif delta_T/T > 1:
            #dt = max(dt / 1e2, dt_min)
            dt =dt/2 # decreasing the time_step
        
        #updating the temperature time for next iteration
        T=T_next
        t+=dt
        #print("Time",t)
        count+=1
    
    plt.plot(rho_list,T_list, label='Classical_ideal_law', color='tab:orange')
    plt.plot(rho_list, T_adiabatic, label='Adiabatic (γ=5/3)', color='tab:blue', linestyle='--')


    plt.xlabel("rho [gm/cm^3]")
    plt.ylabel("Temperature [K]")
    plt.title("Temperature vs Rho in classical case no nuclear burning")
    #plt.yscale('log')  # Optional: log scale is useful for comparing orders of magnitude
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig("Temperature_vs_Rho_classical_ideal_law"+str(mass)+".png")
    return T_list,t_list,rho_list
