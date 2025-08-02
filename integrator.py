import matplotlib
matplotlib.use('Agg')  # Non-GUI backend
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

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
        if(T_after-T_before<=1e5):
            dt = min(dt*1e8,dt_max)
        elif(T_after-T_before>=1e6):
            dt = max(dt/1e8,dt_min)
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
            
