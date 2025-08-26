# Stellar Temperature Evolution of White Dwarfs

This project simulates the **temperature evolution of white dwarfs of constant mass**, including:  

1. **Classical Evolution** – heating occurs **only due to compressional work** (no nuclear reactions).  
2. **Full Nuclear Evolution** – heating includes **both compressional work and nuclear reactions**.  

---

## **Project Structure**

|main.py # Entry point: choose evolution type and initial parameters
├── integrator.py # Contains temperature evolution methods and integrators
├── compute.py # Contains utility functions (e.g., guess_time_list)
├── spin_vs_j.py # Analyzes RNS data and plots spin-down vs angular momentum
├── call_self_heat.py # Wrapper for self_heating network calls
└── self_heating_network/ # Nuclear network data and output files


---

### **File Descriptions**

- **main.py**  
  - Entry point for the simulation.  
  - Accepts `evolution_type` as a parameter (classical or nuclear).  
  - Initializes input parameters such as initial abundances, temperature, and maximum evolution time.  

- **integrator.py**  
  - Contains methods to integrate temperature evolution over time:  
    1. **`stepwise_coupled_integrator`**  
       - Used for **nuclear evolution**.  
       - Inputs: `time_list`, `density_list`, `drho_dt_list`, `pressure_list`, initial abundances, and total evolution time.  
       - Computes temperature using compressional heating + nuclear reactions.  
       - Implements **adaptive timestep** and checks for **runaway events** (heating timescale < dynamic timescale).  
    2. **`classical_integration_case`**  
       - Used for **classical evolution (ideal gas law)**.  
       - Inputs are similar to the nuclear integrator.  
       - Computes temperature assuming no nuclear burning; heat capacity `cv` is calculated from the ideal gas law.  
       - Compares numerical evolution to **adiabatic evolution** (\(T \propto \rho^{\gamma-1}\), with \(\gamma=5/3\)).  

- **compute.py**  
  - Contains utility functions such as `guess_time_list`.  
  - Generates a **time array** based on `density_list` and `drho_dt_list`.  
  - Calculates **compression timescales** for each density element and accumulates to get the total evolution time list.  

- **spin_vs_j.py**  
  - Processes data from the **RNS code**.  
  - Extracts physical quantities such as density, central pressure, angular momentum, and radius for a given mass.  
  - Plots **log(Ω) vs log(J)** to visualize spin-down effects across different masses.  

- **call_self_heat.py**  
  - Provides wrapper functions for the **nuclear network**:  
    1. **`run_self_heat`** → evolves temperature, returns time-averaged nuclear energy and neutrino losses, updated abundances, and heat capacity.  
    2. **`run_self_heat_classical`** → returns **initial temperature** (where nuclear energy and neutrino losses are zero) and heat capacity, without evolving the temperature.  

- **self_heating_network/**  
  - Directory containing **nuclear network data files** and output from `self_heat.py`.  
  - Includes `helm_table.dat` and `abundances_vs_time.dat`.  

---

## **Running the Simulation**

### Classical Evolution (Ideal Gas Law)
python main.py -evolution_type classical
### Nuclear Evolution 
pyhton main.py
