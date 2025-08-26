import subprocess
import numpy as np

def run_self_heat(
    dt, rho, T=1e9, xp=0, xhe4=0, xc12=0.5, xo16=0.5, xne20=0, xna23=0, xmg24=0, other_args=None
):
    """
    Runs the self_heat.py script to simulate nuclear heating and neutrino losses.

    Parameters:
        dt (float): Time duration for the simulation.
        rho (float): Density in g/cm^3.
        T (float): Temperature in K (default 1e9 K).
        xp, xhe4, xc12, xo16, xne20, xna23, xmg24 (floats): Abundances of species.
        other_args (list): Additional command-line arguments.

    Returns:
        tuple: Time-averaged nuclear energy generation, neutrino loss, heat capacity,
               final abundances, and effective simulation duration.
    """
    script_path = 'self_heat.py'
    work_dir = './self_heating_network'  # directory containing helm_table.dat and output files

    # Build command to run external script
    cmd = [
        'python3', script_path,
        '-rho', str(rho),
        '-T', str(T),
        '-abundance_list', str(xp), str(xhe4), str(xc12), str(xo16), str(xne20), str(xna23), str(xmg24),
        '-tmax', str(dt)
    ]

    if other_args is not None:
        cmd.extend(other_args)

    print(cmd)

    # Run the simulation
    result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
    if result.returncode != 0:
        print("Error output from self_heat.py:\n", result.stderr)
        print("STDOUT:\n", result.stdout)
        raise RuntimeError("self_heat.py failed")

    # Output file produced by self_heat.py (must match its naming convention)
    output_file = f"self_heating_network/self_heating_evolution{rho:.4e}-{T:.2e}.dat"

    # Read and parse output data
    with open(output_file, 'r') as f:
        lines = f.readlines()
    data_lines = [line for line in lines if not line.startswith("#") and line.strip()]

    # Initialize sums for time averaging
    sum_eps_nuc = sum_eps_nu = sum_cv = 0.0
    E_nuc, E_nu, times = [], [], []
    count = 0

    # Extract quantities from data
    for line in data_lines:
        parts = line.split()
        if len(parts) < 7:
            continue
        time_now = float(parts[0])
        eps_nuc = float(parts[3])
        eps_nu = float(parts[4])
        cv = float(parts[6])

        E_nuc.append(eps_nuc)
        E_nu.append(eps_nu)
        times.append(time_now)

        sum_eps_nuc += eps_nuc
        sum_eps_nu += eps_nu
        sum_cv += cv
        count += 1

    # Average values
    avg_eps_nuc = sum_eps_nuc / count if count > 0 else 0.0
    avg_eps_nu = sum_eps_nu / count if count > 0 else 0.0
    avg_cv = sum_cv / count if count > 0 else 0.0

    # Time-weighted averages
    E_nuc_np = np.array(E_nuc)
    E_nu_np = np.array(E_nu)
    dt_np = np.array(times)
    dt_steps = np.diff(dt_np, prepend=0.0)
    E_nuc_time_avg = np.sum(E_nuc_np * dt_steps) / (dt_np[-1] if len(dt_np) > 0 else 1)
    E_nu_time_avg = np.sum(E_nu_np * dt_steps) / (times[-1] if len(times) > 0 else 1)

    # Read final abundances from abundance file
    output_file2 = "self_heating_network/abundances_vs_time.dat"
    x_p = x_he4 = x_c12 = x_o16 = x_ne20 = x_na23 = x_mg24 = 0
    with open(output_file2, 'r') as f:
        lines = f.readlines()
        data_lines = [line for line in lines if not line.startswith("#") and line.strip()]
        if data_lines:
            last_line = data_lines[-1].split()
            x_p, x_he4, x_c12, x_o16, x_ne20, x_na23, x_mg24 = map(float, last_line[1:8])

    return (
        E_nuc_time_avg, E_nu_time_avg, avg_cv,
        x_p, x_he4, x_c12, x_o16, x_ne20, x_na23, x_mg24,
        times[-1] - dt if times else 0.0
    )

def run_self_heat_classical(
    dt, rho, T=1e9, xp=0, xhe4=0, xc12=0.5, xo16=0.5, xne20=0, xna23=0, xmg24=0, epsilon_gravity=0, other_args=None
):
    """
    Runs the self_heat_classical.py script for classical integration mode.

    Returns:
        tuple: Initial temperature, heat capacity, and pressure from self_heat_classical doesn't evolve temperatutre.
    """
    script_path = 'self_heat_classical.py'
    work_dir = './self_heating_network'

    cmd = [
        'python3', script_path,
        '-rho', str(rho),
        '-T', str(T),
        '-abundance_list', str(xp), str(xhe4), str(xc12), str(xo16), str(xne20), str(xna23), str(xmg24),
        '-tmax', str(dt),
        '-epsilon_gravity', str(epsilon_gravity)
    ]

    print(cmd)

    # Run the classical simulation
    result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
    if result.returncode != 0:
        print("Error output from self_heat_classical.py:\n", result.stderr)
        print("STDOUT:\n", result.stdout)
        raise RuntimeError("self_heat_classical.py failed")

    # Extract results from stdout (last few lines contain needed values)
    lines = result.stdout.split("\n")
    cv = lines[-4]
    pressure = lines[-3]
    T_initial = lines[-2]

    return T_initial, cv, pressure
