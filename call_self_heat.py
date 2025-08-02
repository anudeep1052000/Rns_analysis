import subprocess
import numpy as np

def run_self_heat(rho, T=1000000000.0, other_args=None):
    script_path = 'self_heat.py'
    work_dir = './self_heating_network'  # where helm_table.dat should be
    cmd = ['python3', script_path, '-rho', str(rho),'-xhe4',str(0),'-xc12',str(0.5),'-xo16',str(0.5)]
    
    result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)

    if result.returncode != 0:
        print("Error output from self_heat.py:\n", result.stderr)
        print("STDOUT:\n", result.stdout)
        raise RuntimeError("self_heat.py failed")
        



    if other_args is not None:
        cmd.extend(other_args)

    # Run the command, wait for it to finish
    
    # After it runs, read the output file (make sure self_heat.py writes a predictable filename)
    #output_file = f"/home/anudeep/self_heating_network/self_heating_evolution{rho:.4e}-{T:.2e}.dat"
    output_file = f"self_heating_network/self_heating_evolution{rho:.4e}-{T:.2e}.dat"

    # Parse the output file to get the last line or the relevant line(s)
    with open(output_file, 'r') as f:
        lines = f.readlines()
    # Filter out comment lines and empty lines
        data_lines = [line for line in lines if not line.startswith("#") and line.strip()]

# Initialize sums and count
    sum_eps_nuc = 0.0
    sum_eps_nu = 0.0
    sum_cv = 0.0
    count = 0

    for line in data_lines:
        parts = line.split()
    # Make sure line has enough columns
        if len(parts) < 7:
            continue
        eps_nuc = float(parts[3])
        eps_nu = float(parts[4])
        cv = float(parts[6])

        sum_eps_nuc += eps_nuc
        sum_eps_nu += eps_nu
        sum_cv += cv
        count += 1

    if count > 0:
        avg_eps_nuc = sum_eps_nuc / count
        avg_eps_nu = sum_eps_nu / count
        avg_cv = sum_cv / count
    else:
    # Handle empty file or no data case
        avg_eps_nuc = avg_eps_nu = avg_cv = 0.0

    return avg_eps_nuc, avg_eps_nu, avg_cv


#eps_nuc ,eps_nu,cv = run_self_heat(2.5e9,1e9)
#print(eps_nuc,eps_nu,cv)
