import subprocess
import numpy as np

def run_self_heat(dt,rho, T=1000000000.0,xp=0,xhe4=0,xc12=0.5,xo16=0.5,xne20=0,xna23=0,xmg24=0,other_args=None):
    script_path = 'self_heat.py'
    work_dir = './self_heating_network'  # where helm_table.dat should be
    cmd = [
    'python3', script_path,
    '-rho', str(rho),
    '-T', str(T),
    '-abundance_list',str(xp), str(xhe4) ,str(xc12) ,str(xo16), str(xne20),str(xna23),str(xmg24),
    '-tmax', str(dt)
]
    print(cmd)
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
    x_c12 =0
    x_o16 = 0
    x_he4 = 0
    x_p = 0
    x_ne20 = 0
    x_na23 = 0
    x_mg24 = 0

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
    output_file2 = "self_heating_network/abundances_vs_time.dat"
    with open(output_file2,'r') as f:
        lines = f.readlines()
        data_lines = [line for line in lines if not line.startswith("#") and line.strip()]
        if data_lines:
            last_line = data_lines[-1]
            parts = last_line.split()
            #if len(parts) < 8:
            x_p = float(parts[1])
            x_he4 = float(parts[2])
            
            x_c12 = float(parts[3])
            x_o16 = float(parts[4])
            x_ne20 = float(parts[5])
            x_na23 = float(parts[6])
            x_mg24 = float(parts[7])
            
            # Do your processing here
            

    return avg_eps_nuc, avg_eps_nu, avg_cv,x_p,x_he4,x_c12,x_o16,x_ne20,x_na23,x_mg24


#eps_nuc,eps_nu,cv,x_p,x_he4x_c12,x_o16,x_ne20,x_na23,x_mg24 = run_self_heat(1e5,1.5e9,1e9)
#print(eps_nuc,eps_nu,cv,x_p,x_he4x_c12,x_o16,x_ne20,x_na23,x_mg24)
