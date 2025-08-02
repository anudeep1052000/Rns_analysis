# Corrected import
import os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

# Mock EOS class for Helmholtz tabulated EOS (needs real implementation)
class HelmholtzEOS:
    def __init__(self, filepath):
        self.rho_vals, self.eps_vals, self.P_vals = [], [], []
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 4:
                    self.rho_vals.append(float(parts[1]))
                    self.eps_vals.append(float(parts[2]))
                    self.P_vals.append(float(parts[3]))
        self.rho_vals = np.array(self.rho_vals)
        self.P_vals = np.array(self.P_vals)
        self.eps_vals = np.array(self.eps_vals)

        self.eps_of_P = interp1d(self.P_vals, self.eps_vals, kind='cubic', fill_value='extrapolate')
        self.rho_of_P = interp1d(self.P_vals, self.rho_vals, kind='cubic', fill_value='extrapolate')

    def pressure(self, rho):  # Inverse interpolation if needed
        return np.interp(rho, self.rho_vals, self.P_vals)

    def energy_density(self, P):
        return self.eps_of_P(P)

    def get_interp_dict(self):
        return {'eps': self.energy_density}

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
"""
SMAX =0.9999
SDIV = 129
s_gp = []
for s in range(1,SDIV+1):
    s_gp[s] = SMAX*(s-1.0)/(SDIV-1.0)
""" 
def tov_rhs(r, y, eos_interp):
    m, P, nu = y
    if P <= 0:
        return [0, 0, 0]
    eps = eos_interp['eps'](P)
    dmdr = 4 * np.pi * r**2 * eps
    dPdr = -(eps + P) * (m + 4*np.pi*r**3*P) / (r * (r - 2*m))
    dnudr = 2 * (m + 4*np.pi*r**3*P) / (r * (r - 2*m))
    return [dmdr, dPdr, dnudr]

def integrate_tov(e_center, P_c, eos_interp, r_max=2e6):
    y0 = [4/3 * np.pi * eos_interp['eps'](P_c) * 1e-6**3, P_c, 0.0]
    sol = solve_ivp(
        tov_rhs, [1e-6, r_max], y0, args=(eos_interp,), 
        events=lambda r, y,eos: y[1], dense_output=True, max_step=100
    )
    r_arr = sol.t
    m_arr, P_arr, nu_arr = sol.y
    return r_arr, m_arr, P_arr, nu_arr, r_arr[-1], m_arr[-1]


def sphere(s_gp, eos_interp, e_center, P_c, P_surface, e_surface, white_dwarf_mode=False):
    SDIV = len(s_gp) - 1
    RDIV = 501

    # Integrate TOV three times (placeholder behavior, can refine later)
    r_arr, m_arr, P_arr, nu_arr, r_final, m_final = integrate_tov(e_center, P_c, eos_interp)

    # Define r_is_gp and lambda/nu interpolation arrays
    r_is_gp = np.linspace(1e-6, r_final, RDIV)
    lambda_gp = np.log((1 + m_final / (2*r_is_gp))**2)
    nu_gp_func = interp1d(r_arr, nu_arr, kind='cubic', fill_value='extrapolate')
    nu_gp = nu_gp_func(r_is_gp)

    # Prepare metric fields
    gama = np.zeros((SDIV+1, SDIV+1))
    rho = np.zeros((SDIV+1, SDIV+1))
    alpha = np.zeros((SDIV+1, SDIV+1))
    omega = np.zeros((SDIV+1, SDIV+1))  # zero in non-rotating case

    gama_mu_0 = np.zeros(SDIV+1)
    rho_mu_0 = np.zeros(SDIV+1)

    for s in range(1, SDIV+1):
        r_is_s = r_final * (s_gp[s] / (1.0 - s_gp[s]))
        if r_is_s < r_final:
            lambda_s = np.interp(r_is_s, r_is_gp, lambda_gp)
            nu_s = np.interp(r_is_s, r_is_gp, nu_gp)
        else:
            lambda_s = 2 * np.log(1 + m_final / (2*r_is_s))
            nu_s = np.log((1 - m_final/(2*r_is_s)) / (1 + m_final/(2*r_is_s)))

        gama_s = nu_s + lambda_s
        rho_s = nu_s - lambda_s

        gama[s, 1] = gama_s
        rho[s, 1] = rho_s
        for m in range(1, SDIV+1):
            gama[s, m] = gama_s
            rho[s, m] = rho_s
            alpha[s, m] = 0.5 * (gama_s - rho_s)

        gama_mu_0[s] = gama_s
        rho_mu_0[s] = rho_s

    # Compute isotropic radius at equator
    s_e = 0.5
    gama_eq = np.interp(s_e, s_gp[1:], gama_mu_0[1:])
    rho_eq = np.interp(s_e, s_gp[1:], rho_mu_0[1:])
    r_e = r_final * np.exp(0.5 * (rho_eq - gama_eq))

    return gama, rho, alpha, omega, r_e



# Extract values from result file
def extract_nonrotating_star_info(filename):
    central_pressure = None
    first_model_line = None
    with open(filename, 'r') as f:
        for line in f:
            if "DEBUG: Central pressure computed:" in line:
                parts = line.strip().split(':')[-1].strip().split('=')
                if len(parts) == 2:
                    central_pressure = float(parts[1])
            if first_model_line is None and line.strip().startswith("1.000"):
                first_model_line = line.strip()
    return central_pressure, first_model_line

# Main routine
def compute_star_from_helmholtz(directory, eos_filename,result_directory, result_filename, rho_c):
    eos_path = os.path.join(directory, eos_filename)
    result_path = os.path.join(result_directory, result_filename)
    eos = HelmholtzEOS(eos_path)

    # Calculate central pressure and energy density
    P_c = eos.pressure(rho_c)           # Pressure at central density
    e_center = eos.energy_density(P_c)  # Energy density at central pressure
    
    # Assume surface pressure and energy density near zero (vacuum)
    P_surface = 0.0
    e_surface = 0.0
    
    # Create interpolation dict for EOS (needed by sphere)
    eos_interp = {'eps': eos.energy_density, 'pressure': eos.pressure}
    
    # Build s_gp array for sphere (SDIV+1 points from 0 to ~1)
    SDIV = 129
    SMAX = 0.9999
    s_gp = np.zeros(SDIV+1)
    for s in range(1, SDIV+1):
        s_gp[s] = SMAX * (s - 1) / (SDIV - 1)
    
    # Call sphere with correct parameters
    gama, rho, alpha, omega, r_e = sphere(
        s_gp=s_gp,
        eos_interp=eos_interp,
        e_center=e_center,
        P_c=P_c,
        P_surface=P_surface,
        e_surface=e_surface,
        white_dwarf_mode=True
    )
    
    # Extract central pressure and first line from result file
    central_pressure, model_line = extract_nonrotating_star_info(result_path)
    
    # NOTE: You also need to get radius and mass from the TOV integration if you want to return them
    # You can integrate TOV here or modify sphere() to return radius and mass if needed.
    # For now, let's just set placeholders:
    radius = r_e   # isotropic equatorial radius from sphere()
    mass = None    # Not computed here
    print(r_e)
    return {
        "radius_cm": radius,
        "mass_g": mass,
        "central_pressure": central_pressure,
        "model_line": model_line,
    }


# Example usage
directory = "/home/anudeep/rns_wd_backup/rns_wd-update1/source/eos"
eos_filename = "wd_eos_helmholtz.dat"
result_directory = "/home/anudeep/rns_wd_backup/rns_wd-update1/source/rns.v2.0/outputs2"
result_filename = "rho_1.500e+09.txt"
rho_c = 1.5e9

compute_star_from_helmholtz(directory, eos_filename,result_directory, result_filename, rho_c)
