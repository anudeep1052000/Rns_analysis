import math
import numpy as np
from call_self_heat import run_self_heat

# Physical constants
M_sun = 1.989e33          # Solar mass in grams
C = 2.99792458e10         # Speed of light in cm/s
G = 6.67430e-8            # Gravitational constant in cm^3/g/s^2
seconds_in_year = 365.25 * 24 * 3600

mass_tolerance = 1e-4


def compute_breakup_omega(mass, radius):
    """
    Compute the breakup angular velocity (rad/s) for a star.
    """
    M = mass * M_sun
    return math.sqrt(G * M / (radius * 1e5)**3)


def drho_dt(mass_slope, mass_dj_dt, e_nuc, e_neutrino):
    """
    Compute d(rho)/dt for each mass sequence.
    """
    mass_drho_dt = {}
    for mass_key in mass_slope.keys():
        slope_drho_dt = []
        for j in range(len(mass_slope[mass_key])):
            drho = mass_slope[mass_key][j] * mass_dj_dt[mass_key][j]
            slope_drho_dt.append(drho)
        mass_drho_dt[mass_key] = slope_drho_dt
    return mass_drho_dt


def guess_timelist(mass_drho_dt, rho_mass):
    """
    Estimate an evolutionary time list based on density change rates.
    """
    mass_time_list = {}
    for mass_key in mass_drho_dt.keys():
        time_list = [0.0]
        for j in range(1, len(mass_drho_dt[mass_key])):
            # Avoid division by zero or near-zero derivatives
            if abs(mass_drho_dt[mass_key][j]) < 1e-30:
                dt = 1e3
            else:
                dt = (rho_mass[mass_key][j] - rho_mass[mass_key][j-1]) / \
                     (mass_drho_dt[mass_key][j] - mass_drho_dt[mass_key][j-1])

            dt = abs(dt) if dt != 0 else 1e3
            time_list.append(time_list[-1] + dt)

        mass_time_list[mass_key] = time_list
    return mass_time_list


def remove_near_duplicates_sorted(rho, J, radius, omega, pc, tol=1e4):
    """
    Remove near-duplicate rho values from sorted arrays by clustering points
    within 'tol' and keeping representative values.
    """
    reduced_rho, reduced_J, reduced_radius, reduced_omega, reduced_pc = [], [], [], [], []
    cluster_rho, cluster_J, cluster_radius, cluster_omega, cluster_pc = \
        [rho[0]], [J[0]], [radius[0]], [omega[0]], [pc[0]]

    for i in range(1, len(rho)):
        if abs(rho[i] - cluster_rho[-1]) <= tol:
            cluster_rho.append(rho[i])
            cluster_J.append(J[i])
            cluster_radius.append(radius[i])
            cluster_omega.append(omega[i])
            cluster_pc.append(pc[i])
        else:
            reduced_rho.append(np.mean(cluster_rho))
            reduced_J.append(np.max(cluster_J))
            reduced_radius.append(np.max(cluster_radius))
            reduced_omega.append(np.max(cluster_omega))
            reduced_pc.append(np.max(cluster_pc))
            cluster_rho, cluster_J, cluster_radius, cluster_omega, cluster_pc = \
                [rho[i]], [J[i]], [radius[i]], [omega[i]], [pc[i]]

    # Save final cluster
    reduced_rho.append(np.mean(cluster_rho))
    reduced_J.append(np.max(cluster_J))
    reduced_radius.append(np.max(cluster_radius))
    reduced_omega.append(np.max(cluster_omega))
    reduced_pc.append(np.max(cluster_pc))

    return (np.array(reduced_rho), np.array(reduced_J),
            np.array(reduced_radius), np.array(reduced_omega), np.array(reduced_pc))


def epsilon_grav(mass_drho_dt, mass_pc, rho_mass):
    """
    Compute gravitational energy release rate.
    """
    mass_E_gravity = {}
    for mass_key in mass_drho_dt.keys():
        E_gravity = (np.array(mass_drho_dt[mass_key]) * np.array(mass_pc[mass_key])) / \
                    np.array(rho_mass[mass_key])**2
        mass_E_gravity[mass_key] = E_gravity
    return mass_E_gravity


def epsilon_nuc_nu(rho, mass):
    """
    Compute nuclear and neutrino energy generation rates.
    """
    eps_nuc, eps_nu, cv = run_self_heat(rho, 1e9)
    print(eps_nuc, eps_nu, cv)
    print(f"Mass = {mass}, rho = {rho} successfully done")
    return eps_nuc, eps_nu, cv


def evolution(mass_E_gravity, mass_rho):
    """
    Compute temperature evolution rate for each mass sequence.
    """
    mass_Temp_dt = {}
    for mass_key in mass_E_gravity.keys():
        dTemp_dt = []
        for j, rho_value in enumerate(mass_rho[mass_key]):
            eps_nuc, eps_nu, cv = epsilon_nuc_nu(rho_value, mass_key)
            dTemp_dt.append((eps_nuc - eps_nu + mass_E_gravity[mass_key][j]) / cv)
        mass_Temp_dt[mass_key] = dTemp_dt
    return mass_Temp_dt


def check_lengths_equal(dicts):
    """
    Check if all dictionaries have lists of the same length for each mass key.
    """
    keys = dicts[0].keys()
    for key in keys:
        lengths = [len(d[key]) for d in dicts if key in d]
        if len(set(lengths)) != 1:
            print(f"❌ Mismatch at mass={key}: lengths = {lengths}")
            return False
    print("✅ All dictionaries have lists of equal length for each mass key.")
    return True

