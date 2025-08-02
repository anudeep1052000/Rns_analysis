import re

G = 6.67430e-8  # Gravitational constant (cm^3 g^-1 s^-2)
M_sun = 1.989e33  # Solar mass in grams

def non_rel_virial_check(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("ratio") or "DEBUG" in line or "About" in line:
                continue
            parts = re.split(r'\s+', line)
            if len(parts) != 9:
                continue
            try:
                ratio = float(parts[0])
                M_sun_val = float(parts[2])
                M0_sun_val = float(parts[3])
                r_star_km = float(parts[4])
                Omega = float(parts[6])
                J_over_M2 = float(parts[8])
                I  = float(parts[7])

                M = M_sun_val * M_sun
                M0 = M0_sun_val * M_sun
                R = r_star_km * 1e5  # km to cm

                J = J_over_M2 * (M**2)
                T = 0.5 * I * Omega**2
                #T = 0.5*J*Omega
                W = -(3/5) * G * M**2 / R
                print("T = "+str(T))
                print("W = "+str(W))
                virial_error = abs(2*T + W) / abs(W)

                data.append((ratio, virial_error))
            except Exception as e:
                continue

    # Sort by ratio for better visualization
    data.sort(key=lambda x: x[0])
    return data

# Usage example
data = non_rel_virial_check('rho_9.969e+09.txt')
for ratio, verr in data:
    print(f"ratio={ratio:.3f}, non-rel virial error={verr:.3e}")
