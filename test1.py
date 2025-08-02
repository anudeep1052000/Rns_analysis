M_sun = 1.989e33
mass_values = []  # Array to store mass values
spin_values = []  # Array to store spin values
momentum_values = []
    # Open the output file and read the contents
with open('rho_2.844e+09.txt', 'r') as file:
    lines = file.readlines()

    # Find the line where the table starts (looking for 'ratio' in the header)
start_idx = 0
for i, line in enumerate(lines):
    if line.startswith('ratio'):
        start_idx = i + 1  # Skip the header line
        break

    # Iterate over the lines after the header to extract mass and spin
for line in lines[start_idx:]:
        # Skip empty lines or lines that don't look like data rows
    if not line.strip():
        continue
    parts = line.split()
# Ensure there are enough columns (should have at least 7 c
    if len(parts) < 7:
        continue

    try:
            # Extract mass (3rd column) and spin (6th column)
        mass = float(parts[2])  # Mass is in the 3rd column (index 2)
        spin = float(parts[5])
        amomentum = float(parts[8])
        print(amomentum)  # Spin is in the 6th column (index 5)
        J = amomentum * (mass * M_sun)**2
        print(J)
            # Append values to respective arrays
        mass_values.append(mass)
        spin_values.append(spin)
        momentum_values.append(amomentum * (mass * M_sun)**2)

    except:
          print("error")
