
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# load columns: time(s) ... t_heat ... t_dyn  (adjust indices if file different)
data = np.loadtxt("stellar_evolution_output.dat", comments="#")
time = data[:,0]
t_heat = data[:,7]    # column index with t_heat (0-based as in your header)
t_dyn  = data[:,8]    # t_dyn column

# use only rows where t_heat is finite and positive
mask = (t_heat > 0) & np.isfinite(t_heat)
time_fit = time[mask]
logth = np.log10(t_heat[mask])

# linear fit log10(t_heat) = m * time + b
m, b, r, p, se = stats.linregress(time_fit, logth)

# solve for t where t_heat = t_dyn (use median t_dyn or last value)
t_dyn_val = np.median(t_dyn)   # or t_dyn[-1]
logt_target = np.log10(t_dyn_val)
if m != 0:
    t_cross = (logt_target - b) / m
else:
    t_cross = np.inf

print("fit slope m =", m, "intercept b =", b)
print("predicted cross time (s):", t_cross)

# plot
plt.figure()
plt.plot(time, t_heat, '.-', label='t_heat (data)')
plt.plot(time_fit, 10**(m*time_fit + b), '--', label='fit')
plt.axhline(t_dyn_val, color='k', lw=1, label=f"t_dyn={t_dyn_val:.2e}s")
if np.isfinite(t_cross):
    plt.axvline(t_cross, color='r', lw=1, label=f"predicted cross {t_cross:.2e}s")
plt.yscale('log')
plt.xlabel('time (s)')
plt.ylabel('t_heat (s)')
plt.legend()
plt.savefig("Trend.png")
#plt.show()
