from scipy.signal import savgol_filter
import numpy as np
def find_slope_savgol(rho_mass, J_mass, window_length=5, polyorder=2,mode='nearest'):
    mass_slope = {}
    for i in rho_mass.keys():
        a = np.array(rho_mass[i])
        b = np.array(J_mass[i])
        # Sort by independent variable
        sorted_idx = np.argsort(b)
        a_sorted = a[sorted_idx]
        b_sorted = b[sorted_idx]
        # Compute derivative using Savitzky-Golay
        if len(a_sorted) >= window_length:
            slope = robust_numerical_derivative(b_sorted,a_sorted)
            mass_slope[i] = slope

    return mass_slope
def robust_numerical_derivative(x, y, window=7, poly=2):
    """
    Compute dy/dx from x and y with:
    - Savitzky-Golay smoothing on y
    - Central + forward/backward differences
    """
    x = np.array(x)
    y = np.array(y)

    # Smooth y (optional but helps for noisy signals)
    if len(y) >= window:
        y_smooth = savgol_filter(y, window_length=window, polyorder=poly)
        x_smooth = savgol_filter(x,window_length=window,polyorder=poly)
    else:
        y_smooth = y.copy()
        x_smooth = x.copy()

    dy_dx = np.zeros_like(y_smooth)

    for i in range(len(y)):
        if i == 0:
            dy = y_smooth[i+1] - y_smooth[i]
            dx = x_smooth[i+1] - x_smooth[i]
        elif i == len(y) - 1:
            dy = y_smooth[i] - y_smooth[i-1]
            dx = x_smooth[i] - x_smooth[i-1]
        else:
            dy = y_smooth[i+1] - y_smooth[i-1]
            dx = x_smooth[i+1] - x_smooth[i-1]

        dy_dx[i] = dy / dx if dx != 0 else 0.0

    return dy_dx
