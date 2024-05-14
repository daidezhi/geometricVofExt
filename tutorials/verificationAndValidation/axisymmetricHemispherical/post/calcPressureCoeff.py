"""
Description
    This script is designed to extract the dimensionless distance ($s/D$)
    and pressure coefficients ($c_p$) along the structure surface.

Developer
    Dezhi Dai, NSE, ANL
"""

import numpy as np
import sys

# Operating conditions --------------------------------------------------------
run_name = str(sys.argv[1]) ## Case name
sigma = float(sys.argv[2])  ## Cavitation number $\sigma$,      [-]
u_in = 5.45                 ## Inlet velocity $u_in$,           [m/s]
# -----------------------------------------------------------------------------


# Water properties ------------------------------------------------------------
p_v = 2339.0        ## Saturation pressure $p_v$,               [Pa]
rho_water = 998.0   ## Water density $\rho_{water}$,            [kg/m^3]
# -----------------------------------------------------------------------------


# Compute outlet pressure -----------------------------------------------------
p_dyn = 0.5 * rho_water * u_in * u_in   ## Dynamic pressure,    [Pa]
p_0 = sigma * p_dyn + p_v               ## Outlet pressure,     [Pa]
# -----------------------------------------------------------------------------


# Geometry info of the hemispherical structure
x_0 = 0.0       ## x coordinate of structure nose,              [m]
D = 0.025       ## Diameter,                                    [m]
R = 0.5 * D     ## Radius,                                      [m]
# -----------------------------------------------------------------------------


# Read numerical data ---------------------------------------------------------
raw_file = str(sys.argv[3]) ## Path of 'structureSurface_pMean.xy'
x, p_mean = np.loadtxt(raw_file, delimiter=None, unpack=True)
# -----------------------------------------------------------------------------


# Extract distance $s$ along the structure surface ----------------------------
s = []

for x_i in x:
    dist_i = abs(x_i - x_0)

    if dist_i <= R:
        theta_i = np.arccos((R - dist_i) / R)
        s.append(theta_i * R)
    else:
        s.append(0.5*np.pi*R + (dist_i-R))

s = np.array(s)
s_by_D = s / D      ## Dimensionless distance $s/D$
# -----------------------------------------------------------------------------


# Calculate pressure coefficients $c_p$ ---------------------------------------
c_p = []

for p_mean_i in p_mean:
    c_p.append((p_mean_i - p_0) / p_dyn)

c_p = np.array(c_p)
# -----------------------------------------------------------------------------


# Save numerical data to file -------------------------------------------------
np.savetxt(
    str(sys.argv[4])+'/'+run_name+'.dat',
    np.transpose([s_by_D, c_p]),
    comments='# ',
    header='s_by_D                 c_p'
)
# -----------------------------------------------------------------------------