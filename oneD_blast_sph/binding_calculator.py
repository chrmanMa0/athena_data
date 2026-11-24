import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.integrate import cumulative_trapezoid, simpson
filename = "s18/s18.0"
data_frame = pd.read_fwf(filename, skprows=1)
data = data_frame.to_numpy()

G = 6.6743e-8
M_sun = 1.989e33
M_cut = 1.2 * M_sun

r = data[:,2]
rho = data[:,4]
integrand_M = 4 * np.pi * r**2 * rho
M_r = cumulative_trapezoid(integrand_M, r, initial=0)

envelope_mask = M_r >= M_cut

r_env = r[envelope_mask]
rho_env = rho[envelope_mask]
M_env = M_r[envelope_mask]

integrand_E = 4 * np.pi * G * M_env * rho_env * r_env
E_binding_envelope = simpson(integrand_E, x=r_env)

print(f"Envelope Base Radius (r_env): {r_env[0]:.4e} cm")
print(f"Total Envelope Mass {(M_r[-1] - M_cut)/M_sun: .2f} M_sun")
print(f"Gravitational Binding Energy: {E_binding_envelope:.4} ergs")
