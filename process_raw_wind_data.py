import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import numpy as np
import os
os.environ["NUMEXPR_MAX_THREADS"] = "272"  # Sets the max allowance high enough
os.environ["NUMEXPR_NUM_THREADS"] = "8"    # Restricts actual usage so it stays stable
#We will need to read in the data from Potsdam site
#https://www.astro.physik.uni-potsdam.de/PoWR/powrgrid1.php#wrgridstart
# For wolf-rayet data

a = 7.5657e-15

def construct_enclosed_mass(radius_data, density_data):
    enclosed_mass_data = []
    for i in range(len(radius_data)):
        enclosed_mass = 0
        for j in range(i):
            enclosed_mass += 4.0/3.0 * np.pi * (radius_data[j+1]**3 - radius_data[j]**3) * density_data[j]
        enclosed_mass_data.append(enclosed_mass)
    return np.array(enclosed_mass_data)


initial_data_filename = "wne_09-12_stratification.txt"
#[0] Depth Index
#[1] R_1
#[2] Log(R_1)
#[3] Temp (K)
#[4] Log Particle Density
#[5] Log Electron Density
#[6] Tau Thompson (NLTE)
#[7] Tau Rosseland (NLTE)
#[8] Tau Rosselandcont
#[9] Velocity (km/s)
#[10] Gradient
#[11] Clummping
#[12] Mean Partcle Mass

df = pd.read_csv(initial_data_filename, engine="python", sep=r'\s+', skiprows=6, skipinitialspace=True)
print(df.iloc[:, 2])

radius_data = np.flip(df.iloc[:,1].to_numpy())
density_data = np.flip(df.iloc[:, 12].to_numpy()) * 10**np.flip(df.iloc[:, 4])
pressure_data = (1.0/3.0)*a*np.power(np.flip(df.iloc[:, 3].to_numpy()), 4) #Using radiation dominated temperature
#enclosed_mass_data = construct_enclosed_mass(radius_data, density_data)
full_data = np.column_stack((radius_data, density_data, pressure_data)) #, enclosed_mass_data))
print(pressure_data)
np.savetxt("wind_data_filename.txt", full_data, delimiter=" ")
x_vals = np.linspace(radius_data[0], radius_data[-1], 4000)
cs = CubicSpline(radius_data, pressure_data, bc_type='natural')
y_vals = cs(x_vals)
plt.yscale("log")
plt.xscale("log")
plt.plot(x_vals, y_vals)
plt.show()
