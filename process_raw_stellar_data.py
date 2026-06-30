import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import numpy as np


def construct_enclosed_mass(radius_data, density_data):
    enclosed_mass_data = []
    for i in range(len(radius_data)):
        enclosed_mass = 0
        for j in range(i):
            enclosed_mass += 4.0/3.0 * np.pi * (radius_data[j+1]**3 - radius_data[j]**3) * density_data[j]
        enclosed_mass_data.append(enclosed_mass)
    return np.array(enclosed_mass_data)


initial_data_filename = "s20.0_stripped_and_cored.txt"
df = pd.read_csv(initial_data_filename, sep=r'\s+')
print(df.iloc[:, 2])

radius_data = df.iloc[:,2].to_numpy()
density_data = df.iloc[:, 4].to_numpy()
pressure_data = df.iloc[:, 6].to_numpy()
enclosed_mass_data = construct_enclosed_mass(radius_data, density_data)
full_data = np.column_stack((radius_data, density_data, pressure_data, enclosed_mass_data))

np.savetxt("data_filename.txt", full_data, delimiter=" ")
x_vals = np.linspace(radius_data[0], radius_data[-1], 4000)
cs = CubicSpline(radius_data, density_data, bc_type='natural')
y_vals = cs(x_vals)
plt.yscale("log")
plt.xscale("log")
plt.plot(x_vals, y_vals)
plt.show()
