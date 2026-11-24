#!/usr/bin/python3

import pandas as pd
import numpy as np


data_pd = pd.read_csv("Blast.block0.out3.13300-cgs.dat", sep=r"\s+", skiprows=2, header=None)
integral_start = 400
radius = data_pd.iloc[:,1]
dr = np.gradient(radius)
dV = 4 * np.pi * (radius**2) * dr
mass = dV*data_pd.iloc[:, 2]
ke = .5 * mass * data_pd.iloc[:, 4]**2
total_ke = np.sum(ke[184:])
print(total_ke)
