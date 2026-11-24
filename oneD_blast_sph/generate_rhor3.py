import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt

# ApJ Formatting

import matplotlib as mpl

import cycler

fntsz = 8

mpl.rcParams['text.usetex'] = False  # FIXME

mpl.rcParams['font.size'] = fntsz

mpl.rcParams['lines.color'] = 'k'

mpl.rcParams['figure.dpi'] = 300

mpl.rcParams['font.family'] = 'serif'

mpl.rcParams['savefig.dpi'] = 600

mpl.rcParams['xtick.direction'] = 'in'

mpl.rcParams['ytick.direction'] = 'in'

mpl.rcParams['lines.markersize'] = 3.0

mpl.rcParams['lines.linewidth'] = 1.2

mpl.rcParams['mathtext.default'] = 'regular'

colorcycle = mpl.rcParams['axes.prop_cycle'].by_key()['color']

linestyles = ('-', (0, (10, 4)), (0, (3, 2)), (0, (10, 4, 4, 4)), (0, (10, 2, 2, 2, 2, 2)), (0, (2, 1)),

              (0, (10, 2, 1, 2, 10, 2, 5, 2)), (0, (10, 2, 10, 2, 3, 2)))

propcycle = cycler.cycler(color=colorcycle[:8], linestyle=linestyles)
filename = "s12/Blast.block0.out3.00000-cgs.dat"
data_frame = pd.read_csv(filename, header=None, skiprows=2, sep=r'\s+')
data = data_frame.to_numpy()
print(data_frame)

rho = data[:, 2]
r = data[:, 1]

r3 = r**3
rhor3 = rho*r3

plt.plot(r, rhor3, label=r"12$\mathrm{M}_{\odot}$")
plt.axvline(x=1.98e11, color='b', linestyle="--", label=r"12$\mathrm{M}_{\odot}$ He/H Layer Transition Location")

filename = "s15/Blast.block0.out3.00000-cgs.dat"

data_frame = pd.read_csv(filename, header=None, skiprows=2, sep=r'\s+')
data = data_frame.to_numpy()
print(data[:, 2])

rho = data[:, 2]
r = data[:, 1]
r3 = r**3
rhor3 = rho*r3

plt.xscale("log")
plt.yscale("log")
plt.plot(r, rhor3, label=r"15$\mathrm{M}_{\odot}$")
plt.axvline(x=2.14e11, color='r', linestyle="--", label=r"15$\mathrm{M}_{\odot}$ He/H Layer Transition Location")

filename = "s18/Blast.block0.out3.00000-cgs.dat"
data_frame = pd.read_csv(filename, header=None, skiprows=2, sep=r'\s+')
data = data_frame.to_numpy()
print(data_frame)

rho = data[:, 2]
r = data[:, 1]

r3 = r**3
rhor3 = rho*r3

plt.plot(r, rhor3, label=r"18$\mathrm{M}_{\odot}$")
plt.axvline(x=2.54e11, color='g', linestyle="--", label=r"18$\mathrm{M}_{\odot}$ He/H Layer Transition Location")
plt.xlabel("Radius (cm)")
plt.ylabel(r"$\rho r^3$ (g)")
plt.ylim( bottom=1e29 )
plt.legend()
plt.title(r"$\rho r^3$ vs Radius")



plt.show()
