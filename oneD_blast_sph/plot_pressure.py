import numpy as np
from matplotlib import pyplot as plt
import argparse
from scipy.integrate import cumulative_trapezoid

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



day_s = 86400

def getDataAndTime(filename):
    with open(filename, 'r') as f:
        time_line_list = f.readline().strip().split(" ")
        #time_string = time_line_list[1] + " " + str(float(time_line_list[2])/day_s)
        time_string = time_line_list[1] + " " + time_line_list[2]

        data_array = np.loadtxt(f, skiprows=1, dtype=float)
        return data_array, time_string
parser = argparse.ArgumentParser(description="Plots data from columns x, y")
parser.add_argument("-x", type=int, default=1)
parser.add_argument("-y", type=int, default=3)
args = parser.parse_args()
def calc_entropy_data(data_matrix, gamma):
    r = data_matrix[:, 1]
    rho = data_matrix[:, 2]
    integrand = 4 * np.pi * (r**2)*rho
    enclosed_mass = cumulative_trapezoid(integrand, x=r, initial=0)
    P = data_matrix[:, 3]
    m_B = 1.6726219e-24
    k_B = 1.380649e-16
    a = 7.565733e-15
    multiplier = (4.0 *m_B / (3.0 *k_B))*(3.0 *a)**(0.25)
    entropy_per_baryon = multiplier * (P**0.75) / rho
    return (enclosed_mass, entropy_per_baryon, rho)

if(args.y == 6):
    data_matrix, time_string = getDataAndTime("s12/Blast.block0.out3.00000-cgs.dat")
    plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=rf"$12M:${time_string:.9} (s)")
    data_matrix, time_string = getDataAndTime("s15/Blast.block0.out3.00002-cgs.dat")
    plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=time_string)
    data_matrix, time_string = getDataAndTime("s18/Blast.block0.out3.00001-cgs.dat")
    plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=time_string)
###data_matrix, time_string = getDataAndTime("Blast.block0.out3.00001-cgs.dat")
#plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=time_string)
#data_matrix, time_string = getDataAndTime("Blast.block0.out3.00005-cgs.dat")
#plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=time_string)
#data_matrix, time_string = getDataAndTime("Blast.block0.out3.00010-cgs.dat")
#plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=time_string)
#data_matrix, time_string = getDataAndTime("Blast.block0.out3.01000-cgs.dat")
#plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=time_string)
#data_matrix, time_string = getDataAndTime("Blast.block0.out3.02000-cgs.dat")
#plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=time_string)
#data_matrix, time_string = getDataAndTime("Blast.block0.out3.03000-cgs.dat")
#plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=time_string)
#data_matrix, time_string = getDataAndTime("Blast.block0.out3.00500-cgs.dat")
#plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=time_string)
#data_matrix, time_string = getDataAndTime("Blast.block0.out3.01000-cgs.dat")
#plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=time_string)
#data_matrix, time_string = getDataAndTime("Blast.block0.out3.02000-cgs.dat")
#plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=time_string)
if(args.y != 7):
    data_matrix, time_string = getDataAndTime("s12/Blast.block0.out3.00142-cgs.dat")
    plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=rf"$12M:$ {time_string:.12} (s)")
    data_matrix, time_string = getDataAndTime("s15/Blast.block0.out3.00075-cgs.dat")
    plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=rf"$15M:$ {time_string:.12} (s)")
    data_matrix, time_string = getDataAndTime("s18/Blast.block0.out3.00193-cgs.dat")
    plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=rf"$18M:$ {time_string:.12} (s)")
else:
    #need to plot entropy vs enclosed mass
    data_matrix, time_string = getDataAndTime("s12/Blast.block0.out3.01010-cgs.dat")
    enclosed_mass, entropy, rho = calc_entropy_data(data_matrix, 4.0/3.0)
    plt.plot(enclosed_mass, entropy, label=rf"$12M:$ {time_string:.12} (s)")

    data_matrix, time_string = getDataAndTime("s15/Blast.block0.out3.00556-cgs.dat")
    enclosed_mass, entropy, rho = calc_entropy_data(data_matrix, 4.0/3.0)
    plt.plot(enclosed_mass, rho, label=rf"$15M:$ {time_string:.12} (s)")
    #plt.plot(data_matrix[:, args.x], rho, label=rf"$15M:$ {time_string:.12} (s)")

    data_matrix, time_string = getDataAndTime("s18/Blast.block0.out3.00193-cgs.dat")
    enclosed_mass, entropy, rho = calc_entropy_data(data_matrix, 4.0/3.0)
    plt.plot(enclosed_mass, entropy, label=rf"$18M:$ {time_string:.12} (s)")


    #plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=rf"$18M:$ {time_string:.12} (s)")


#data_matrix, time_string = getDataAndTime("s12/Blast.block0.out3.12000-cgs.dat")
#plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=rf"$12M:$ {time_string:.12} (s)")
#data_matrix, time_string = getDataAndTime("s15/Blast.block0.out3.12000-cgs.dat")
#plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=rf"$15M:$ {time_string:.12} (s)")
#data_matrix, time_string = getDataAndTime("s18/Blast.block0.out3.10770-cgs.dat")
#plt.plot(data_matrix[:, args.x], data_matrix[:, args.y], label=rf"$18M:$ {time_string:.12} (s)")










plt.xlabel("Radius (cm)")
if(args.y == 3):
    plt.ylabel("Pressure (dynes)")
    plt.ylim(1, 1e27)
    plt.yscale("log")
elif(args.y == 2):
    plt.ylabel(r"Density ($\frac{g}{cm^3}$)")
    plt.ylim(1e-10, 10)
    #plt.xlim(1, 1e14)
    plt.yscale("log")
elif(args.y == 4):
    plt.ylabel(r"Velocity ($\frac{cm}{s}$)")
    plt.ylim(1e7,3.0e9)
elif(args.y == 7):
    plt.ylabel(r"Density")
    plt.xlabel("Enclosed Mass (g)")
    plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.show()
