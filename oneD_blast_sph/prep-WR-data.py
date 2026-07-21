#! /usr/bin/python3
import os
os.environ["NUMEXPR_MAX_THREADS"] = "3232"

import pandas as pd
import numpy as np
import uuid

pd.set_option("compute.use_numexpr", False)





# This program will turn a 2sn end of life star and strip and compress the star to
# simulate a WR star formed from this size star.


# For a WRN Early
# X = .35 for bottom of hydrogen shell.

# For WRN Late
# X = .1 for bottom of hydrogren shell


# Useful constants

M_dot = 1.989e33
R_dot = 6.957e10
a = 7.5657e-15
h1_early = .35
h1_late = .1
G = 6.674e-8

#useful indices
h1_fraction_column = 15
cell_enclosed_mass = 1
cell_outer_radius = 2
cell_density = 4
cell_velocity = 3
cell_pressure = 6
cell_energy_density = 7
cell_temperature = 5
"""
Usefull info about the columns of the 2sn data files

First row is the column header list

Column# Desc.
0       Row number
1       Cell outer total mass
2       Cell outer radius
3       Cell outer velocity
4       Cell density
5       cell temperature
6       cell pressure
7       cell specific energy
8       cell specific entropy
9       cell angular velocity
10      cell A_bar
11      cell Y_e
12      stability
13      NETWORK
14      neutrons
15      h1
16      he3
...

"""





def strip_shell(full_star_np, h1_fraction_cutoff):
   """
   Strips the shell of the star data off at the h1_fraction_cutoff
   leaving only the interior core reflecting the results of the mass
   loss on the content of the star.  Note does not actually relax the
   system to form a realistic WR star, just strips the outer envelope

   Parameters:
       full_star_np (2D np array of 2sn formated data with no headers): The starting data for the star
       h1_fraction_cutoff (float): The fraction of the elements (X) that is made up of hydrogen.

       Returns:
           2D Numpy Array containing only those rows with less than or equal to the cutoff value
           total remaining mass

   """

   # In the 2sn data file column 15 is the hydrogen fraction
   h1_fraction_col = 15

   stripped_star_np = full_star_np[full_star_np[:, h1_fraction_col] <= h1_fraction_cutoff]
   print(len(stripped_star_np))
   return stripped_star_np, stripped_star_np[-1,cell_enclosed_mass], stripped_star_np[-1,cell_outer_radius]

def load_and_prep_data(filename):
    """
    Reads the file and produces a numpy datafile

    Parameters:
        filename (string): Name of the file that will be opened and processed

    Returns: A numpy array of the data in the file missing header row
    """
    data_frame = pd.read_fwf(filename, skiprows=1)
    return data_frame.to_numpy()


def contract_stripped_data(stripped_data, pre_stripped_radius, stripped_radius):
    """
    We will apply a homologous contraction to the data due to the change in radius.
    We will change density and then adjust the cell mass to match the density change
    for each cell

    Parameters:
        stripped_data (numpy array): The stripped but uncontracted data
        pre_stripped_radius (real): The original radius of the star
        stripped_radius (real): The new outer radius of the star

        Returns:
            The contracted data as a 2D numpy array

    """
    #We first calculate the ratio due to the change in radius
    R_ratio = pre_stripped_radius/stripped_radius
    #Now we calculate the new outer radius of each cell
    stripped_data[:,cell_outer_radius] = stripped_data[:, cell_outer_radius]* 1.0/R_ratio
    #Now we recalculate the density for each of the cells due to the new density
    stripped_data[:, 4] = stripped_data[:, 4] * R_ratio**3
    #Now we recalculate the mass of each cell based on new density and radius
    cell_mass = np.diff(stripped_data[:, cell_enclosed_mass])
    print(cell_mass)
    print(stripped_data[:, cell_outer_radius]/R_dot)


def collapse_inner_core(stripped_data, core_size):
    print(core_size*M_dot)
    without_core_np = stripped_data[stripped_data[:, cell_enclosed_mass] >= core_size*M_dot]
    without_core_np[:,1] = without_core_np[:, 1] - core_size*M_dot
    return (without_core_np, without_core_np[-1, 1])

def convert_to_dimensionless(stripped_data, R_star, M_star):
    """
    This function converts the data to be dimensionless program units, allowing
    us to use G as a value of 1 and making sure we don't run into floating point errors
    when the exponents get too big.
    """
    L_0 = R_star #we want the data to be based on the radius of the star
    M_0 = M_star #Same for the mass
    rho_0 = M_0/(L_0)**3

    v_0 = 1e8
    t_0 = L_0/v_0#(((L_0)**3)/(G*M_0))**.5
    P_0 = rho_0*v_0**2
    T_0 = (P_0/a)**.25
    G_0 = (G * rho_0*L_0**2)/v_0**2

    scale_values = (L_0, M_0, rho_0, t_0, v_0, G_0, P_0, T_0)

    stripped_data[:, cell_outer_radius] = stripped_data[:, cell_outer_radius]/L_0
    stripped_data[:, cell_enclosed_mass] = stripped_data[:, cell_enclosed_mass]/M_0
    stripped_data[:, cell_density] = stripped_data[:, cell_density]/rho_0
    stripped_data[:, cell_pressure] = stripped_data[:, cell_pressure]/(rho_0 * v_0**2)
    stripped_data[:, cell_energy_density] = stripped_data[:, cell_energy_density] / (rho_0 * v_0**2)
    stripped_data[:, cell_temperature] = stripped_data[:, cell_temperature] / T_0
    return stripped_data, scale_values

def store_data(filename, scale_values, stripped_data):
    """
    This function will store the two pieces of the input data in two separate files
    with matchin uuids.  One piece is the scale values used to generate the
    dimensionless values used in Athena++ and the other is the actual
    dimensionless input data.

    Parameters
    filename (str): The custom name of the file
    scale_values (np array): The scale values needed to regenerate the data
    stripped_data (np array): The scaled values that will be needed to run athena

    """
    my_uuid = uuid.uuid4()
    filename = filename + "-" + str(my_uuid)
    np.savetxt(filename+".txt", stripped_data, fmt="%f", delimiter=" ")
    np.savetxt(filename+"-scales.txt", scale_values)


def main():
    filename = "s25.0"
    full_data_np = load_and_prep_data(filename)
    pre_stripped_radius = full_data_np[-1, cell_outer_radius]
    stripped_data, remaining_mass, outer_radius = strip_shell(full_data_np, h1_early)
    print(remaining_mass/M_dot)
    print(outer_radius/R_dot)
    #contract_stripped_data(stripped_data, pre_stripped_radius, 2.0*R_dot)
    stripped_data, mass_no_core = collapse_inner_core(stripped_data, 1.2)
    print(stripped_data)
    scaled_data, scale_values = convert_to_dimensionless(stripped_data, outer_radius, mass_no_core)
    save_cols = [cell_outer_radius, cell_enclosed_mass, cell_density, cell_pressure, cell_energy_density, cell_temperature]
    print(str(scale_values))
    store_data("save_test", scale_values, scaled_data[:, save_cols])
    print(scale_values)
    print(scaled_data)

if __name__ == "__main__":
    main()
