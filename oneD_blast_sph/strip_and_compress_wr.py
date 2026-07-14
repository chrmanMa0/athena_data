#! /usr/bin/python3
import os
os.environ["NUMEXPR_MAX_THREADS"] = "3232"

import pandas as pd
import numpy as np

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


#useful indices
h1_fraction_column = 15
cell_enclosed_mass = 1
cell_outer_radius = 2
cell_density = 4
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



def main():
    filename = "s25.0"
    full_data_np = load_and_prep_data(filename)
    pre_stripped_radius = full_data_np[-1, cell_outer_radius]
    stripped_data, remaining_mass, outer_radius = strip_shell(full_data_np, h1_early)
    print(remaining_mass/M_dot)
    print(outer_radius/R_dot)



if __name__ == "__main__":
    main()
