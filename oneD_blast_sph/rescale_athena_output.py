#! /usr/bin/python3
import os
os.environ["NUMEXPR_MAX_THREADS"] = "3232"

import pandas as pd
import numpy as np
import argparse
pd.set_option("compute.use_numexpr", False)





# Useful constants

M_dot = 1.989e33
R_dot = 6.957e10
a = 7.5657e-15
h1_early = .35
h1_late = .1
G = 6.674e-8


def store_scaled_data(out_fn, data_df, time):
    cgs_out_fn = out_fn.rsplit(".", 1)[0]+"-cgs.dat"
    with open(cgs_out_fn, "w") as f:
        f.write("# Time: " + str(time) + "\n")
    data_df.to_csv(cgs_out_fn, sep=" ", index=False, mode="a")

def read_scale_values(filename):
    scale_values = np.loadtxt(filename)
    row_names = ['L_0', 'M_0', 'rho_0', 't_0', 'v_0', 'G_0', 'P_0', 'T_0', 'E_0', "a_0", "sigma_0", "kappa_0"]
    scale_values_df = pd.DataFrame(data=scale_values, index=row_names)
    return scale_values_df

def load_program_scale_data(filename):
    with open(filename, 'r') as f:
        time_info = f.readline().strip()
        time_info_field = time_info.split(" ")
        program_time_val = time_info_field[4].split('=')
        program_time_val = float(program_time_val[1])
        column_info = f.readline().strip()
        column_names = column_info.split()
        data = np.loadtxt(f)
        return data, program_time_val

def rescale_data(unscaled_data, scale_values, unscaled_time):
    unscaled_data[:,1] = unscaled_data[:,1] * scale_values.loc['L_0', 0]
    unscaled_data[:,2] = unscaled_data[:,2] * scale_values.loc['rho_0', 0]
    unscaled_data[:,3] = unscaled_data[:,3] * scale_values.loc['P_0', 0]
    unscaled_data[:,4] = unscaled_data[:,4] * scale_values.loc['v_0', 0]
    unscaled_data[:,5] = unscaled_data[:,5] * scale_values.loc['v_0', 0]
    unscaled_data[:,6] = unscaled_data[:,6] * scale_values.loc['v_0', 0]
    columns = ['cell id', 'xv1', 'rho', 'press', 'vel1', 'vel2', 'vel3']
    scaled_time = unscaled_time * scale_values.loc['t_0', 0]
    scaled_data_df = pd.DataFrame(data = unscaled_data, columns=columns)
    return scaled_data_df, scaled_time

def main():
    parser = argparse.ArgumentParser(description="Returns athena++ program scaled data to cgs units")
    parser.add_argument("filename", type=str, help="Name of the file containing scale values for a program scaled dataset e.g. XXXXXXX-scales.txt")
    parser.add_argument("data_filename", type=str, help="Athena++ output file in program scale")

    args = parser.parse_args()

    filename = args.filename
    out_fn = args.data_filename


    scales = read_scale_values(filename)
    unscaled_data, unscaled_time = load_program_scale_data(out_fn)
    scaled_data, scaled_time = rescale_data(unscaled_data, scales, unscaled_time)
    store_scaled_data(out_fn, scaled_data, scaled_time)


if __name__ == "__main__":
    main()
