#!/usr/bin/python3

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import argparse


a_cgs = 7.56577e-15
sigma_cgs = 5.6704e-5
G_cgs = 6.6743e-8

class WR_Star_With_Wind:


    def __init__(self, stellar_data, scales):
        self.L_0 = None
        self.M_0 = None
        self.rho_0 = None
        self.t_0 = None
        self.v_0 = None
        self.G_0 = None
        self.P_0 = None
        self.T_0 = None
        self.E_0 = None



        self.stellar_data = None
        self.v_0 = None
        self.v_esc = None
        self.v_inf = None
        self.density_distribution = None
        self.velocity_distribution = None

        self.stellar_data = stellar_data
        print(self.stellar_data.columns)
        self.stellar_data.columns = ['cell outer radius', 'cell enclosed mass',
                                     'cell density', 'cell pressure', 'cell energy density',
                                     'cell temperature', 'helium fraction', 'hydrogen fraction']
        self.L_0 = scales[0]
        self.M_0 = scales[1]
        self.rho_0 = scales[2]
        self.t_0 = scales[3]
        self.v_0 = scales[4]
        self.G_0 = scales[5]
        self.P_0 = scales[6]
        self.T_0 = scales[7]
        self.E_0 = scales[8]
        self.a_0 = scales[9]
        self.sigma_0 = scales[10]
        self.kappa_e0 = scales[11]

    """
    This function will generate the internal wind structure from
    the previously ingested data frame.  It should not be called
    before the data frame has been constructed.
    """
    def generate_wind(self):
        self.extract_basic_parameters()
        if(self.v_0 == None or self.v_esc == None or self.v_inf == None):
            print("Wind generation failed.  No standard parameters set")
            return
        else:
            #Time to construct the distribution data from the basic parameters
            print("Generating the distributions\n")
    @classmethod
    def from_file(cls, filename, scale_filename):
        #we need to set up and return a cls
        #with stellar data
        stellar_df = pd.read_csv(filename, sep=" ", header=None)
        scales = np.loadtxt(scale_filename)
        return cls(stellar_df, scales)


    """
    This function will extract the basic parameters from a 2sn processed file
    that is intended to be processed with Athena++.
    calculate v_0 - estimate as the sound speed (P/rho)**.5
    calculate v_esc - v_esc =((2GM_star)/R_star)**.5
    calculate v_inf

    """
    def extract_basic_parameters(self):
        #First we will need to extract the last cell info
        #for pressure and density temperature, radius, etc.  This will determine the
        #base gas pressure at the start of the wind and all are in program units
        P = self.stellar_data['cell pressure'].iloc[-1]
        rho = self.stellar_data['cell density'].iloc[-1]
        M_star = self.stellar_data['cell enclosed mass'].iloc[-1]
        T_star = self.stellar_data['cell temperature'].iloc[-1]
        R_star = self.stellar_data['cell outer radius'].iloc[-1]
        #calculate the gas pressure velocity at the surface --Not this is NOT
        #the v_0 tha is a conversion factor.
        self.v_pres = (P/rho)**.5
        #uncorrected escape calculation
        #Using outer cell of core as the radius and enclosed mass
        #remember there is actually a photosphere further in the
        #wind which is part of the complexity
        self.G_prog = G_cgs/self.G_0
        self.v_esc = ((2.0*self.G_prog*M_star)/R_star)**.5
        #Need to correct for proper electron scaling
        #v_esc,eff = v_esc * (1 - gamma_e)**.5
        #gamma_e = (sigma_e * L_*)/(4*pi*G*M)
        #sigma_e = 0.4*(1+z_ion * Y)/(1 + 4Y)
        #Y is helium mass fraction z is how many free electrons (2 for these temps)

        self.sigma_prog = sigma_cgs/self.sigma_0
        print(f"Sigma_prog: {self.sigma_prog}")
        print(f"Temp: {T_star}")
        self.L_star = 4 * np.pi * R_star**2*self.sigma_prog*T_star**4
        print(f"Luminosity: {self.L_star}")
        self.X = self.stellar_data['hydrogen fraction'].iloc[-1]
        print(f"Hydrogen fraction: {self.X}")
        print(f"Helium fraction: {self.stellar_data['helium fraction'].iloc[-1]}")
        self.Y = self.stellar_data['helium fraction'].iloc[-1]#(self.stellar_data['helium fraction'].iloc[-1]/4)/self.X
        print(f"Y: {self.Y}")
        self.z_ion = 2.0 #mostly he4 in this region so 2 electrons

        kappa_e = 0.4 * ((1.0-self.Y-.02) + .5*self.Y + .5*(.02*.5)) #(0.4*(1.0+self.z_ion*self.Y)/(1.0+4.0*self.Y))#/self.kappa_e0
        print(f"kappa_e: {kappa_e}")
        kappa_e_prog = kappa_e/self.kappa_e0
        print(f"kappa_e_prog: {kappa_e_prog}")
        self.G_prog = G_cgs/self.G_0
        c_prog = 2.99792e10/self.v_0
        self.gamma_e = (kappa_e_prog * self.L_star)/(4*c_prog*np.pi*self.G_prog*M_star)
        print(f"Gamma: {self.gamma_e}")
        #program scale luminosity of the outer shell of the core

        self.v_esc_eff = self.v_esc * (1.0-self.gamma_e)**.5

        #calculate v_inf using empirical scaling factor Howarth & Prinja 1989
        self.v_inf = 2.65 * self.v_esc_eff


    def construct_wind_parameters(self):
        self.v_0 = 0
        self.v_esc = 0
        self.v_inf = 0

    def __str__(self):
        return f"radius: {self.L_0/6.957e10}, v_pres: {self.v_pres}, v_esc_eff: {self.v_esc_eff}, v_inf: {self.v_inf}"

def main():
    parser = argparse.ArgumentParser(description="Generates a wolf rayet compatible wind for a given core in program units")
    parser.add_argument("data_filename", type=str, help="Name of the datafile containing the wolfrayet core")
    parser.add_argument("scales_filename", type=str, help="Name of the scales file used to convert the data to program data")
    parser.add_argument("-p", "--plot", type=bool, default=False)
    args = parser.parse_args()
    filename = args.data_filename
    filename_scales = args.scales_filename

    wind = WR_Star_With_Wind.from_file(filename, filename_scales)
    wind.generate_wind()
    print(wind)

if __name__ == "__main__":
    main()

WR_mass = 7.1441 #in solar masses

def nl_mass_loss_from_mass(mass_sm):
    return -5.99 + 1.06*np.log10(mass_sm)


print(nl_mass_loss_from_mass(WR_mass))

