import pandas as pd
import numpy as np

# 1. Load your downloaded 2sn.org CSV data
#df = pd.read_csv('s20.0_stripped_and_cored.txt', skiprows=1, sep=r'\s+{2,}', engine='python')
df = pd.read_fwf('s25.0', skiprows=1)
print(df.keys())

# 2. Enforce MESA structural logic
# (MESA typically expects profile data ordered from surface to center)
df = df.sort_values(by='cell outer total mass', ascending=False).reset_index(drop=True)

# 3. Create target columns and apply unit transformations if needed
# Example: Convert radius to log10 if required, or ensure density is log10
df['logRho'] = np.log10(df['cell density'])
df['logT'] = np.log10(df['cell temperature'])
df['grid'] = range(len(df))
df['lum'] = 4.0*np.pi*df['cell outer radius']**2 * df['cell temperature']**4
df['grid'] = df['grid']+1
# Add your custom isotope columns here (must match your MESA reaction network)
#mesa_columns = ['zone', 'mass', 'radius', 'logRho', 'logT', 'lumin', 'velocity', 'h1', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24']
df_mesa = df[['grid', 'cell outer radius', 'cell outer total mass', 'lum', 'cell pressure', 'cell temperature', 'cell density']]
nabla_num = np.diff(df['logT'])
nabla_den = np.diff(df['cell pressure'])
nabla = nabla_num/nabla_den
nabla_ad = np.zeros(len(df_mesa)-1)+0.4
gamma = np.zeros(len(df_mesa)-1)+1.33
n2 = np.zeros(len(df_mesa)-1)
kap_kap_t = np.zeros(len(df_mesa)-1)
kap_kap_rho = np.zeros(len(df_mesa)-1)
zeros = np.zeros(len(df_mesa)-1)
kap = np.zeros(len(df_mesa)-1)+0.34
mesa_np = np.column_stack((df_mesa['grid'].iloc[:-1]-693, df_mesa['cell outer radius'].iloc[:-1], df_mesa['cell outer total mass'].iloc[:-1], df_mesa['lum'].iloc[:-1], df_mesa['cell pressure'].iloc[:-1], df_mesa['cell temperature'].iloc[:-1], df_mesa['cell density'].iloc[:-1], nabla, n2, gamma, nabla_ad, zeros, kap, kap_kap_t, kap_kap_rho, zeros, zeros, zeros, zeros, zeros))
print(df_mesa)
#df_mesa.columns = mesa_columns
# Select and rename to official MESA structural variables
# 4. Write out the standard MESA plain-text structure

output_filename = '2sn_start_model.mod'
num_zones = len(df_mesa)
stellar_mass = df_mesa["cell outer total mass"].iloc[0]
r_star = df_mesa["cell outer radius"].iloc[0]
r_lum = df_mesa["lum"].iloc[0]
with open(output_filename, 'w') as f:
    # MESA header metadata block
    f.write(str(num_zones-1) + " " + str(stellar_mass) + " " + str(r_star) + " " + str(r_lum) + " " + str(120) + "\n")  # Total grid zones
    #f.write(" ".join(mesa_columns) + "\n")  # Variable names

    # Save the space-separated numerical grid arrays
    fmt_list = ['%d'] + ['%f'] * 19
    np.savetxt(f, mesa_np[693:], fmt=fmt_list, delimiter=' ')

print(f"Successfully generated {output_filename} for MESA integration.")
