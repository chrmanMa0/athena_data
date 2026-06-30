import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
filename = "s15.0_ready.txt"
data_frame = pd.read_csv(filename, sep=r'\s+')
np_data = data_frame.iloc[:, 2]
np_data_rho = data_frame.iloc[:, 4]
helium_data = data_frame[data_frame.index <= 584]
helium_data_reduced = helium_data.iloc[:, [2, 4, 6]]
helium_data_reduced.to_csv('helium_data.txt', sep=' ', index=False)

