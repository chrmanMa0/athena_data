from scipy.interpolate import CubicSpline
import numpy as np
from matplotlib import pyplot as plt
matrix = np.loadtxt("s25.0-CS1.2-H1C0.1-1784660003.txt")
print(matrix)



s = CubicSpline(matrix[:,1], matrix[:,2], bc_type='natural')
splined_values = s(matrix[:,1])
print(splined_values)
plt.plot(matrix[:, 1], splined_values)
plt.scatter(matrix[:, 1], matrix[:, 2], color='red', s=5)
plt.xscale("log")
plt.yscale("log")
plt.show()
