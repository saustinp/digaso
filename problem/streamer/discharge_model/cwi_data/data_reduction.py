import numpy as np
from scipy.io import savemat

# Doing them in sequence to avoid dealing with the different numbers of zeros
for i in np.arange(9)+1:
    snapshot = np.loadtxt(f'comp_1e13_f_line_00000{i}.txt', skiprows=1)
    mdic = {'snapshot':snapshot}
    savemat(f'snapshot{i}.mat', mdic)

for i in np.arange(9,23)+1:
    snapshot = np.loadtxt(f'comp_1e13_f_line_0000{i}.txt', skiprows=1)
    mdic = {'snapshot':snapshot}
    savemat(f'snapshot{i}.mat', mdic)
