# # %%
# %%html
# <style>
# .cell-output-ipywidget-background {
#    background-color: transparent !important;
# }
# .jp-OutputArea-output {
#    background-color: transparent;
# }  
# .dataframe th {
#     font-size: 6px;
# }
# .dataframe td {
#     font-size: 6px;
# }
# </style>

# %%
# %matplotlib widget
# %matplotlib inline

import ipympl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as style
# from ing_theme_matplotlib import mpl_style
# from qbstyles import mpl_style

import os
import glob

from scipy.interpolate import griddata


from matplotlib.widgets import Cursor

import h5py

plt.style.use("dark_background")
# mpl_style(dark=True)


#################################################################
foldername = './output/'

# PARAMETERS:

# fields to plot
fields = ['theta']

# number of points in each direction
# Grid parameters (user-defined)
nx = 256  # number of points in x
ny = 128  # number of points in y
nz = 200  # number of points in z

nx = 2  # number of points in x
ny = 2  # number of points in y
nz = 1000  # number of points in z

Lx = 6.0  # length of domain in x
Ly = 3.0      # length of domain in y
Lz = 2.0      # length of domain in z

# compute the derivative of the fields (show them instead of the neormal fields)
# 0: no derivative
# 1: x derivative
# 2: y derivative
# 3: z derivative
# list more flag to compute consecutive derivatives (forder 1 FD)
derivative_vec = [0]

# # normal direction of the 2D slice:
# 1: x-direction
# 2: y-direction
# 3: z-direction
slice_dir = 2

# index to take the slice (from 1 to nx_i, choose -1 for computing the average)
slice_idx = 0

# slice_idx = 222

# time_steps to plot 
ts_vec = range(0,230000,10000)

# ts_vec = [10000]

# set 1 to compute time averaged quantities, 0 otherwise
timeaverage = 0

# set 1 to compute fluctuating components, 0 otherwise (expensive)
fluct = 0

# value for the fontsize:
fontsize_val = 10

x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
z = np.linspace(0, Lz, nz)

# Define the dimensions of the reshaped arrays (nvec) [y z x]

nvec = (nx, ny, nz)  # [y, z, x] order in your data file

# nvec = (512, 513, 512)  # Update with the actual dimensions
# nvec = (256, 257, 256)
# # nvec = (128, 129, 128)
# nvec = (0, 0, 0)

nx = nvec[0]
ny = nvec[1]
nz = nvec[2]

id_fnames = -1

for fld in fields:

    id_fnames = id_fnames+1


    plt.figure(figsize=(10, 9))

    for n_step in ts_vec:
        file_names = []

        file_names.append(fld + '_{:08d}.dat')


        # Read the data from each file and reshape

        file_name = f"{fld}_{n_step:08d}.dat"
        file_path = foldername + file_name
        
        # Check if file exists
        if not os.path.exists(file_path):
            print(f"Warning: File {file_path} not found, skipping...")
            continue
            
        with open(file_path, 'rb') as file:
            total_elements = np.prod(nvec)
            data = np.memmap(file, dtype=np.float64, mode='r', shape=(total_elements,))

            data = data.reshape(np.flip(nvec))*1.0
            
            # Validate data
            print(f"Loaded {file_name}: shape={data.shape}, min={data.min():.6f}, max={data.max():.6f}, mean={data.mean():.6f}")

            prof = np.mean(data, axis=(1, 2)) #.transpose((1,0,2))

            # prof = np.flip(np.mean(data, axis=(0, 2)))

            plt.plot(prof, z, label=f'{fld}_{n_step}', alpha=0.7)

    plt.title(fld+' Profiles',fontsize=fontsize_val)
    plt.xlabel(fld,fontsize = fontsize_val)
    plt.ylabel("z",fontsize = fontsize_val)
    # plt.legend(loc ="best",fontsize = fontsize_val)
    plt.legend(bbox_to_anchor=(1.00, 1.02), loc="upper left",fontsize = fontsize_val)
    plt.subplots_adjust(right=0.70)
    plt.xticks(fontsize=fontsize_val, rotation=0)
    plt.yticks(fontsize=fontsize_val, rotation=0)
    plt.show()
