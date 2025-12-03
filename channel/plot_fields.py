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
# ts_vec = [0]
# ts_vec = range(800000,840500,10000)
# ts_vec = range(0,900001,100000)
ts_vec = [60000]  # Test with just one timestep

# ts_vec = [10000]

# set 1 to compute time averaged quantities, 0 otherwise
timeaverage = 0

# set 1 to compute fluctuating components, 0 otherwise (expensive)
fluct = 0

# value for the fontsize:
fontsize_val = 10

# show heatmaps
showmaps_flag = 1

# slice of slice (leave -1 to compute the mean)
meanprof_slice = 0

# value of the figure size
figsize_val = 8

# save figure in png format
savefig_flag = 0

# save data in h5 format
savedata_flag = 0

# Aspect ratio of the heatmaps (-1 for auto, i.e. proportional, WARNING:can't handle huge AR)
AR_plot = -1

# Set your desired color limits here
vmin = None  # e.g., vmin = 0.0
vmax = None  # e.g., vmax = 1.0

# vmin = -1  # e.g., vmin = 0.0
# vmax = 1  # e.g., vmax = 1.0

#################################################################
# AUXILIARY FUNCTIONS
def get_spectrum(x, signal, title='Spectrum Analysis'):

    # Compute the Fourier Transform
    spectrum = np.fft.fft(signal)
    
    # Calculate the amplitude and phase
    amplitude = np.abs(spectrum)
    phase = np.angle(spectrum)

    # Frequency values corresponding to the FFT components
    frequencies = np.fft.fftfreq(len(x), d=(x[1] - x[0]))

    # cut negative freq
    positive_freq_indices = frequencies >= 0
    frequencies = frequencies[positive_freq_indices]
    amplitude = amplitude[positive_freq_indices]
    phase = phase[positive_freq_indices]

    # Sort frequencies and corresponding arrays
    sorted_indices = np.argsort(frequencies)
    frequencies = frequencies[sorted_indices]
    amplitude = amplitude[sorted_indices]
    phase = phase[sorted_indices]

    # crop 2/3 of freq
    amplitude = amplitude[1:int(np.floor(len(frequencies)*2/3))]
    phase = phase[1:int(np.floor(len(frequencies)*2/3))]
    frequencies = frequencies[1:int(np.floor(len(frequencies)*2/3))]

    wavenum = frequencies*2*np.pi
    return [wavenum,amplitude,phase]

#################################################################

# font = {'family' : 'normal',
#         'weight' : 'normal',
#         'size'   : 6}

plt.rcParams.update({'font.size': 5})

# params_txt = {'axes.labelsize': 5,'axes.titlesize':5, 'text.fontsize': 5, 'legend.fontsize': 5, 'xtick.labelsize': 5, 'ytick.labelsize': 5}
# plt.rcParams.update(params_txt)

# solve conflicts: if fluctuations is compute mean is needed
if fluct == 1:
    timeaverage = 1

if timeaverage == 1:
    meansuff = '_MEAN'
else:
    meansuff = ''

done = 0
# GENERATED ON:
from datetime import datetime

# Print current date and time briefly
print("Current date and time:", datetime.now())

glob_profiles = []

FLUCARRAYS_time = []

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

for n_step in ts_vec:
    file_names = []

    for fld in fields:
        file_names.append(fld + '_{:08d}.dat')

    # Initialize an empty list to store the data arrays
    data_arrays = []

    # initialize arrays with sums if want to calculate time averages
    if timeaverage == 1:
        flucarray = []
        if n_step == ts_vec[0]:
            sumarray = []

    # Read the data from each file and reshape
    id_fnames = -1
    for fld in fields:
        file_name = f"{fld}_{n_step:08d}.dat"
        file_path = foldername + file_name
        id_fnames = id_fnames+1
        
        # Check if file exists
        if not os.path.exists(file_path):
            print(f"Warning: File {file_path} not found, skipping...")
            continue
            
        with open(file_path, 'rb') as file:
            total_elements = np.prod(nvec)
            data = np.memmap(file, dtype=np.float64, mode='r', shape=(total_elements,))
            data = data.reshape(np.flip(nvec))*1.0

            data = data.transpose((2, 1, 0)) # Permute first and third index to match the convection [x,z,y]
            
            # Validate data
            print(f"Loaded {file_name}: shape={data.shape}, min={data.min():.6f}, max={data.max():.6f}, mean={data.mean():.6f}")

            dersuff = ''
            if derivative_vec[0] != 0:
                dersuff = '_'
                for ider in derivative_vec:
                    derivative_x, derivative_z, derivative_y = np.gradient(data, x, z, y)
                    if ider == 1:
                        data = derivative_x
                        dersuff = dersuff+'x'
                    elif ider == 2:
                        data = derivative_y
                        dersuff = dersuff+'y'
                    elif ider == 3:
                        data = derivative_z
                        dersuff = dersuff+'z'

            # define axes:
            if slice_dir == 1:
                hor_name = 'y'
                ver_name = 'z'
                hor = y
                nhor = ny
                ver = z
                nver = nz
            elif slice_dir == 2:
                hor_name = 'x'
                ver_name = 'z'
                hor = x
                nhor = nx
                ver = z
                nver = nz
            elif slice_dir == 3:
                hor_name = 'x'
                ver_name = 'y'
                hor = x
                nhor = nx
                ver = y
                nver = ny

            if slice_idx == -1:
                type_name = 'average'
                if slice_dir == 1:
                    mean_array = np.mean(data, axis=0)
                elif slice_dir == 2:
                    mean_array = np.mean(data, axis=1)
                elif slice_dir == 3:
                    mean_array = np.mean(data, axis=2)

            else:
                type_name = 'slice'
                if slice_dir == 1:
                    mean_array = data[slice_idx-1,:,:]
                elif slice_dir == 2:
                    mean_array = data[:,slice_idx-1,:]
                elif slice_dir == 3:
                    mean_array = data[:,:,slice_idx-1]

            data_arrays.append(mean_array)

            if timeaverage == 1:
                if fluct == 1:
                    flucarray.append(mean_array)
                if n_step == ts_vec[0]:
                    sumarray.append(mean_array)
                else:
                    sumarray[id_fnames] = sumarray[id_fnames]+mean_array

                if n_step == ts_vec[-1]:
                    data_arrays[id_fnames] = sumarray[id_fnames]/len(ts_vec)

                    FLUCARRAYS_time.append(flucarray)

    # Plot each array as a heat map with coordinates

    if timeaverage == 0 or n_step == ts_vec[-1]:
        if showmaps_flag == 1:
            for i, array in enumerate(data_arrays):
                if array.ndim != 2:
                    print(f"Invalid shape for {file_name}")
                    continue
                
                N = 500j
                extent = (hor.min(),hor.max(),ver.min(),ver.max())

                HOR,VER = np.meshgrid(hor,ver)
                hors,vers = np.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]

                plt.figure(figsize=(figsize_val, figsize_val*0.12+figsize_val))
                
                ax = plt.gca()
                # im = ax.imshow(np.flip(resampled.T, axis = 0), cmap='jet', origin='lower', extent=extent, aspect='auto')#extent=[y.min(), y.max(), z.min(), z.max()])
                # im = ax.imshow(resampled.T, cmap='jet', origin='lower', extent=extent, aspect='auto')#extent=[y.min(), y.max(), z.min(), z.max()])
                # Use pcolormesh for non-uniform grids
                X, Y = np.meshgrid(hor,ver)

                im = ax.pcolormesh(X, Y, array.T, cmap='jet', shading='gouraud', vmin=vmin, vmax=vmax)  # smoother shading

                ax.set_aspect('equal')  # <-- Add this line to set axis equal

                plt.xlabel(hor_name, fontsize=fontsize_val)
                plt.ylabel(ver_name, fontsize=fontsize_val)
                plt.xticks(fontsize=fontsize_val, rotation=0)
                plt.yticks(fontsize=fontsize_val, rotation=0)

                # create an axes on the right side of ax. The width of cax will be 5%
                # of ax and the padding between cax and ax will be fixed at 0.05 inch.
                cax = plt.axes([0.15, 0.95, 0.7, 0.03])
                plt.colorbar(im, cax=cax,orientation = "horizontal")
                cax.xaxis.set_ticks_position('top')
                plt.xticks(fontsize=fontsize_val, rotation=0)
                plt.title(file_name+dersuff,fontsize=fontsize_val)

                plt.show()

                # save arrays and figure in the case of temporal average
                if timeaverage == 1:
                    if savefig_flag == 1:
                        plt.savefig('./'+fields[i]+dersuff+'_'+str(ts_vec[0])+'_'+str(ts_vec[-1])+'.png', dpi=800)
                    
                    # Save the arrays to an HDF5 file using h5py
                    if savedata_flag == 1:
                        with h5py.File('./'+fields[i]+dersuff+'_'+str(ts_vec[0])+'_'+str(ts_vec[-1])+'timeav.h5', 'w') as hf:
                            hf.create_dataset('data', data=array)
                            hf.create_dataset('x', data=hor)
                            hf.create_dataset('y', data=ver)


    # # # # # PLOT ALL LINES IN THE SAME FIGURE
    # # Plot the vertical profiles in the opposite direction with different colors
    # plt.figure(figsize=(8, 6))

    # for i, array in enumerate(data_arrays):
    #     if array.ndim != 2:
    #         print(f"Invalid shape for {file_name}")
    #         continue

    #     # Compute the vertical profile in the opposite direction
    #     vertical_profile = np.mean(array, axis=0)[::-1]

    #     # Plot the vertical profile with a different color for each array
    #     plt.plot(vertical_profile, ver, label=file_name, alpha=0.7)

    # plt.title(fields[i]+'Profiles',fontsize=fontsize_val)
    # plt.xlabel(fields[i] +'Mean',fontsize = fontsize_val)
    # plt.ylabel(ver_name,fontsize = fontsize_val)
    # plt.legend(loc ="best",fontsize = fontsize_val)
    # plt.xticks(fontsize=fontsize_val, rotation=0)
    # plt.yticks(fontsize=fontsize_val, rotation=0)
    # plt.show()

    # Plot the vertical profiles in the opposite direction with different colors

    vertical_profile = []

    for i, array in enumerate(data_arrays):
        if array.ndim != 2:
            print(f"Invalid shape for {file_name}")
            continue

        # Compute the vertical profile in the opposite direction
        if meanprof_slice == -1:
            vertical_profile.append(np.mean(array, axis=0))
            horproflabelaux = ' Mean'
        else:
            vertical_profile.append(np.transpose(array[meanprof_slice-1,:]))
            horproflabelaux = ' '
            

        # plt.figure(figsize=(8, 6))

        # # Plot the vertical profile with a different color for each array
        # plt.plot(vertical_profile[i], ver, label=file_name, alpha=0.7)

        # plt.title(fields[i]+' Profiles',fontsize=fontsize_val)
        # plt.xlabel(fields[i] +' Mean',fontsize = fontsize_val)
        # plt.ylabel(ver_name,fontsize = fontsize_val)
        # plt.legend(loc ="best",fontsize = fontsize_val)
        # plt.xticks(fontsize=fontsize_val, rotation=0)
        # plt.yticks(fontsize=fontsize_val, rotation=0)
        # plt.show()

    glob_profiles.append(vertical_profile)

prf_plt = []

if fluct == 1:
    for i, array in enumerate(data_arrays):
        mean = data_arrays[id_fnames]
        rmsfluc = np.zeros_like(mean)
        for j in range(len(ts_vec)):
            flucs = FLUCARRAYS_time[j]
            singfield = flucs[i]
            rmsfluc = rmsfluc+np.square(singfield-mean)
        rmsfluc = np.sqrt(rmsfluc/len(ts_vec))

        # PLOT RMS FIELD
        N = 500j
        extent = (hor.min(),hor.max(),ver.min(),ver.max())

        HOR,VER = np.meshgrid(hor, ver)
        hors,vers = np.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
        aaa = np.ones((len(HOR.flatten()),2))
        aaa[:,0] = HOR.flatten()
        aaa[:,1] = VER.flatten()

        resampled = griddata( aaa, rmsfluc.T.flatten(), (hors, vers), method='linear')

        if done == 0:
            done = 1
            if AR_plot == -1:
                AR_plot = 1/(hor.max()-hor.min())*(ver.max()-ver.min())
                AR_plt_flag = 0
            else:
                AR_plt_flag = 1

        plt.figure(figsize=(figsize_val, figsize_val*0.12+figsize_val*AR_plot))
        
        ax = plt.gca()
        if AR_plt_flag == 1:
            # im = ax.imshow(np.flip(resampled.T, axis = 0), cmap='jet', origin='lower', extent=extent, aspect='auto')#extent=[y.min(), y.max(), z.min(), z.max()])
            im = ax.imshow(resampled.T, cmap='jet', origin='lower', extent=extent, aspect='auto')#extent=[y.min(), y.max(), z.min(), z.max()])
        else:
            # im = ax.imshow(np.flip(resampled.T, axis = 0), cmap='jet', origin='lower', extent=extent, aspect=1)#extent=[y.min(), y.max(), z.min(), z.max()])
            im = ax.imshow(resampled.T, cmap='jet', origin='lower', extent=extent, aspect=1)#extent=[y.min(), y.max(), z.min(), z.max()])

        # plt.imshow(np.flip(resampled.T, axis = 0), cmap='jet', origin='lower', extent=extent)#extent=[y.min(), y.max(), z.min(), z.max()])
        # plt.colorbar(orientation = "horizontal")

        # plt.title(file_name,fontsize=fontsize_val)
        plt.xlabel(hor_name,fontsize=fontsize_val)
        plt.ylabel(ver_name,fontsize=fontsize_val)
        plt.xticks(fontsize=fontsize_val, rotation=0)
        plt.yticks(fontsize=fontsize_val, rotation=0)

        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax,orientation = "horizontal")
        cax.xaxis.set_ticks_position('top')
        plt.xticks(fontsize=fontsize_val, rotation=0)
        plt.title(file_name+dersuff+'_RMS',fontsize=fontsize_val)

        plt.show()

        # SAVE FLUCTUATION FIELDS
        if savedata_flag == 1:
            with h5py.File('./'+fields[i]+dersuff+'_'+str(ts_vec[0])+'_'+str(ts_vec[-1])+'time_rms.h5', 'w') as hf:
                hf.create_dataset('data', data=rmsfluc)
                hf.create_dataset('x', data=hor)
                hf.create_dataset('y', data=ver)

numfields = len(glob_profiles[0])

for i in range(numfields):

    plt.figure(figsize=(figsize_val, 9))

    for j in range(len(ts_vec)):

        profile_sets = glob_profiles[j]
        single_profile = profile_sets[i]
                
        # plt.figure(figsize=(8, 6))

        # Plot the vertical profile with a different color for each array
        plt.plot(single_profile, ver, label=fields[i]+'_'+str(ts_vec[j]), alpha=0.7)


    if fld == 'u' and slice_dir != 3:
        zmin = np.min(ver)
        zmax = np.max(ver)
        zc = 0.5 * (zmin + zmax)
        zlen = zmax - zmin

        # Parabolic profile: u(z) = umax * (1 - ((z - zc)/H)^2), with H = zlen/2
        H = zlen / 2.0
        umax = 1.0 / (2.0 * 1.0)  # assuming mu = rho = 1 → umax = 1/(2μ)
        zvals = np.linspace(zmin, zmax, 256)
        u_analytical = umax * (1.0 - ((zvals - zc)/H)**2)

        plt.plot(u_analytical, zvals, 'r--', linewidth=2, label='Poiseuille (analytical)')

    plt.title(fields[i]+dersuff+' Profiles',fontsize=fontsize_val)
    plt.xlabel(fields[i]+dersuff + horproflabelaux,fontsize = fontsize_val)
    plt.ylabel(ver_name,fontsize = fontsize_val)
    # plt.legend(loc ="best",fontsize = fontsize_val)
    plt.legend(bbox_to_anchor=(1.00, 1.02), loc="upper left",fontsize = fontsize_val)
    plt.subplots_adjust(right=0.70)
    plt.xticks(fontsize=fontsize_val, rotation=0)
    plt.yticks(fontsize=fontsize_val, rotation=0)
    plt.show()

for i in range(numfields):

    plt.figure(figsize=(figsize_val, 8))

    for j in range(len(ts_vec)):

        profile_sets = glob_profiles[j]
        single_profile = profile_sets[i]
                
        # plt.figure(figsize=(8, 6))

        # Plot the vertical profile with a different color for each array
        [freq,amp,phase] = get_spectrum(ver,single_profile)

        # Plot amplitude on log scale
        plt.subplot(2, 1, 1)
        plt.plot(freq, amp,label=fields[i]+'_'+str(ts_vec[j]), alpha=0.7)
        plt.xscale('log')
        if np.any(amp > 0):  # Only use log scale if there are positive values
            plt.yscale('log')
        plt.title('Amplitude')
        plt.xlabel('k_'+ver_name)
        plt.ylabel('Amplitude')

        # Plot phase on log scale
        plt.subplot(2, 1, 2)
        plt.plot(freq, phase,label=fields[i]+'_'+str(ts_vec[j]), alpha=0.7)
        plt.xscale('log')
        plt.title('Phase')
        plt.xlabel('k_'+ver_name)
        plt.ylabel('Phase (radians)')

    # Adjust layout and show the plot
    plt.subplot(2, 1, 1)
    plt.legend(bbox_to_anchor=(1.00, 1.02), loc="upper left",fontsize = fontsize_val)
    plt.tight_layout()
    plt.show()
