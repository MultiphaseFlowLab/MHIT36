# %%
import numpy as np
import matplotlib.pyplot as plt

nx, ny, nz = 256, 128, 200

data = np.fromfile('output/u_00005000.dat', dtype=np.float64)  # or float64
print("Data size:", data.size)
data = data.reshape((nz, ny, nx))  

slice_index = ny // 2  # middle of z
slice_data = data[ :,slice_index, :]  # xy slice at fixed z

print("Max:", np.max(data))
print("Min:", np.min(data))
plt.figure(figsize=(14,6))
plt.imshow(slice_data, cmap='jet', origin='lower', aspect='3.14')
plt.colorbar(label='Velocity')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('scaled')
plt.show()
#1D plot
#line_data = data[1, 1, : ]  # middle column
#plt.subplot(1, 2, 2)
#plt.plot(line_data, np.arange(ny), color='black')
#plt.xlabel('Variable')
#plt.ylabel('y')
#plt.title('1D Profile along y (x = nx/2)')
#plt.grid(True)

#plt.tight_layout()
#plt.show()

# %%
