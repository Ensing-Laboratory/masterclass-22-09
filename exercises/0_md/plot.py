import numpy as np
import matplotlib.pyplot as plt

# Read the Muller-Brown data file
filename_pes = 'potential.dat'
with open(filename_pes, 'r') as f:
    lines = f.readlines()

# Parse the header to get the number of grid points
header = lines[0].strip()
nx, ny = map(int, header.split()[1:3])

# Initialize arrays for x, y, and z
x_pes = np.zeros(nx * ny)
y_pes = np.zeros(nx * ny)
z_pes = np.zeros(nx * ny)

# Read the data
data_lines = [line for line in lines[1:] if line.strip()]

for i, line in enumerate(data_lines):
    parts = line.split()
    x_pes[i] = float(parts[0])
    y_pes[i] = float(parts[1])
    z_pes[i] = float(parts[2])

# Reshape the data into a 2D grid
x_grid = x_pes.reshape((nx, ny))
y_grid = y_pes.reshape((nx, ny))
z_grid = z_pes.reshape((nx, ny))

# Read the trajectory file
filename_traj = 'colvar.out'
traj_data = np.loadtxt(filename_traj, skiprows=1)
x_traj = traj_data[:, 1]
y_traj = traj_data[:, 2]

# Create the contour plot
plt.figure(figsize=(8, 6))
contour = plt.contourf(x_grid, y_grid, z_grid, levels=50, cmap='viridis')
plt.colorbar(contour, label='Potential Energy')
plt.title('Muller-Brown Potential Energy Surface with Trajectory')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)

# Plot the trajectory as a line with points
plt.plot(x_traj, y_traj, color='red', linewidth=2, label='Trajectory')
plt.scatter(x_traj, y_traj, color='red', s=20, label='Trajectory Points')

plt.legend()
plt.show()

