import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# USER INPUT SECTION
# -----------------------------
Lx = float(input("Enter length in X direction (m): "))
Ly = float(input("Enter length in Y direction (m): "))
Nx = int(input("Enter number of grid points in X: "))
Ny = int(input("Enter number of grid points in Y: "))
alpha = float(input("Enter thermal diffusivity (m^2/s): "))
dt = float(input("Enter time step (s): "))
Nt = int(input("Enter number of time steps: "))
Tb = float(input("Enter boundary temperature (°C): "))
Th = float(input("Enter initial hot spot temperature (°C): "))

# -----------------------------
# GRID SETUP
# -----------------------------
dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)

x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)

# -----------------------------
# INITIAL CONDITION
# -----------------------------
T = np.zeros((Nx, Ny))

# Hot spot at center
T[Nx//2, Ny//2] = Th

# -----------------------------
# TIME INTEGRATION (FTCS METHOD)
# -----------------------------
for n in range(Nt):
    T_new = T.copy()
    
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
            T_new[i, j] = T[i, j] + alpha * dt * (
                (T[i+1, j] - 2*T[i, j] + T[i-1, j]) / dx**2 +
                (T[i, j+1] - 2*T[i, j] + T[i, j-1]) / dy**2
            )
    
    # Boundary conditions
    T_new[0, :] = Tb
    T_new[-1, :] = Tb
    T_new[:, 0] = Tb
    T_new[:, -1] = Tb

    T = T_new

# -----------------------------
# VISUALIZATION
# -----------------------------
plt.figure(figsize=(6,5))
plt.contourf(x, y, T.T, 50)
plt.colorbar(label="Temperature (°C)")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.title("2D Heat Transfer PDE (User Input Based)")
plt.show()
