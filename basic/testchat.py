import numpy as np
import matplotlib.pyplot as plt

# Define the initial conditions
r = 6371.0  # Earth radius in kilometers
mu = 1.0  # Initial refractive index
rho_r = 0.3  # Initial derivative of refractive index with respect to r
rho_theta = 0.1  # Initial derivative of refractive index with respect to theta
rho_phi = 0.3  # Initial derivative of refractive index with respect to phi
dr_dt = 0.1  # Initial derivative of r with respect to time
dtheta_dt = 0.0  # Initial derivative of theta with respect to time
dphi_dt = 0.2  # Initial derivative of phi with respect to time

# Define the time step and maximum number of steps
dt = 0.01  # Time step size
max_steps = 200000  # Maximum number of steps

# Define arrays to store the trajectory
r_values = [r]
theta_values = [0.0]
phi_values = [0.0]

# Perform the ray tracing
for step in range(max_steps):
    # Calculate the derivatives
    denominator = r * mu**2 * np.sin(theta_values[-1])
    
    if denominator != 0:
        dphi_dt = (1 / denominator) * (rho_phi - mu * np.cos(theta_values[-1]) * np.sin(theta_values[-1]) * rho_r -
                                       mu * r * np.cos(theta_values[-1]) * np.sin(theta_values[-1]) * rho_theta)
    else:
        dphi_dt = 0.0
    
    # Update the variables
    dr_dt = (1 / mu**2) * (rho_r - mu * np.cos(theta_values[-1]) * np.sin(theta_values[-1]) * rho_theta -
                           mu * np.sin(theta_values[-1]) * np.sin(theta_values[-1]) * rho_phi)
    dtheta_dt = (1 / (r * mu**2)) * (rho_theta - mu * np.cos(theta_values[-1]) * np.sin(theta_values[-1]) * rho_r -
                                     mu * r * np.sin(theta_values[-1]) * np.sin(theta_values[-1]) * rho_phi)
    
    # Update the variables
    r += dr_dt * dt
    theta = theta_values[-1] + dtheta_dt * dt
    phi = phi_values[-1] + dphi_dt * dt
    
    # Append the new values to the arrays
    r_values.append(r)
    theta_values.append(theta)
    phi_values.append(phi)
    
    # Check if the ray has reached the end or escaped to infinity
    if r <= 0 or r >= 2 * r:
        break

# Convert the spherical coordinates to Cartesian coordinates
x_values = r_values * np.sin(theta_values) * np.cos(phi_values)
y_values = r_values * np.sin(theta_values) * np.sin(phi_values)
z_values = r_values * np.cos(theta_values)

# Plot the ray trace
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_values, y_values, z_values)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Ray Tracing in Earth\'s Magnetic Field')
plt.show()