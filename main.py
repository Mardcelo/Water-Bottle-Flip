# main.py

# Import modules
import numpy as np
import matplotlib.pyplot as plt
import kinetics
from constants import *

# Set kinetic parameters
water_mass = 0.1  # Mass of water in the bottle, unit: kg
filling_fraction = water_mass / WATER_MASS_MAX  # Filling fraction of the bottle
water_height = filling_fraction * BOTTLE_HEIGHT  # Height of the water level from the bottom of the bottle at rest, unit: m
epsilon = BOTTLE_MASS / (BOTTLE_MASS + water_mass)  # correction factor
alpha = 0.01  # Empirical drag coefficient, unit: kg/m
beta = -0.3  # Empirical restitution coefficient, unit: 1
omega_0 = 20  # Initial angular velocity, unit: rad/s
theta_0 = np.pi / 2  # Initial angle, unit: rad

# Set numerical parameters for the simulation
n_steps = 1000  # Number of time steps
t_max = 0.6  # Total duration of the flipping motion, unit: s
t, dt = np.linspace(0, t_max, n_steps, retstep=True)  # Time vector and time step
n_slices = 100  # Number of water slices
slice_height = water_height / n_slices  # Height of a slice, unit: m
slice_mass = water_mass / n_slices  # Mass of a slice, unit: kg

# Initialize the variables
omega = omega_0  # Initial angular velocity, unit: rad/s
theta = theta_0
slices_positions = BOTTLE_HEIGHT - slice_height / 2 - np.arange(0, n_slices) * slice_height
slices_velocities = 0 * slices_positions  # Initial radial velocity of the slices
center_of_mass_system = kinetics.find_center_of_mass(slices_positions, epsilon)
rotational_inertia_water = kinetics.rotational_inertia_water(slices_positions, center_of_mass_system, water_mass, n_slices)
rotational_inertia_bottle = kinetics.rotational_inertia_bottle(center_of_mass_system, BOTTLE_MASS, BOTTLE_RADIUS, BOTTLE_HEIGHT)
rotational_inertia = rotational_inertia_water + rotational_inertia_bottle

# Create vectors to store the results
omega_values = np.zeros(n_steps)
theta_values = np.zeros(n_steps)
slices_positions_values = np.zeros((n_steps, n_slices))
slices_velocities_values = np.zeros((n_steps, n_slices))
center_of_mass_system_values = np.zeros(n_steps)
rotational_inertia_values = np.zeros(n_steps)

# Store initial conditions
omega_values[0] = omega
theta_values[0] = theta
slices_positions_values[0] = slices_positions
slices_velocities_values[0] = slices_velocities
center_of_mass_system_values[0] = center_of_mass_system
rotational_inertia_values[0] = rotational_inertia

# Numerical integration
for k in range(1, len(t)):
    # Recall values of the variables
    omega = omega_values[k - 1]
    theta = theta_values[k - 1]
    slices_positions = slices_positions_values[k - 1]
    slices_velocities = slices_velocities_values[k - 1]
    center_of_mass_system = center_of_mass_system_values[k - 1]
    rotational_inertia = rotational_inertia_values[k - 1]

    # Update values of the variables using the previous values
    slices_positions_tmp, slices_velocities_tmp = kinetics.update_slice_positions(slices_positions, center_of_mass_system, omega, theta, slices_velocities, alpha, slice_mass, dt, G)
    new_slices_positions_tmp, new_slices_velocities_tmp = kinetics.check_boundary_conditions(slices_positions_tmp, slices_positions, slices_velocities_tmp, beta, slice_height / 2, BOTTLE_HEIGHT - slice_height / 2)
    center_of_mass_system_tmp = kinetics.find_center_of_mass(slices_positions, epsilon)
    rotational_inertia_water_tmp = kinetics.rotational_inertia_water(slices_positions, center_of_mass_system, water_mass, n_slices)
    rotational_inertia_bottle_tmp = kinetics.rotational_inertia_bottle(center_of_mass_system, BOTTLE_MASS, BOTTLE_RADIUS, BOTTLE_HEIGHT)
    rotational_inertia_tmp = rotational_inertia_water_tmp + rotational_inertia_bottle_tmp
    omega_tmp = omega_values[0] * rotational_inertia_values[0] / rotational_inertia_tmp
    theta_tmp = theta + omega * dt

    # Store new values of the variables
    omega_values[k] = omega_tmp
    theta_values[k] = theta_tmp
    slices_positions_values[k] = new_slices_positions_tmp
    slices_velocities_values[k] = new_slices_velocities_tmp
    center_of_mass_system_values[k] = center_of_mass_system_tmp
    rotational_inertia_values[k] = rotational_inertia_tmp

# Plot figures
fig, ax = plt.subplots(3, 1, figsize=(5, 10))
ax[0].plot(t, theta_values)
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Angle')
ax[1].plot(t, center_of_mass_system_values, '--')
ax[1].plot(t, slices_positions_values[:, 0], 'k')
ax[1].plot(t, slices_positions_values[:, -1], 'k')
ax[1].plot(t, slices_positions_values[:, n_slices // 2], 'k')
ax[1].set_xlabel('Time')
ax[1].set_ylabel('Center of Mass')
ax[1].set_ylim(0, 1.1 * BOTTLE_HEIGHT)
ax[2].plot(t, omega_values)
ax[2].set_xlabel('Time')
ax[2].set_ylabel('Angular Velocity')
ax[2].set_ylim(0, 1.1 * omega_0)
plt.tight_layout()
plt.show()

# Create data file (t, position center of mass, theta, omega)
# data = np.array([t, center_of_mass_system_values, theta_values, omega_values]).T
# np.savetxt('data.txt', temporary_data, delimiter=' ')
