# kinetics.py

# Import modules
import numpy as np
from constants import *

# Compute the rotational inertia of the water with respect to the center of mass
def rotational_inertia_water(positions, center_of_mass, water_mass, n_slices):
    # j = np.sum(mass * (positions - center_of_mass) ** 2)
    j = water_mass / n_slices * np.sum((positions - center_of_mass) ** 2)
    return j

# Compute the rotational inertia of the bottle with respect to the center of mass
def rotational_inertia_bottle(center_of_mass, mass_bottle=BOTTLE_MASS, radius=BOTTLE_RADIUS, length=BOTTLE_HEIGHT):
    j_bottle_axis = mass_bottle * (radius ** 2 / 2 + length ** 2 / 12)
    j_new_axis = j_bottle_axis + mass_bottle * (length / 2 - center_of_mass) ** 2  # Parallel axis theorem
    return j_new_axis

# Compute the position of the center of mass
def find_center_of_mass(positions, eps, length_bottle=BOTTLE_HEIGHT):
    r_center_of_mass = eps * length_bottle / 2 + (1 - eps) * np.mean(positions)
    return r_center_of_mass

# Update slice positions and velocities using previous values
def update_slice_positions(positions, center_of_mass, angular_velocity, angle, velocities, alpha, mass, dt, gravity):
    tau = (2 * mass) / alpha
    factor_1 = positions - center_of_mass + (gravity / (angular_velocity ** 2)) * np.cos(angle)
    factor_2 = velocities / angular_velocity + (1 / (tau * angular_velocity)) * (positions - center_of_mass + (gravity / (angular_velocity ** 2)) * np.cos(angle))
    new_positions = center_of_mass - (gravity / (angular_velocity ** 2)) * np.cos(angle) + factor_1 * np.exp(- dt / tau) * np.cosh(angular_velocity * dt) + factor_2 * np.exp(- dt / tau) * np.sinh(angular_velocity * dt)
    factor_1 = ((angular_velocity ** 2 - (1 / (tau ** 2))) / angular_velocity) * (positions - center_of_mass + (gravity / (angular_velocity ** 2)) * np.cos(angle)) + (1 / (tau * angular_velocity)) * velocities
    factor_2 = velocities
    new_velocities = factor_1 * np.exp(- dt / tau) * np.sinh(angular_velocity * dt) + factor_2 * np.exp(- dt / tau) * np.cosh(angular_velocity * dt)
    return new_positions, new_velocities

# Take into account the rebound on the boundaries
def check_boundary_conditions(positions, previous_positions, velocities, restitution, l_min, l_max):
    for i in range(len(positions)):
        if positions[i] < l_min:
            velocities[i] = restitution * velocities[i]
            positions[i] = previous_positions[i].copy()  # Just use the previous position
            # positions[i] = l_min + (l_min - positions[i])  # Would make more sense to have a symmetrical rebound
        elif positions[i] > l_max:
            velocities[i] = restitution * velocities[i]
            positions[i] = previous_positions[i].copy()  # Just use the previous position
            # positions[i] = l_max - (positions[i] - l_max)  # Would make more sense to have a symmetrical rebound
    return positions, velocities
