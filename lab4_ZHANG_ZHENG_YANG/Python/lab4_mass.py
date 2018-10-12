""" Force-Velocity Setup """
import numpy as np
import biolog


def mass_equation(pos, vel, force, mass_params):
    """ Mass equation"""
    d2x = -force/mass_params.mass + mass_params.g
    return d2x


def mass_system(pos, vel, force, mass_params):
    """ Muscle-Mass System"""
    return np.array(
        [vel,
         mass_equation(pos, vel, force, mass_params)])  # xdd)

