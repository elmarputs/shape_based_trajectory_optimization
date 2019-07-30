import numpy as np


class Exposin:

    def __init__(self, scaling_factor, dynamic_range, winding_parameter, phase_angle, number_of_revs, psi_angle):
        self.k0 = scaling_factor
        self.k1 = dynamic_range
        self.k2 = winding_parameter
        self.phi = phase_angle
        self.N = number_of_revs
        self.psi = psi_angle

    def get_radial_distance_from_theta(self, number_of_points):

        theta = np.linspace(0.0, self.psi + 2*np.pi*self.N, number_of_points)
        r = self.k0 * np.exp(self.k1 * np.sin(self.k2 * theta + self.phi))

        return theta, r

