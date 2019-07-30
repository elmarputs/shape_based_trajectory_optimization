import numpy as np
import scipy.integrate as sp_int


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

    def get_time_of_flight(self, gravitational_parameter, number_of_points):
        # Calculate integrand
        theta, r = self.get_radial_distance_from_theta(number_of_points)
        s = np.sin(self.k2 * theta + self.phi)
        c = np.cos(self.k2 * theta + self.phi)
        tan_gamma = self.k1*self.k2*c
        integrand = np.sqrt(np.power(r, 3)*(np.power(tan_gamma, 2) + self.k1*np.power(self.k2, 2)*s + 1)
                            / gravitational_parameter)

        time_of_flight = sp_int.simps(integrand, theta)
        return time_of_flight

