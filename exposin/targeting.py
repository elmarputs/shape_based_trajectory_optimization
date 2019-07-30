import numpy as np
from shape_based.exposin.exposin import Exposin

class ExposinTargeting:

    def __init__(self, windingParameter, radialDistanceDep, radialDistanceArr, psiAngle, numberOfRevs):
        self.k2 = windingParameter
        self.r1 = radialDistanceDep
        self.r2 = radialDistanceArr
        self.psi = psiAngle
        self.N = numberOfRevs
        self.thetaBar = psiAngle + 2*np.pi*numberOfRevs

    def find_flight_path_angle_bounds(self):

        delta = (2*(1 - np.cos(self.k2 * self.thetaBar)))/np.power(self.k2, 4) - np.power(np.log(self.r1 / self.r2), 2)
        common_term = (-1 * np.log(self.r1 / self.r2)) / np.tan(self.k2 * self.thetaBar / 2)
        lower_limit = np.arctan(0.5 * self.k2 * (common_term - np.sqrt(delta)))
        upper_limit = np.arctan(0.5 * self.k2 * (common_term + np.sqrt(delta)))

        return lower_limit, upper_limit

    def find_exposin_parameters(self, initial_flight_path_angle):
        g1 = initial_flight_path_angle
        common_term = np.log(self.r1 / self.r2) + np.tan(g1) / self.k2 * np.sin(self.k2 * self.thetaBar)

        k1 = np.sqrt(np.power(common_term / (1 - np.cos(self.k2 * self.thetaBar)), 2)
                     + np.power(np.tan(g1) / self.k2, 2))
        sign = np.sign(common_term)
        k1 = k1 * sign

        phi = np.arccos(np.tan(g1) / (k1 * self.k2))

        k0 = self.r1 / (np.exp(k1 * np.sin(phi)))

        return k0, k1, phi

    def find_exposin(self, initial_flight_path_angle):

        k0, k1, phi = self.find_exposin_parameters(initial_flight_path_angle)

        return Exposin(k0, k1, self.k2, phi, self.N, self.psi)

