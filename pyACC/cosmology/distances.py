import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import trapezoid, simpson, quad
from scipy.constants import c

class Distances:
    """
    Class to calculate cosmological distances given the Hubble function
    """

    def __init__(self, hubble_function):
        """
        Parameters
        ----------
        hubble_function : callable
            Hubble function of the form H(z) = H0 * E(z)
        """
        self.hubble_function = hubble_function

    def comoving_distance(self, z, n_nodes):
        """
        Comoving distance in Mpc
        Parameters
        ----------
        z : float
            Redshift
        n_nodes : int
            Nodes number for the redshift integration array    
        Returns
        -------
        float
            Comoving distance in Mpc    
        """
        zz = np.linspace(0, z, n_nodes) # integration interval array
        Integrand = 1. / self.hubble_function(zz) # integrand function
        
        Integral_trapz = (c / 1000.) * trapezoid(Integrand, zz) # trapezoid integral
        
        return Integral_trapz