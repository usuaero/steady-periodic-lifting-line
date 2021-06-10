import json

import numpy as np
import matplotlib.pyplot as plt


class Wing:
    """Defines a wing for steady-periodic lifting-line analysis.

    Parameters
    ----------
    span : float
        Wingspan in feet.

    chord : float or list
        Chord length in feet. May be ["elliptic", ##] to specify an elliptic planform with root chord of ## feet.

    section : dict
        Airfoil section data. Should be

        {   
            "CLa" : <lift slope in rad^-1>,
            "aL0" : <zero-lift angle of attack in rad^-1>
        }

    N : int, optional
        Number of discrete sections to use for this wing. Defaults to 20.
    """

    def __init__(self, **kwargs):


        # Get kwargs
        self.b = kwargs.get("span")
        self.c = kwargs.get("chord")
        self.CLa = kwargs.get("section").get("CLa")
        self.aL0 = kwargs.get("section").get("aL0")
        self.N = kwargs.get("N", 20)

        # Initialize grid
        self._setup_grid()
        plt.figure()
        plt.plot(self.cp_span_locs, np.zeros(self.N), 'ro')
        plt.plot(self.node_span_locs, np.zeros(self.N+1), 'bo')
        plt.show()


    def _setup_grid(self):
        # Creates cosine-clustered grid

        # Initialize span location storage
        node_span_locs = [0.0]
        cp_span_locs = []

        # Create node distribution
        theta = list(np.linspace(0.0, np.pi, self.N+1))
        self.node_span_locs = -0.5*np.cos(theta)

        # Create control point distribution
        theta = np.linspace(np.pi/self.N, np.pi, self.N)-np.pi/(2*self.N)
        self.cp_span_locs = -0.5*np.cos(theta)

        # Initialize points
        self.y_node = self.b*self.node_span_locs
        self.y_cp = self.b*self.cp_span_locs