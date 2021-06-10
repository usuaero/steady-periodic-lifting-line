import json

import numpy as np
import matplotlib.pyplot as plt


def vec_cross(x, y):
    """Calculates the cross product of the last dimensions of x and y."""
    xT = x.T
    yT = y.T
    return np.array([xT[1]*yT[2]-xT[2]*yT[1],
                     xT[2]*yT[0]-xT[0]*yT[2],
                     xT[0]*yT[1]-xT[1]*yT[0]]).T


def vec_norm(x):
    """Calculates the norm of the last dimension of x."""
    xT = x.T
    return np.sqrt(xT[0]*xT[0]+xT[1]*xT[1]+xT[2]*xT[2]).T


def vec_inner(x, y):
    """Calculates the inner product of the last dimensions of x and y."""
    xT = x.T
    yT = y.T
    return (xT[0]*yT[0]+xT[1]*yT[1]+xT[2]*yT[2]).T


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

    max_iter : int
        Maximum iterations to compute.

    planar : bool
        Whether to make the wake planar.

    V_init : ndarray
        Initial velocity for setting the wake direction.

    V_inf : float
        Freestream velocity.
    """

    def __init__(self, **kwargs):


        # Get kwargs
        self.b = kwargs.get("span")
        self.CLa = kwargs.get("section").get("CLa")
        self.aL0 = kwargs.get("section").get("aL0")
        self.N = kwargs.get("N", 20)
        self.planar = kwargs.get("planar")

        # Initialize grid
        self._setup_grid()

        # Initialize chord
        self._setup_chord(kwargs.get("chord"))

        # Initialize wake
        self._setup_wake(kwargs.get("max_iter"))

        self.k = 1


    def _setup_chord(self, c):
        # Initializes the chord distribution

        # Check for elliptic
        if isinstance(c, list):
            self.c = c[1]*np.sqrt(1.0-4.0*self.cp_span_locs**2)
            self.c_root = c[1]
        
        # Constant
        else:
            self.c = np.ones_like(self.cp_span_locs)*c
            self.c_root = c

        # Calculate dS
        self.dS = self.c*vec_norm(self.dl)


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
        
        # Calculate dl
        self.dl = np.zeros((self.N, 3))
        self.dl[:,1] = self.y_node[1:]-self.y_node[:-1]


    def _setup_wake(self, max_iter):
        # Initializes the wake points

        # Initialize wake storage
        self.wake_points = np.zeros((max_iter+1, self.N+1, 3))

        # Store trailing points for initial ring
        if self.planar:
            self.wake_points[0,:,0] = -100.0*self.c_root
            self.wake_points[0,:,1] = self.y_node
        else:
            raise IOError("Non-planar wakes are not allowed.")


    def get_influences(self, cp, k):
        """Determines the influences of the vortex rings on the points cp for the k-th wake iteration assuming unit vortex strength.

        Parameters
        ----------
        cp : ndarray
            Points at which to calculate the influences.

        k : int
            Iteration index.

        Returns
        -------
        v_mji : ndarray
            Velocity induced on the i-th control point (first axis) by a ring from the m-th iteration (second axis), j-th spanwise station (third axis). Fourth axis is the vector component.
        """

        # Get radial distances
        points = self.wake_points[:k+1,:,:]
        r = cp[:,np.newaxis,np.newaxis,:]-points[np.newaxis,:,:,:]
        r_mag = vec_norm(r)

        # Set up for transverse calculation
        r_1 = r[:,:,:-1,:]
        r_2 = r[:,:,1:,:]
        r_1_mag = r_mag[:,:,:-1]
        r_2_mag = r_mag[:,:,1:]

        # Get transverse vortex influences
        v_tran = vec_cross(r_1, r_2)*((r_1_mag+r_2_mag)/(4.0*np.pi*r_1_mag*r_2_mag*(r_1_mag*r_2_mag+vec_inner(r_1, r_2))))[:,:,:,np.newaxis]

        # Set up for longitudinal calculation
        r_1 = r[:,:k,:,:]
        r_2 = r[:,1:k+1,:,:]
        r_1_mag = r_mag[:,:k,:]
        r_2_mag = r_mag[:,1:k+1,:]

        # Get longitudinal vortex influences
        v_long = vec_cross(r_1, r_2)*((r_1_mag+r_2_mag)/(4.0*np.pi*r_1_mag*r_2_mag*(r_1_mag*r_2_mag+vec_inner(r_1, r_2))))[:,:,:,np.newaxis]

        # Initialize influence matrix
        v_mji = np.zeros((self.N, k, self.N, 3))

        # Add downstream vortex
        v_mji -= v_tran[:,:k,:,:]

        # Add right vortex
        v_mji -= v_long[:,:,1:,:]

        # Add upstream vortex
        v_mji[:,:-1,:,:] += v_tran[:,1:-1,:,:]

        # Add left vortex
        v_mji += v_long[:,:,:-1,:]

        return v_mji