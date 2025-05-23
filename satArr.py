import numpy as np
import numba as nb

from numba.experimental import jitclass
from collections import namedtuple

import KeplerTools as KT

# Constants
tau = 2*np.pi
REarth = 6378.135e3 

# Data Structures
satSpec = (
    ('names', nb.from_dtype(np.dtype('U30'))[:]),   # Name
    ('a', nb.float64[:]),                           # Semi-major axis
    ('ma', nb.float64[:]),                          # Mean anomaly
    ('omega', nb.float64[:]),                       # Angular velocity
    ('Omega', nb.float64[:]),                       # Node angle
    ('Omega_dot', nb.float64[:]),                   # Nodal precession rate
    ('n', nb.float64[:]),                           # Inverse period
    ('e', nb.float64[:]),                           # Eccentricity
    ('I', nb.float64[:]),                           # Inclination
    ('pos', nb.float64[:, :]),                      # Position
    ('vel', nb.float64[:, :]),                      # Velocity
    ('t', nb.float64),                              # Current time
    ('NSat', nb.int64),                           # Total number of objects
    ('NConj', nb.int64)                           # Total number of conjunctions at a given time
)

conjEvent = namedtuple('conjunctionEvent', ['t', 'dist', 'vel', 'alt', 'id1', 'id2', 'name1', 'name2'])

@jitclass(satSpec)
class satArray(object):
    """
    Class to hold an array of satellite objects with methods to propagate them along a Keplerian orbit and extract conjunction data.

    Attributes:
        names (str array): Names of all satellite objects
        a (array): Semi-major axes
        ma (array): Mean anomalies
        omega (array): Angular velocities
        Omega (array): Node angles
        Omega_dot (array): Nodal precession rates
        n (array): Rate of sweep (mean angular motion)
        e (array): Eccentricities
        I (array): Orbital inclination angles
        pos (array): 3D Cartesian position coordinates
        vel (array): 3D Cartesian velocity components
        t (float): Current time (calculated as the sum of all timesteps propagated)
        NSat (int): Total number of satellite objects
        NConj (int): The total number of conjunctions below a specified distance threshold
    """

    def __init__(self, names: np.ndarray, a: np.ndarray, ma: np.ndarray, omega: np.ndarray, 
                 Omega: np.ndarray, Omega_dot: np.ndarray, n: np.ndarray, e: np.ndarray, I: np.ndarray, pos: np.ndarray, vel: np.ndarray):
        
        self.names = names
        self.a = a
        self.ma = ma
        self.omega = omega
        self.Omega = Omega
        self.Omega_dot = Omega_dot
        self.n = n
        self.e = e
        self.I = I
        self.pos = pos
        self.vel = vel

        self.t = 0.
        self.NSat = int(len(self.names))
        self.NConj = 0

    # Simulation methods
    def updateKinematics(self) -> None:
        """
        Recalculates all satellite positions and velocities using the currently stored orbital elements.
        """
        self.pos, self.vel = KT.getXYZVVV(self.ma, self.a, self.omega, self.e, self.Omega, self.I, self.n)

    def updateOrbit(self, dt: float) -> None:
        """
        Advances the orbit of all satellites by a specified time-step dt, including positions and velocities.
        """
        self.ma = (self.ma + self.n*dt) % tau               # Advance anomaly
        self.Omega = (self.Omega + self.Omega_dot*dt) % tau # Advance node
        self.t += dt                                        # Advance time

        self.updateKinematics()

    def getConjunctionData(self, arr: np.ndarray) -> tuple:
        """
        Returns conjunction data for selected satellites at the current time-step.

        Parameters:
            arr (array): An array containing arrays with pairs of indices for which conjunction data should be calculated.

        Returns:
            conjEvent (namedTuple): A tuple of arrays containing data for all conjunctions accessible via keyword:
                Time: ['t']
                Relative distance: ['dist']
                Relative velocity: ['vel']
                Altitude of first satellite: ['alt']
                Index of first satellite: ['id1']
                Index of second satellite: ['id2']
                Name of first satellite: ['name1']
                Name of second satellite: ['name2']
        """
        self.NConj = len(arr)

        # Create arrays of the indices of the involved satellites
        i1 = arr[:, 0]
        i2 = arr[:, 1]

        # Get satellite names
        satOne = self.names[i1]
        satTwo = self.names[i2]

        # Compute relative distances and velocities
        relD = KT.getNorm(self.pos[i1] - self.pos[i2])
        relV = KT.getNorm(self.vel[i1] - self.vel[i2])

        # Compute altitudes
        alt = np.sqrt(self.pos[0]**2 + self.pos[1]**2 + self.pos[2]**2) - REarth
        time = [self.t for i in range(len(arr))]

        return conjEvent(time, relD, relV, alt, i1, i2, satOne, satTwo)
    
    def randomizeOrbits(self) -> None:
        """
        Distributes mean anomalies and nodal angles randomly according to a uniform distribution. Simulates behaviour on long timescales.
        """
        for i in range(self.NSat):
            self.ma[i]    = np.random.uniform(0, tau - 1e-8) # Subtract a little bit to avoid overlap
            self.Omega[i] = np.random.uniform(0, tau - 1e-8)

    