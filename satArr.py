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
    ('names', nb.from_dtype(np.dtype('U10'))[:]),   # Name
    ('a', nb.float64[:]),                           # Semi-major axis
    ('ma', nb.float64[:]),                          # Mean anomaly
    ('omega', nb.float64[:]),                       # Angular velocity
    ('Omega', nb.float64[:]),
    ('Omega_dot', nb.float64[:]),                   # Precession rate
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
        self.pos, self.vel = KT.getXYZVVV(self.ma, self.a, self.omega, self.e, self.Omega, self.I, self.n)

    def updateOrbit(self, dt: float) -> None:
        self.ma = (self.ma + self.n*dt) % tau               # Advance anomaly
        self.Omega = (self.Omega + self.Omega_dot*dt) % tau # Advance node
        self.t += dt                                        # Advance time

        self.updateKinematics()

    def getConjunctionData(self, arr: np.ndarray) -> tuple:
        self.NConj = len(arr)

        # Create arrays of the indices of the involved satellites
        i1 = arr[:, 0]
        i2 = arr[:, 1]

        # Get satellite names
        satOne = self.names[i1]
        satTwo = self.names[i2]

        # Compute relative distances and velocities
        relD_vec = self.pos[i1] - self.pos[i2]
        relV_vec = self.vel[i1] - self.vel[i2]

        relD = KT.getNorm(relD_vec)
        relV = KT.getNorm(relV_vec)

        # Compute altitudes
        alt = np.sqrt(self.pos[0]**2 + self.pos[1]**2 + self.pos[2]**2) - REarth
        time = [self.t for i in range(len(arr))]

        return conjEvent(time, relD, relV, alt, i1, i2, satOne, satTwo)

    