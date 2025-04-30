import numpy as np
import numba as nb

from numba.experimental import jitclass

import KeplerTools as KT

# Constants
tau = 2*np.pi

satSpec = (
    ('names', nb.from_dtype(np.dtype('U10'))[:]),
    ('a', nb.float64[:]),
    ('ma', nb.float64[:]),
    ('omega', nb.float64[:]),
    ('Omega', nb.float64[:]),
    ('Omega_dot', nb.float64[:]),
    ('n', nb.float64[:]),
    ('e', nb.float64[:]),
    ('I', nb.float64[:]),
    ('pos', nb.float64[:, :]),
    ('vel', nb.float64[:, :]),
    ('t', nb.float64),
    ('NSat', nb.float64),
    ('NConj', nb.float64)

)

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
        self.NSat = len(self.names)
        self.NConj = 0

    # Simulation methods
    def updateKinematics(self) -> None:
        self.pos, self.vel = KT.getXYZVVV(self.ma, self.a, self.omega, self.e, self.Omega, self.I, self.n)

    def updateOrbit(self, dt: float) -> None:
        self.ma = (self.ma + self.n*dt) % tau               # Advance anomaly
        self.Omega = (self.Omega + self.Omega_dot*dt) % tau # Advance node
        self.t += dt

    def getConjunctionData(self, arr: np.ndarray) -> tuple:
        self.NConj = len(arr)
        pairs = np.array([(arr[i, 0], arr[i, 1]) for i in range(self.NConj)])

        # Pickup here
        relD = self.pos[pairs[:, 0]] - self.pos[pairs[:, 1]]
        relV = self.vel[pairs[:, 0]] - self.vel[pairs[:, 1]]

    