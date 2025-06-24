## Written by Aaron Boley. June 2021
## Super klunky. 
## Awesome licence stuff 
## No guarantees

## Made even cooler by Skye Heiland. February 2025.

# External dependencies
import numpy as np
import pandas as pd
import sgp4.api as sgp4
import csv
import datetime

from scipy.spatial import KDTree

# Internal
import KeplerTools as KT
from satArr import satArray

# Select TLE Catalogue
# Also set epoch information for TLEs
TLES_FILENAME = "./in/full_cat_linked.tles"
# TLES_FILENAME = "./in/full_cat_linked.tles"
OUT_SUFFIX = 'dev_test'
CHECKPOINT_FILENAME = f'./out/checkpoint_{OUT_SUFFIX}.dat'
JD = 2460591.5
FR = 0.0
JDOFFSET = 100
EPOCHLIM = 24280

# Constants and parameters
ECCSCALE = 0.0002         # Make eccentric enough to fill shells for any artificial systems
SMASCALE = 1000           # km to metres
DT = 0.05                 # Time step in seconds
TTIME = 30           # How long to run sim (s)
NTIME = int(TTIME / DT)   # Number of steps

EXAMINE_PHS = False                   # Flag to determine if we examine phase-space mixing and produce histograms
PHS_INT_S = 3600                      # How many seconds to wait between calculating phase-space coords,
PHS_INT = int(PHS_INT_S / DT)        # and how many time steps.
PHS_FRAMES = int(NTIME / PHS_INT)    # How many frames of our phase space plot we'll have.

DCLOSE_METRE = 10000.        # track if closer than this

PLOT = False           # If we want to create the plots directly in satEvol
RANDOM_ORBITS = False   # Randomizes nodes and mean anomalies before running
LINK = True          # Links objects that are close together before running

twopi = np.pi*2             # How many pi??
MEarth = 5.97e24            # Mass of Earth (kg)
REarth = 6378.135e3         # Radius of Earth (m)
REkm = REarth/SMASCALE      # Radius of Earth (km)

P_THRESH = 2000e3 + REarth  # do not include objects with pericentres above this

aukm = 1.496e8            # Astronomical unit (km)
au = aukm*SMASCALE        # Astronomical unit (m)
G = 6.6743e-11            # Big G (SI)
J2 = 1082.64e-6           # Arcane kinematic parameters
muE = 3.986004418e14

# Helper I/O functions
# Read infile of TLEs and return a structured array of satellites
def readTLEs(TLE_path: str, PHASE_path: str, check_linked: bool = True) -> tuple:
    """
    Reads a list of two or three-line element sets from a text file.

    Parameters:
        TLE_path (str): Path to text file containing TLEs
        PHASE_path (str): Path for phase outfile containing names and initial kinematic variables
        check_linked (bool, optional): Whether or not to check for objects that should be linked

    Returns:
        tuple: Tuple of arrays containing satellite properties where the ith element corresponds to the ith satellite:
            Name
            Semi-major axis         [m]
            Mean anomaly            [rad]
            Angular velocity        [rad/s]
            Nodal precession angle  [rad]
            Nodel precession rate   [rad/s]
            Inverse period
            Orbital eccentricity
            Orbital inclination     [deg]
            Initial position        [m]
    """
    woh = open(PHASE_path,"w")

    with open(TLE_path,'r') as f:
        lines = f.readlines()
        TotalSat = int(len(lines)/3) # Total number of satellites, assuming three lines per

        names       = np.zeros(TotalSat, dtype='U30')
        a           = np.zeros(TotalSat, dtype='f8')
        ma          = np.zeros(TotalSat, dtype='f8')
        omega       = np.zeros(TotalSat, dtype='f8')
        omega_dot   = np.zeros(TotalSat, dtype='f8')
        Omega       = np.zeros(TotalSat, dtype='f8')
        Omega_dot   = np.zeros(TotalSat, dtype='f8')
        n           = np.zeros(TotalSat, dtype='f8')
        ecc         = np.zeros(TotalSat, dtype='f8')
        inc         = np.zeros(TotalSat, dtype='f8')
        pos         = np.zeros(TotalSat, dtype=('f8', 3))
        vel         = np.zeros(TotalSat, dtype=('f8', 3))

        satIndex = 0
        for i, line in enumerate(lines):
            if len(line) < 1: break
            if line[0] == '0':
                sname = line.rstrip().partition(' ')[2]

            elif line[0] == '1':
                s = line.rstrip()

            elif line[0] == '2':
                t = line.rstrip()
                sat = sgp4.Satrec.twoline2rv(s, t)

                if sat.jdsatepoch < JD - JDOFFSET:
                    print(f"ISSUE WITH JDs: WANTING {JD} GOT {sat.jdsatepoch}")
                    continue

                err, r_km, v_km = sat.sgp4(JD, FR)
                sat_vel = np.array([v_km[0],v_km[1],v_km[2]])*SMASCALE
                sat_pos = np.array([r_km[0],r_km[1],r_km[2]])*SMASCALE
                sat_a, sat_ecc, sat_omega, sat_inc, sat_Omega, sat_nu = KT.getORBELM(sat_pos, sat_vel, muE)

                sat_ma = KT.meanAnom(sat_ecc, sat_nu)  # Mean anomaly
                sat_n = KT.meanAngularMotion(sat_a) # Mean angular motion
                sat_omega_dot, sat_Omega_dot = KT.getPrecessionRate(sat_a, sat_ecc, sat_n, sat_inc)

                if sat_a*(1 - sat_ecc) < P_THRESH:

                    names[satIndex]     = sname
                    a[satIndex]         = sat_a
                    ma[satIndex]        = sat_ma
                    omega[satIndex]     = sat_omega
                    omega_dot[satIndex] = sat_omega_dot
                    Omega[satIndex]     = sat_Omega
                    Omega_dot[satIndex] = sat_Omega_dot
                    n[satIndex]         = sat_n
                    ecc[satIndex]       = sat_ecc
                    inc[satIndex]       = sat_inc
                    pos[satIndex]       = sat_pos
                    vel[satIndex]       = sat_vel

                    woh.write("{},{},{},{},{},{},{}\n".format(sname, r_km[0], r_km[1], r_km[2], v_km[0], v_km[1], v_km[2]))
                    satIndex += 1

    woh.close()

    if check_linked:
        print('Checking for linked objects...')

        woh_df = pd.read_csv(PHASE_path, header=None, encoding='us-ascii')
        woh_df.columns = ['name','x','y','z','vx','vy','vz']

        woh_df = woh_df[['x', 'y', 'z']].round(2)   # This is a jank solution I need to clean up

        woh_keep = woh_df.drop_duplicates(['x','y','z'], keep='first').index.to_numpy()

        names       = names[woh_keep]
        a           = a[woh_keep]
        ma          = ma[woh_keep]
        omega       = omega[woh_keep]
        omega_dot   = omega_dot[woh_keep]
        Omega       = Omega[woh_keep]
        Omega_dot   = Omega_dot[woh_keep]
        n           = n[woh_keep]
        ecc         = ecc[woh_keep]
        inc         = inc[woh_keep]
        pos         = pos[woh_keep]
        vel         = vel[woh_keep]

        print(f'No. of satellites linked: {TotalSat - len(woh_keep)}')

    return names, a, ma, omega, omega_dot, Omega, Omega_dot, n, ecc, inc, pos, vel

def writeOutfile(conjList: tuple, fname: str, type='close-approach') -> None:
    """
    Writes a list of conjunction events to an outfile.
    """
    with open(fname, 'a', newline='') as f:
        for i in range(len(conjList)):
            conj = conjList[i]._asdict()
            writer = csv.writer(f, delimiter=',')
            writer.writerows(zip(conj['t'], conj['dist'], conj['vel'], conj['alt'],
                                 conj['id1'], conj['id2'], conj['name1'], conj['name2']))
            
def createCheckpoint(sat: satArray, fname: str) -> None:
    """
    Writes the orbital elements of all satellites at the current time-step to a checkpoint file.
    """
    with open(fname, 'w', newline='') as f:
        f.write(f'## Checkpoint created for run {OUT_SUFFIX} at t = {sat.t} on {datetime.datetime.now()} \n')
        f.write('## Name, semi-major axis, mean anom., arg. of peri., apsidal precession, RAAN, nodal precession, mean ang. motion, ecc., inc. \n')

        writer = csv.writer(f, delimiter=',')
        writer.writerows(zip(sat.names, sat.a, sat.ma, sat.omega, sat.omega_dot, sat.Omega, 
                             sat.Omega_dot, sat.n, sat.e, sat.I))

        
# Distance searching function (can't be jitted because of SciPy)
def threshQuery(pos: np.ndarray, d: float, pairwise: bool = True):
    """
    Constructs a k-dimensional tree out of satellite positions to search for nearest neighbours.

    Parameters:
        pos (array): Array of 3D satellite positions.
        d (float): The distance threshold to count to objects as close.
        pairwise (bool): Whether or not to calculate this distance pairwise, or for all satellites.

    Returns:
        array: An array containing arrays with pairs of indices corresponding to satellites below the distance threshold.
    """
    posTree = KDTree(pos, leafsize=15)
    if pairwise: return posTree.query_pairs(d, output_type='ndarray')
    else: return posTree.query_ball_tree(posTree, d)

# ================= MAIN BLOCK ================= #

def main() -> None:

    # Nerd stuff
    F_CHUNK = int(1e1)        # How many time steps we wait between dumping to outfile
    NEWFILE_CHUNK = int(1e5)  # How many time steps we wait before starting a new output file
    NEWFILE_INDEX = 1         # Output file indexing

    PHASE_FOUT = f"./out/catdata_xyzvxvyvz_{OUT_SUFFIX}.dat"
    OUTFILE = f"./out/track_out_{OUT_SUFFIX}.dat"

    # Read in TLEs
    print('Reading TLEs...')
    satData = readTLEs(TLES_FILENAME, PHASE_FOUT, check_linked=LINK)

    Satellites = satArray(satData[0], satData[1], satData[2], satData[3], satData[4], satData[5],
                          satData[6], satData[7], satData[8], satData[9], satData[10], satData[11])  # For some reason this doesn't construct if I just pass the tuple in
    
    print(f'Total no. of satellites: {Satellites.NSat}')

    # Initialize satellite positions
    if RANDOM_ORBITS: 
        print('Randomizing orbits...')
        Satellites.randomizeOrbits()

    Satellites.updateKinematics()

    simTime = 0
    conjunctionData = []
    
    # Main integration
    print(f'Beginning main integration for {TTIME}s')
    for itime in range(NTIME):
        Satellites.updateOrbit(DT)

        closeIndices = threshQuery(Satellites.pos, DCLOSE_METRE, pairwise=True)
        conjunctionData.append(Satellites.getConjunctionData(closeIndices))

        # Dump to outfile
        if (itime + 1) % F_CHUNK == 0: 
            writeOutfile(conjunctionData, OUTFILE)
            conjunctionData.clear()

        # Tick up the outfile index
        if (itime + 1) % NEWFILE_CHUNK == 0:
            OUTFILE = f"./out/track_out_{OUT_SUFFIX}_{NEWFILE_INDEX}.dat"
            NEWFILE_INDEX += 1

        if (itime*DT) % 3600 == 0:
            createCheckpoint(Satellites, CHECKPOINT_FILENAME)

        print(f'Time: {simTime:.3f}s  |  No. of conjunctions: {Satellites.NConj}', end='\r')
        simTime += DT

    writeOutfile(conjunctionData, OUTFILE)
    print('Done!')


if __name__ == '__main__':
    main()
    quit()
        
