## Written by Aaron Boley. June 2021
## Super klunky. 
## Awesome licence stuff 
## No guarantees

## Made even cooler by Skye Heiland. February 2025.

# External dependencies
import numpy as np
import numba as nb
import matplotlib.pylab as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sgp4.api as sgp4
import csv

from scipy.spatial import KDTree
# Internal
import KeplerTools as KT
from satArr import satArray

# Select TLE Catalogue
# Also set epoch information for TLEs
# TLES_FILENAME = "./in/starlink_07OCT2024.tles"
TLES_FILENAME = "./in/full_cat.tles"
OUT_SUFFIX = '_test'
PHASE_FOUT = "./out/catdata_xyzvxvyvz" + OUT_SUFFIX +".dat"
OUTFILE = './out/track_out' + OUT_SUFFIX + '.dat'
JD = 2460591.5
FR = 0.0
JDOFFSET = 100
EPOCHLIM = 24280

# Constants and parameters
ECCSCALE = 0.0002         # Make eccentric enough to fill shells for any artificial systems
SMASCALE = 1000           # km to metres
DT = 0.05                 # Time step in seconds
TTIME = 60           # How long to run sim (s)
NTIME = int(TTIME / DT)   # Number of steps

EXAMINE_PHS = False                   # Flag to determine if we examine phase-space mixing and produce histograms
PHS_INT_S = 3600                      # How many seconds to wait between calculating phase-space coords,
PHS_INT = int(PHS_INT_S / DT)        # and how many time steps.
PHS_FRAMES = int(NTIME / PHS_INT)    # How many frames of our phase space plot we'll have.

DCLOSE_METRE = 1000.        # track if closer than this
LINKED_THRESH = 5.          # Count objects initialized closer than 5m together as the same object

PLOT = True     # If we want to create the plots directly in satEvol

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

# Nerd stuff
F_CHUNK = int(50)   # How many time steps we wait between dumping to outfile

# Helper I/O functions
# Read infile of TLEs and return a structured array of satellites
def readTLEs(path: str, phase_outfile: bool = False, check_linked: bool = True) -> tuple:
    if phase_outfile: woh = open(PHASE_FOUT,"w")

    with open(path,'r') as f:
        lines = f.readlines()
        NSat = int(len(lines)/3) # Total number of satellites, assuming three lines per

        # names = np.zeros(NSat, dtype='U10')
        # a = np.zeros(NSat, dtype='f8')
        # ma = np.zeros(NSat, dtype='f8')
        # omega = np.zeros(NSat, dtype='f8')
        # Omega = np.zeros(NSat, dtype='f8')
        # Omega_dot = np.zeros(NSat, dtype='f8')
        # n = np.zeros(NSat, dtype='f8')
        # ecc = np.zeros(NSat, dtype='f8')
        # inc = np.zeros(NSat, dtype='f8')

        names, a, ma, omega, Omega, Omega_dot, n, ecc, inc = [], [], [], [], [], [], [], [], []

        satIndex = 0
        for i, line in enumerate(lines):
            if len(line) < 1: break
            if line[0] == '0':
                id, junk, sname = line.rstrip().partition(' ')

            elif line[0] == '1':
                s = line.rstrip()

            elif line[0] == '2':
                t = line.rstrip()
                sat = sgp4.Satrec.twoline2rv(s, t)

                if sat.jdsatepoch < JD - JDOFFSET:
                    print("ISSUE WITH JDs: WANTING {} GOT {}".format(JD, sat.jdsatepoch))
                    continue

                err, r_km, v_km = sat.sgp4(JD, FR)
                v = np.array([v_km[0],v_km[1],v_km[2]])*SMASCALE
                r = np.array([r_km[0],r_km[1],r_km[2]])*SMASCALE
                sat_a, sat_ecc, sat_omega, sat_inc, sat_Omega, sat_nu = KT.getORBELM(r, v, muE)

                sat_ma = KT.meanAnom(sat_ecc, sat_nu)  # Mean anomaly
                sat_n = np.sqrt(G*MEarth/sat_a**3) # Mean angular motion
                sat_Omega_dot = -1.5*(REarth)**2 / (sat_a*(1 - sat_ecc**2))**2 * J2 * sat_n * np.cos(sat_inc) # Precession rate

                if sat_a*(1 - sat_ecc) < P_THRESH:
                    # names[satIndex] = sname
                    # a[satIndex] = sat_a
                    # ma[satIndex] = sat_ma
                    # omega[satIndex] = sat_omega
                    # Omega[satIndex] = sat_Omega
                    # Omega_dot[satIndex] = sat_Omega_dot
                    # n[satIndex] = sat_n
                    # ecc[satIndex] = sat_ecc
                    # inc[satIndex] = sat_inc

                    names.append(sname)
                    a.append(sat_a)
                    ma.append(sat_ma)
                    omega.append(sat_omega)
                    Omega.append(sat_Omega)
                    Omega_dot.append(sat_Omega_dot)
                    n.append(sat_n)
                    ecc.append(sat_ecc)
                    inc.append(sat_inc)

                    if phase_outfile: woh.write("{},{},{},{},{},{},{}\n".format(sname, r_km[0], r_km[1], r_km[2], v_km[0], v_km[1], v_km[2]))
                satIndex += 1

    pos, vel = KT.getXYZVVV(np.array(ma), np.array(a), np.array(omega), np.array(ecc), 
                            np.array(Omega), np.array(inc), np.array(n))
    
    pos = pos.tolist()
    vel = vel.tolist()

    if phase_outfile: woh.close()
    if check_linked:
        print('Checking for linked objects...')
        linkedIndices = threshQuery(pos, LINKED_THRESH, pairwise=False)
        removedIndices = []

        for i in range(len(linkedIndices)):
            if len(linkedIndices[i]) == 1 or i in removedIndices: continue
            for j in linkedIndices[i]:
                names[i] += '_' + names[j]

                del names[j]
                del a[j]
                del ma[j]
                del omega[j]
                del Omega[j]
                del Omega_dot[j]
                del n[j]
                del ecc[j]
                del inc[j]
                del pos[j]
                del vel[j]

                removedIndices.append(j)

        print('Number of satellites linked: {}'.format(len(removedIndices)))

    # Convert everything back to arrays to feed into class
    names = np.array(names, dtype='U10')
    a = np.array(a, dtype='f8')
    ma = np.array(ma, dtype='f8')
    omega = np.array(omega, dtype='f8')
    Omega = np.array(Omega, dtype='f8')
    Omega_dot = np.array(Omega_dot, dtype='f8')
    n = np.array(n, dtype='f8')
    ecc = np.array(ecc, dtype='f8')
    inc = np.array(inc, dtype='f8')
    pos = np.array(pos, dtype=('f8', 3))
    vel = np.array(vel, dtype=('f8', 3))

    return names, a, ma, omega, Omega, Omega_dot, n, ecc, inc, pos, vel

def writeOutfile(conjList: tuple, fname: str, type='close-approach') -> None:
    with open(fname, 'a', newline='') as f:
        for i in range(len(conjList)):
            conj = conjList[i]._asdict()
            writer = csv.writer(f, delimiter=',')
            writer.writerows(zip(conj['t'], conj['dist'], conj['vel'], conj['alt'],
                                 conj['id1'], conj['id2'], conj['name1'], conj['name2']))

        
# Distance searching function (can't be jitted because of SciPy)
def threshQuery(pos: np.ndarray, d: float, pairwise: bool = True):
    posTree = KDTree(pos)
    if pairwise: return posTree.query_pairs(d, output_type='ndarray')
    else: return posTree.query_ball_tree(posTree, d)

# ================= MAIN BLOCK ================= #

def main() -> None:
    # Read in TLEs
    print('Reading TLEs...')
    satData = readTLEs(TLES_FILENAME, phase_outfile=True, check_linked=True)

    Satellites = satArray(satData[0], satData[1], satData[2], satData[3], satData[4], satData[5],
                          satData[6], satData[7], satData[8], satData[9], satData[10])  # For some reason this doesn't construct if I just pass the tuple in
    
    print('Total # of satellites: {}'.format(Satellites.NSat))

    # Initialize satellite positions
    Satellites.updateKinematics()

    simTime = 0
    conjunctionData = []
    
    # Main integration
    print('Beginning main integration... \nTotal time: {}s'.format(TTIME))
    for itime in range(NTIME):
        Satellites.updateOrbit(DT)

        closeIndices = threshQuery(Satellites.pos, DCLOSE_METRE)
        conjunctionData.append(Satellites.getConjunctionData(closeIndices))

        if itime % F_CHUNK == 0: 
            writeOutfile(conjunctionData, OUTFILE)
            conjunctionData.clear()

        print(simTime)
        simTime += DT


if __name__ == '__main__':
    main()
    quit()
        
simtime = 0.

plotProp = {
    't':     [],
    'dist':  [],
    'vel':   [],
    'alt':   [],
    'id1':   [],
    'id2':   [],
    'name1': [],
    'name2': []
}

hist_dist = []      # List of arrays of nearest neighbour distances
hist_vel = []       # List of arrays of nearest neighbour relative velocities

fh=open("./out/track_out" + OUT_SUFFIX +".dat", "w")

# begin main integration
print('Beginning main integration... \nTotal time: {}s'.format(TTIME))
for itime in range(NTIME):
    # advance mean anomaly one step
    sat_ma = sat_ma + sat_n*DT
    flag = sat_ma > twopi
    sat_ma[flag] -= twopi

    # advance node one step
    sat_Omega = sat_Omega + Omega_dot*DT
    flag = sat_Omega > twopi
    sat_Omega[flag] -= twopi

    # NOTE: Get rid of this loop once you vectorize KeplerTools
    for i in range(NSAT):
        x0,y0,z0,vx0,vy0,vz0 = KT.getXYZVVV(sat_ma[i], sat_a[i], sat_omega[i], sat_e[i], sat_Omega[i], sat_I[i], m0=MEarth, m1=0.)
        x[i]=x0
        y[i]=y0
        z[i]=z0
        vx[i]=vx0
        vy[i]=vy0
        vz[i]=vz0

    rsphere = np.sqrt(x*x+y*y+z*z)
    phiSat = np.arctan2(y,x)
    thetaSat = np.arccos(z/rsphere)

    # Assemble a KDTree of satellite positions to query for satellites within the distance threshold
    satPos = np.column_stack((x, y, z)) # Array where the ith element contains the 3d coordinates of the ith satellite
    satVel = np.column_stack((vx, vy, vz))

    # If we want to examine phase-space mixing, calculate the nearest neighbour distance and relative velocity of every satellite every PHASE_INT
    # time steps, otherwise keep track of close approach tracks under DCLOSE_METRE.
    if EXAMINE_PHS and itime % PHS_INT == 0:
        satPosTree = KDTree(satPos) # Only construct this on time-steps we have to
        NNd, NNi = satPosTree.query(satPos, k=2, workers=4)
        relD, relV = np.zeros(NSAT), np.zeros(NSAT)

        for i in range(len(NNi)):
            i1, i2 = NNi[i][0], NNi[i][1]

            relD[i] = np.linalg.norm(satPos[i1] - satPos[i2])
            relV[i] = np.linalg.norm(satVel[i1] - satVel[i2])

        hist_dist.append(relD)
        hist_vel.append(relV)
    else:
        satPosTree = KDTree(satPos)
        pairs = satPosTree.query_pairs(DCLOSE_METRE, output_type='ndarray') # Array of pairs of indices for satellites below the distance threshold
        closeN = len(pairs)

        # Loop over that array to pull out information from the corresponding satellites
        for i in range(len(pairs)):
            i1, i2 = pairs[i, 0], pairs[i, 1] # Pull out both indices

            relD = np.linalg.norm(satPos[i1] - satPos[i2]) # Distance
            relV = np.linalg.norm(satVel[i1] - satVel[i2]) # Velocity

            alt = rsphere[i1] - REarth

            plotProp['t'].append(simtime)
            plotProp['dist'].append(relD)
            plotProp['vel'].append(relV)

            plotProp['alt'].append(alt) # Altitudes

            plotProp['id1'].append(i1)
            plotProp['id2'].append(i2)
            plotProp['name1'].append(sat_sname[i1])
            plotProp['name2'].append(sat_sname[i2])

        if itime % F_CHUNK == 0:   
            writeOutfile(plotProp)

            for k, v in plotProp.items():
                plotProp[k].clear() # Reset the lists to clear space
    
    print("Time: {time}s \t Satellites within close-approach distance: {n}".format(time = str(simtime)[0:7], n = closeN), end='\r')

    simtime+=DT

hist_dist = np.array(hist_dist)
hist_vel = np.array(hist_vel)

print('\nWriting outfile/Plotting...')

writeOutfile(plotProp)  # Write one last time to catch the last time-steps

fh.close()

if EXAMINE_PHS and PLOT:
    # for iFrame in range(PHS_FRAMES):
    #     fig, ax = plt.subplots(1, 1, figsize=(10,10))

    #     ax.hist2d(hist_dist[iFrame], hist_vel[iFrame], bins = 100)
    #     plt.savefig('./out/phs_hist/hist_frame_%s.png' % (iFrame), format='png', dpi=300, facecolor='white')

    for iFrame in range(PHS_FRAMES):
        fig, (axD, axV) = plt.subplots(1, 2, figsize=(20,10))

        axD.hist(hist_dist[iFrame]/SMASCALE, bins=50, range=[0, 800])
        axV.hist(hist_vel[iFrame]/SMASCALE, bins=50, range=[0, 25])

        axD.set_xlabel('Nearest Neighbour Distance [km]')
        axV.set_xlabel('Nearest Neighbour Relative Velocity [km/s]')

        axD.set_ylim(0, 800)
        axV.set_ylim(0, 2000)

        axD.set_ylabel('Frequency')

        plt.savefig('./out/phs_hist/hist_frame_%s.png' % (iFrame), format='png', dpi=300, facecolor='white')
        plt.close()

else:
    phiSat=phiSat*180/np.pi
    thetaSat=90-thetaSat*180/np.pi

    flag = phiSat>180
    phiSat[flag]=phiSat[flag]-360

    if PLOT:
        plt.figure()
        plt.scatter(np.array(plotProp['t'])/60,np.array(plotProp['dist'])/1000,s=1)
        plt.title("Close approach tracks")
        plt.ylabel("Close Approach Distance [km]")
        plt.xlabel("Time (minutes)")
        plt.savefig("./out/close_approach_tracks" + OUT_SUFFIX + ".pdf")

        fig=plt.figure(figsize=[15,8])
        #ax=fig.add_subplot(1,1,1, projection=ccrs.Robinson())
        ax=fig.add_subplot(1,1,1, projection=ccrs.Mollweide())
        #ax.set_global()
        ax.stock_img()
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS)
        ax.gridlines()

        poly = ax.scatter(phiSat,thetaSat,s=0.1,transform=ccrs.PlateCarree(),c='black',alpha=0.75)
        plt.title("Simulated Satellites Projected onto Earth")
        plt.savefig("./out/sats_mollweide" + OUT_SUFFIX + ".pdf")


        fig=plt.figure(figsize=[15,8])
        ax=fig.add_subplot(1,1,1, projection=ccrs.Mollweide())
        ax.gridlines()


        poly = ax.scatter((-phiSat),thetaSat,s=0.2,transform=ccrs.PlateCarree(),c='black')
        plt.title("Simulated Satellites on CS")
        plt.savefig("./out/sats_mollweide_CS" + OUT_SUFFIX + ".pdf")



        fig=plt.figure(figsize=[15,8])
        #ax=fig.add_subplot(1,1,1, projection=ccrs.Robinson())
        ax=fig.add_subplot(1,1,1, projection=ccrs.Mollweide())
        #ax.set_global()
        #ax.stock_img()
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS)
        ax.gridlines()

        poly = ax.scatter(phiSat,thetaSat,s=0.2,transform=ccrs.PlateCarree(),c='blue',alpha=1.00)
        plt.title("Simulated Satellites Projected onto Earth", fontsize=20)
        plt.savefig("./out/sats_65k_mollweide_lineonly" + OUT_SUFFIX + ".pdf")

print('Done!')
