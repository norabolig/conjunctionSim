## Written by Aaron Boley. June 2021
## Super klunky. 
## Awesome licence stuff 
## No guarantees

## Made even cooler by Skye Heiland and Sarah Thiele. Spring 2025

import KeplerTools as KT
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
import sgp4.api as sgp4

from scipy.spatial import KDTree

# Select TLE Catalogue
# Also set epoch information for TLEs
TLES_FILENAME = "./in/full_cat.tles"
OUT_SUFFIX = '_full_7day_10km_cut'
PHASE_FOUT = "./out/catdata_xyzvxvyvz" + OUT_SUFFIX +".dat"
JD = 2460591.5
FR = 0.0
JDOFFSET = 100
EPOCHLIM = 24280

# Constants and parameters
ECCSCALE = 0.0002         # Make eccentric enough to fill shells for any artificial systems
SMASCALE = 1000           # km to metres
DT = 0.05                 # Time step in seconds
TTIME = 604800            # How long to run sim (s)
NTIME = int(TTIME / DT)   # Number of steps

EXAMINE_PHS = True                   # Flag to determine if we examine phase-space mixing and produce histograms
PHS_INT_S = 86400                    # How many seconds to wait between calculating phase-space coords,
PHS_INT = int(PHS_INT_S / DT)        # and how many time steps.
PHS_FRAMES = int(NTIME / PHS_INT)    # How many frames of our phase space plot we'll have.

DCLOSE_METRE = 15000        # track if closer than this

PLOT = False     # If we want to create the plots directly in satEvol

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
F_CHUNK = int(5e4)   # How many time steps we wait between dumping to outfile
NEWFILE_CHUNK = 150000  # How many time steps we wait before starting a new output file
NEWFILE_INDEX = 1  # Output file indexing

sat_sname=[]    # Satellite name
sat_a=[]        # Semi-major axis
sat_ma=[]       # Mean anomaly
sat_omega=[]    # Angle between something I need to figure out
sat_Omega=[]    # Angular velocity
sat_e=[]        # Eccentricity
sat_I=[]        # Inclination

# Helper functions
def writeOutfile(prop: dict, type='close-approach') -> None:
    for i in range(len(prop['t'])):
        fh.write("{},{},{},{},{},{},{},{}\n".format(prop['t'][i], prop['dist'][i], prop['vel'][i], 
                                                    prop['alt'][i], prop['id1'][i], prop['id2'][i], prop['name1'][i], prop['name2'][i]))

# read tle file, open log file
tles_fh = open(TLES_FILENAME,"r")
woh = open(PHASE_FOUT,"w")

while tles_fh:
    line=tles_fh.readline()
    # print(line.rstrip(),len(line))

    if len(line)<1: break
    if line[0]=="0": 
        id, junk, sname = line.rstrip().partition(" ")
    elif line[0]=="1": 
        s = line.rstrip()
    elif line[0]=="2":
        t = line.rstrip()
        satellite = sgp4.Satrec.twoline2rv(s, t)
        jd_last = satellite.jdsatepoch

        if jd_last < JD - JDOFFSET: 
              print("ISSUE WITH JDs: WANTING {} GOT {}".format(JD,jd_last))
              continue 
        
        err,r_km,v_km = satellite.sgp4(JD, FR)
        v=np.array([v_km[0],v_km[1],v_km[2]])*SMASCALE
        r=np.array([r_km[0],r_km[1],r_km[2]])*SMASCALE

        a,ecc,omega,inc,Omega,nu = KT.getORBELM(r,v,muE)

        if a*(1-ecc) < P_THRESH: 

            sat_sname.append(sname)
            sat_a.append(a)
            sat_e.append(ecc)
            sat_omega.append(omega)
            sat_Omega.append(Omega)
            sat_I.append(inc)
            EA=KT.EAnom(ecc,nu)
            MA = EA - ecc*np.sin(EA)
            sat_ma.append(MA)
            woh.write("{},{},{},{},{},{},{}\n".format(sname,r_km[0],r_km[1],r_km[2],v_km[0],v_km[1],v_km[2]))
#    else:
#       except(typeError): print("Warning -- TLE LINE ID NOT FOUND")

tles_fh.close()
woh.close() 

sat_a=np.array(sat_a)
sat_omega=np.array(sat_omega)
sat_Omega=np.array(sat_Omega)
sat_e=np.array(sat_e)
sat_I=np.array(sat_I)
sat_ma=np.array(sat_ma)

woh_df = pd.read_csv(PHASE_FOUT, header=None)
woh_df.columns = ['name','x','y','z','vx','vy','vz']
woh_keep = woh_df.drop_duplicates(['x','y','z','vx','vy','vz'], keep='first').index.values
sat_a=sat_a[woh_keep]
sat_omega=sat_omega[woh_keep]
sat_Omega=sat_Omega[woh_keep]
sat_e=sat_e[woh_keep]
sat_I=sat_I[woh_keep]
sat_ma=sat_ma[woh_keep]
woh_df = []

NSAT=len(sat_a)

sat_n = np.sqrt(G*MEarth/sat_a**3)

# set precession rate
Omega_dot = -1.5*(REarth)**2/(sat_a*(1-sat_e**2))**2*J2*sat_n*np.cos(sat_I)

# sanity check
print("Total sats: {}".format(NSAT))

x=np.zeros(NSAT)
y=np.zeros(NSAT)
z=np.zeros(NSAT)
vx=np.zeros(NSAT)
vy=np.zeros(NSAT)
vz=np.zeros(NSAT)

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

        if itime % NEWFILE_CHUNK == 0:
            writeOutfile(plotProp)

            for k, v in plotProp.items():
                plotProp[k].clear() # Reset the lists to clear space
                
            fh.close()
            fh=open("./out/track_out" + OUT_SUFFIX + '_' + str(NEWFILE_INDEX) +".dat", "w")
            NEWFILE_INDEX += 1

    simtime+=DT

hist_dist = np.array(hist_dist)
hist_vel = np.array(hist_vel)

print('\nWriting outfile/Plotting...')

writeOutfile(plotProp)  # Write one last time to catch the last time-steps

fh.close()

print('moving on to plot')

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

        '''
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
        '''
print('Done!')
