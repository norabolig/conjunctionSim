## Written by Aaron Boley. June 2021
## Super klunky. 
## Awesome licence stuff 
## No guarantees

## Made even cooler by Skye Heiland. February 2025

import KeplerTools as KT
import numpy as np
import matplotlib.pylab as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sgp4.api as sgp4

# Select TLE Catalogue
# Also set epoch information for TLEs
#TLES_FILENAME="starlink_07OCT2024.tles"
TLES_FILENAME="./in/starlink_deb.tles"
PHASE_FOUT="./out/catdata_xyzvxvyvz.dat"
JD=2460591.5
FR=0.0
JDOFFSET=100
EPOCHLIM = 24280


# constants and parameters
twopi=np.pi*2
MEarth = 5.97e24
REarth = 6378.135e3
REkm = 6378.135
ECCSCALE=0.0002 # make eccentric enough to fill shells for any artificial systems
SMASCALE=1000 # km to metres
DT=0.05 # time step in seconds
NTIME=18000 # number of steps
DCLOSE_METRE=10000 # track if closer than this
P_THRESH=2000e3 + REarth # do not include objects with pericentres above this


aukm=1.496e8
au=aukm*1e3
G=6.6743e-11
J2=1082.64e-6
muE = 3.986004418e14

sat_sname=[]
sat_a=[]
sat_ma=[]
sat_omega=[]
sat_Omega=[]
sat_e=[]
sat_I=[]


# read tle file, open log file
tles_fh = open(TLES_FILENAME,"r")
woh = open(PHASE_FOUT,"w")

while tles_fh:
    calculate=0
    line=tles_fh.readline()
    print(line.rstrip(),len(line))
    if len(line)<1: break
    if line[0]=="0":
       id,junk,sname=line.rstrip().partition(" ")
    elif line[0]=="1":
       s=line.rstrip()
    elif line[0]=="2":
       t=line.rstrip()
       calculate=1
#    else:
#       except(typeError): print("Warning -- TLE LINE ID NOT FOUND")

    if calculate==1:
        satellite = sgp4.Satrec.twoline2rv(s, t)
        jd_last = satellite.jdsatepoch
        if jd_last < JD - JDOFFSET: 
              print("ISSUE WITH JDs: WANTING {} GOT {}".format(JD,jd_last))
              continue 
        err,r_km,v_km = satellite.sgp4(JD, FR)
        v=np.array([v_km[0],v_km[1],v_km[2]])*1000
        r=np.array([r_km[0],r_km[1],r_km[2]])*1000

        print(r)

        a,ecc,omega,inc,Omega,nu = KT.getORBELM(r,v,muE)
        print(a,ecc,omega,inc,Omega,nu)

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

woh.close() 

sat_a=np.array(sat_a)
sat_omega=np.array(sat_omega)
sat_Omega=np.array(sat_Omega)
sat_e=np.array(sat_e)
sat_I=np.array(sat_I)
sat_ma=np.array(sat_ma)
NSAT=len(sat_a)

sat_n = np.sqrt(G*MEarth/sat_a**3)

# set precession rate
Omega_dot = -1.5* (REarth)**2/(sat_a*(1-sat_e))**2*J2*sat_n*np.cos(sat_I)

# sanity check
print("Total sats: {}".format(NSAT))

x=np.zeros(NSAT)
y=np.zeros(NSAT)
z=np.zeros(NSAT)
vx=np.zeros(NSAT)
vy=np.zeros(NSAT)
vz=np.zeros(NSAT)

simtime=0.
plot_time=[]
plot_dist=[]
plot_range=[]
plot_id1=[]
plot_id2=[]
plot_name1=[]
plot_name2=[]
plot_vel=[]

# begin main integration
for itime in range(NTIME):
    print("TIME {}".format(simtime))
   
    # advance mean anomaly one step
    sat_ma=sat_ma+sat_n*DT
    flag = sat_ma>twopi
    sat_ma[flag]-=twopi

    # advance node one step
    sat_Omega = sat_Omega + Omega_dot*DT
    flag = sat_Omega > twopi
    sat_Omega[flag]-=twopi


    for i in range(NSAT):
        x0,y0,z0,vx0,vy0,vz0 = KT.getXYZVVV(sat_ma[i],sat_a[i],sat_omega[i],sat_e[i],sat_Omega[i],sat_I[i],m0=MEarth,m1=0.)
        x[i]=x0
        y[i]=y0
        z[i]=z0
        vx[i]=vx0
        vy[i]=vy0
        vz[i]=vz0
#    x,y,z,vx,vy,vz = KT.getXYZVVV(sat_ma,sat_a,sat_omega,sat_e,sat_Omega,sat_I,m0=MEarth,m1=0.)




    rsphere = np.sqrt(x*x+y*y+z*z)
    phiSat = np.arctan2(y,x)
    thetaSat = np.arccos(z/rsphere)

    pairs=[]

    # now bookkeeping. Work through satellies and find close pairs. Log it.
    for i in range(NSAT):

        #if not "STARLINK" in sat_sname[i]:continue

        d = np.sqrt( (x[i]-x)**2+(y[i]-y)**2+(z[i]-z)**2)
        flag = d < DCLOSE_METRE
        indices = np.nonzero(flag)
        #print(indices[0])
        #ncand=len(indices[0])-1
        #if ncand>0:print("NUMBER OF CANDIDATES {}".format(len(indices[0])))

        for j in indices[0]:
           #print(i,j)
           if j==i: continue 
           thisPair="{}+{}".format(j,i)
           #print(thisPair)
           try:
             tagged= np.char.equal(pairs,thisPair)
             #print(tagged)
           except:
             tagged=[False]
           if not np.any(tagged):
             v = np.sqrt( (vx[i]-vx[j])**2 + (vy[i]-vy[j])**2 + (vz[i]-vz[j])**2 )
             pairs.append("{}+{}".format(i,j))
             #print("CLOSE APPROACH {} metres btw {} and {}".format(d,i,j))
             plot_time.append(simtime)
             plot_dist.append(d[j])
             plot_vel.append(v)
             plot_range.append(rsphere[i]-REarth)
             plot_id1.append(i)
             plot_id2.append(j)
             plot_name1.append(sat_sname[i])
             plot_name2.append(sat_sname[j])

    simtime+=DT

plot_time=np.array(plot_time)
plot_dist=np.array(plot_dist)

plt.figure()
plt.scatter(plot_time/60,plot_dist/1000,s=1)
plt.title("Close approach tracks")
plt.ylabel("Close Approach Distance [km]")
plt.xlabel("Time (minutes)")
plt.savefig("close_approach_tracks.pdf")

fh=open("./out/track_out.dat","w")
for i in range(len(plot_time)):
    fh.write("{},{},{},{},{},{},{},{}\n".format(plot_time[i],plot_dist[i],plot_vel[i],plot_range[i],plot_id1[i],plot_id2[i],plot_name1[i],plot_name2[i]))
fh.close()

phiSat=phiSat*180/np.pi
thetaSat=90-thetaSat*180/np.pi

flag = phiSat>180
phiSat[flag]=phiSat[flag]-360

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
plt.savefig("./out/sats_mollweide.pdf")


fig=plt.figure(figsize=[15,8])
ax=fig.add_subplot(1,1,1, projection=ccrs.Mollweide())
ax.gridlines()


poly = ax.scatter((-phiSat),thetaSat,s=0.2,transform=ccrs.PlateCarree(),c='black')
plt.title("Simulated Satellites on CS")
plt.savefig("./out/sats_mollweide_CS.pdf")



fig=plt.figure(figsize=[15,8])
#ax=fig.add_subplot(1,1,1, projection=ccrs.Robinson())
ax=fig.add_subplot(1,1,1, projection=ccrs.Mollweide())
#ax.set_global()
#ax.stock_img()
ax.coastlines()
ax.add_feature(cfeature.BORDERS)
ax.gridlines()

poly = ax.scatter(phiSat,thetaSat,s=0.2,transform=ccrs.PlateCarree(),c='blue',alpha=1.00)
plt.title("Simulated Satellites Projected onto Earth",fontsize=20)
plt.savefig("./out/sats_65k_mollweide_lineonly.pdf")





plt.show()





