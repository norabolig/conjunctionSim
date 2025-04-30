import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm

CLOSE=1e3  # CHOOSE CLOSE APPROACH DISTANCE HERE
DTMAX=0.05
TMIN=0
TMAX=5*60  # CHOOSE MAX TIME FOR PLOT HERE
DTSEP=1000
PLOTSTAR=False
OUTPUTENCOUNTER=True  # CHOOSE IF YOU WANT TO SAVE ENCOUNTER DATA

FILENAME = 'out/track_close_encounters_7day_10km.dat'  # PUT IN YOUR OUTPUT FILE HERE
OUTPUTFILE = 'out/output_track_close_encounters_7day_10km.hdf'  # FOR IF YOU WANT TO SAVE ENCOUNTER DATA

columns = ['time','dist','dv', 'alt', 'index1', 'index2', 'name1', 'name2']
dat = pd.read_csv(FILENAME, delimiter=',', header=None)  # PUT IN YOUR OUTPUT FILE HERE
dat.columns=columns
subset = dat.loc[(dat.dist>0)&(dat.dist<CLOSE)&(dat.time>TMIN)&(dat.time<TMAX)]

time=np.array(subset.time)
dist=np.array(subset.dist)
iid=np.array(subset.index1)
jid=np.array(subset.index2)
v=np.array(subset.dv)
alt=np.array(subset.alt)
iname=np.array(subset.name1)
jname=np.array(subset.name2)

#========================================================================================

dclose=[]
vclose=[]
altclose=[]
tclose=[]
cclose=[]
DTCLOSE=np.array([])

sum=0
rbsatsum=0
debsatsum=0
rbdebsum=0
satsatsum=0
rbstarsum=0
debstarsum=0
satstarsum=0
rbrbsum=0
debdebsum=0

plt.figure(figsize=(12,4))

while True:
    if time.size<1: break
    id=iid[0]
    jd=jid[0]
    flag1 = iid==id
    flag2 = jid==jd
    flag = flag1*flag2
            
    min_d0 = min(dist[flag])
    
    dt = np.append(np.zeros(1), np.array((time[flag][1:] - time[flag][:-1])))
    DTCLOSE = np.append(DTCLOSE, dt[dt>DTSEP])
    nclose = len(dt[dt>DTSEP])+1
    #print('Nclose: ', nclose)

    if nclose == 1:
        sum += 1
        min_d = min(dist[flag])
        id = np.argmin(dist[flag])
        dclose.append(min_d)
        vclose.append(v[flag][id])
        altclose.append(alt[flag][id])
        tclose.append(time[flag][id])
    elif nclose >= 2:
        sum += 1
        tsep = np.where(dt>DTSEP)[0][0]
        flag[np.where(flag==True)[0][tsep:]]=False
        min_d = min(dist[flag])
        id = np.argmin(dist[flag])
        dclose.append(min_d)
        vclose.append(v[flag][id])
        altclose.append(alt[flag][id])
        tclose.append(time[flag][id])

    if iname[flag][id].find("DEB")>-1:
        if jname[flag][id].find("DEB")>-1:
            carray="green"
            ls='-'
            debdebsum+=1
        elif jname[flag][id].find("R/B")>-1:
            carray="purple"
            ls='-'
            rbdebsum+=1
        elif jname[flag][id].find("STARLINK")>-1:
            carray="blue"
            ls='--'
            debstarsum+=1
            debsatsum+=1
        else:
            carray="blue"
            ls='-'
            debsatsum+=1

    elif iname[flag][id].find("R/B")>-1:
        if jname[flag][id].find("DEB")>-1:
            carray="purple"
            ls='-'
            rbdebsum+=1
        elif jname[flag][id].find("R/B")>-1:
            carray="pink"
            ls='-'
            rbrbsum+=1
        elif jname[flag][id].find("STARLINK")>-1:
            carray="orange"
            ls='--'
            rbsatsum+=1
            rbstarsum+=1
        else:
            carray="orange"
            ls='-'
            rbsatsum+=1

    elif iname[flag][id].find("STARLINK")>-1:
        if jname[flag][id].find("DEB")>-1:
            carray="blue"
            ls='--'
            debstarsum+=1
            debsatsum+=1
        elif jname[flag][id].find("R/B")>-1:
            carray="orange"
            ls='--'
            rbstarsum+=1
            rbsatsum+=1
        else:
            carray="red"
            ls='--'
            satsatsum+=1
            satstarsum+=1

    else:
        if jname[flag][id].find("DEB")>-1:
            carray="blue"
            ls='-'
            debsatsum+=1
        elif jname[flag][id].find("R/B")>-1:
            carray="orange"
            ls='-'
            rbsatsum+=1
        elif jname[flag][id].find("STARLINK")>-1:
            carray="red"
            ls='--'
            satstarsum+=1
            satsatsum+=1
        else:
            carray="red"
            ls='-'
            satsatsum+=1

    cclose.append(carray)
    
    if PLOTSTAR:
        plt.plot(time[flag]/60,dist[flag]/1e3,ls=ls,color=carray)
    else:
        plt.plot(time[flag]/60,dist[flag]/1e3,color=carray)
        
    time=np.delete(time,flag)
    dist=np.delete(dist,flag)
    iid=np.delete(iid,flag)
    jid=np.delete(jid,flag)
    v=np.delete(v,flag)
    alt=np.delete(alt,flag)
    iname=np.delete(iname,flag)
    jname=np.delete(jname,flag)
    
dclose = np.array(dclose)
vclose = np.array(vclose)
altclose = np.array(altclose)
tclose = np.array(tclose)
cclose = np.array(cclose)

DF = pd.DataFrame(np.array([tclose,dclose,vclose,altclose,cclose]).T)
DF.columns = ['time','dist','dV','alt','type']
DF.to_hdf(OUTPUTFILE,key='data')


custom_lines = [Line2D([0], [0], color='green', lw=4),
                Line2D([0], [0], color='blue', lw=4),
                Line2D([0], [0], color='red', lw=4),
                Line2D([0], [0], color='orange', lw=4),
                Line2D([0], [0], color='purple', lw=4),
                Line2D([0], [0], color='pink', lw=4)]
plt.legend(custom_lines, ['DEB-DEB: {}'.format(debdebsum), 'SAT-DEB: {}({})'.format(debsatsum,debstarsum), 
                          'SAT-SAT: {}({})'.format(satsatsum,satstarsum),
                          'SAT-R/B: {}({})'.format(rbsatsum,rbstarsum),
                          'DEB-R/B: {}'.format(rbdebsum),
                          'R/B-R/B: {}'.format(rbrbsum)], ncol=3)
plt.xlabel("Time in minutes")
plt.ylabel("Conjunction Distance [km]")
plt.title('Total number of close encounters < {} km is {}'.format(CLOSE/1e3,sum))
plt.savefig("./out/close_approach_tracks_fit.pdf", bbox_inches="tight")



print("Total detected approaches < {} km is {}".format(CLOSE/1e3,sum))

#========================================================================================

fig, ax = plt.subplots(2,3, figsize=(10,7))
ax[0,0].hist2d(dclose[cclose=='green']/1e3,vclose[cclose=='green']/1e3,bins=20)
ax[0,0].set_xlabel('min. close approach distance [km]')
ax[0,0].set_ylabel('dv at min. close approach [km/s]')
ax[0,0].set_title('DEB-DEB')

ax[0,1].hist2d(dclose[cclose=='blue']/1e3,vclose[cclose=='blue']/1e3,bins=20)
ax[0,1].set_xlabel('min. close approach distance [km]')
ax[0,1].set_ylabel('dv at min. close approach [km/s]')
ax[0,1].set_title('SAT-DEB')

ax[0,2].hist2d(dclose[cclose=='red']/1e3,vclose[cclose=='red']/1e3,bins=20)
ax[0,2].set_xlabel('min. close approach distance [km]')
ax[0,2].set_ylabel('dv at min. close approach [km/s]')
ax[0,2].set_title('SAT-SAT')

ax[1,0].hist2d(dclose[cclose=='orange']/1e3,vclose[cclose=='orange']/1e3,bins=20)
ax[1,0].set_xlabel('min. close approach distance [km]')
ax[1,0].set_ylabel('dv at min. close approach [km/s]')
ax[1,0].set_title('SAT-R/B')

ax[1,1].hist2d(dclose[cclose=='purple']/1e3,vclose[cclose=='purple']/1e3,bins=20)
ax[1,1].set_xlabel('min. close approach distance [km]')
ax[1,1].set_ylabel('dv at min. close approach [km/s]')
ax[1,1].set_title('DEB-R/B')

ax[1,2].hist2d(dclose[cclose=='pink']/1e3,vclose[cclose=='pink']/1e3,bins=20)
ax[1,2].set_xlabel('min. close approach distance [km]')
ax[1,2].set_ylabel('dv at min. close approach [km/s]')
ax[1,2].set_title('R/B-R/B')

plt.subplots_adjust(wspace=0.3,hspace=0.4)
plt.savefig("./out/relative_dist_vel.pdf", bbox_inches="tight")

#========================================================================================

plt.figure()
h,bins,_ = plt.hist(vclose[cclose=='red']/1e3,bins=50);
plt.xlabel('dv at min. close approach [km/s]')
plt.title('SAT-SAT, < 1km')
plt.savefig("./out/vclose_sat_hist.pdf", bbox_inches="tight")

#========================================================================================

plt.figure()
plt.hist2d(tclose/60/60, altclose/1e3, norm=LogNorm(), bins=30);
plt.colorbar(label='No. of close encounters < {} km'.format(int(CLOSE/1e3)))
plt.xlabel('Time [hr]')
plt.ylabel('Altitude [km]')
plt.savefig("./out/time_alt_encounters.pdf", bbox_inches="tight")

#========================================================================================

plt.figure()
bins = np.linspace(0,TMAX/60/60,30)
plt.hist(tclose[cclose=='green']/60/60, color='green',histtype='step', bins=bins);
plt.hist(tclose[cclose=='blue']/60/60, color='blue',histtype='step',bins=bins);
plt.hist(tclose[cclose=='red']/60/60, color='red',histtype='step',bins=bins);
plt.hist(tclose[cclose=='orange']/60/60, color='orange',histtype='step',bins=bins);
plt.hist(tclose[cclose=='purple']/60/60, color='purple',histtype='step',bins=bins);
plt.hist(tclose[cclose=='pink']/60/60, color='pink',histtype='step',bins=bins);
custom_lines = [Line2D([0], [0], color='green', lw=4),
                Line2D([0], [0], color='blue', lw=4),
                Line2D([0], [0], color='red', lw=4),
                Line2D([0], [0], color='orange', lw=4),
                Line2D([0], [0], color='purple', lw=4),
                Line2D([0], [0], color='pink', lw=4)]
plt.legend(custom_lines, ['DEB-DEB: {}'.format(debdebsum), 'SAT-DEB: {} ({})'.format(debsatsum,debstarsum), 
                          'SAT-SAT: {} ({})'.format(satsatsum,satstarsum),
                          'SAT-R/B: {} ({})'.format(rbsatsum,rbstarsum),
                          'DEB-R/B: {}'.format(rbdebsum),
                          'R/B-R/B: {}'.format(rbrbsum)], ncol=1, loc='upper left')
plt.xlabel('Time of a close encounter < {} km [hr]'.format(int(CLOSE/1e3)))
plt.ylabel('Counts')
plt.xlim(0,TMAX/60/60)
plt.savefig("./out/time_encounters_hist.pdf", bbox_inches="tight")

#========================================================================================

plt.figure()
altbins = np.arange(200,2000,1)
h,bins = np.histogram(altclose/1e3,bins=altbins)
binmids = (bins[:-1]+bins[1:])/2
plt.plot(binmids,TMAX/(60*h))
plt.xlim(200,2000)
plt.axhline(TMAX/60,ls='--',label='max sim. time')
plt.axhline(20,color='xkcd:red',ls='--',label='20 minutes')
plt.yscale('log')
plt.legend()
plt.xlabel('Altitude [km]')
plt.ylabel('Avg. time between close encounters < {} km [min]'.format(int(CLOSE/1e3)))
plt.savefig("./out/frequency.pdf", bbox_inches="tight")
