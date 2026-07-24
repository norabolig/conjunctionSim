import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
import pdb
import glob

MAXSIMTIME = 60*60*24*5  # full sim run time
CLOSE=1e3 # close encounter distance 
DTMAX=0.05
TMIN=0  # min time for plot
TMAX=MAXSIMTIME  # max time for plot
DTSEP=1000
PLOT_SUFFIX='_random_starlink_5d_1km'  # plot name suffix
PLOTSTAR=False
PLOT_STARLINK = True # choose if you want to note number of starlink-specific encounters
PLOT_ANALYSIS = True  # If False you can just save the encounter data and run analysis_plots.py later
SAVE = True  # choose if you want to save plots
OUTPUTENCOUNTER=True  # choose if you want to save encounter data

FILE_LOC=r'out/Starlink_random/'
FILE_PREFIX=FILE_LOC+r'track_out_starlink_random_dt_0p05_MJD_2460591.5'
columns = ['time','dist','dv', 'alt', 'index1', 'index2', 'name1', 'name2']
NUMFILES=len(glob.glob(FILE_PREFIX+r'*.dat'))
OUTPUTFILE = FILE_LOC+'encounter_data_random_starlink_5d_1km.hdf'  # FOR IF YOU WANT TO SAVE ENCOUNTER DATA
print('num files: ', NUMFILES)
    #========================================================================================

fig, ax = plt.subplots(figsize=(12,4))

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

dclose=[]
vclose=[]
altclose=[]
tclose=[]
cclose=[]
type1close=[]
type2close=[]
DTCLOSE=np.array([])

for IFILE in range(NUMFILES):
    print('file {} of {}'.format(IFILE,NUMFILES))
    if IFILE==0:
        FILENAME = FILE_PREFIX + '.dat'
    else: 
        FILENAME = FILE_PREFIX + '_{}.dat'.format(int(IFILE))

    print('file {} of {}: {}'.format(IFILE,NUMFILES,FILENAME))
    dat = pd.read_csv(FILENAME, delimiter=',', header=None)
    dat.columns=columns
    subset = dat.loc[(dat.dist>0.0)&(dat.dist<CLOSE)&(dat.time>TMIN)&(dat.time<TMAX)]
    print('loaded file')

    time=np.array(subset.time)
    dist=np.array(subset.dist)
    iid=np.array(subset.index1)
    jid=np.array(subset.index2)
    v=np.array(subset.dv)
    alt=np.array(subset.alt)
    iname=np.array(subset.name1)
    jname=np.array(subset.name2)
        
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
            type1 = "DEB"
            if jname[flag][id].find("DEB")>-1:
                carray="green"
                ls='-'
                debdebsum+=1
                type2 = "DEB"
            elif jname[flag][id].find("R/B")>-1:
                carray="purple"
                ls='-'
                rbdebsum+=1
                type2 = "R/B"
            elif jname[flag][id].find("STARLINK")>-1:
                carray="blue"
                ls='--'
                debstarsum+=1
                debsatsum+=1
                type2 = "SL"
            else:
                carray="blue"
                ls='-'
                debsatsum+=1
                type2 = "SAT"
    
        elif iname[flag][id].find("R/B")>-1:
            type1 = "R/B"
            if jname[flag][id].find("DEB")>-1:
                carray="purple"
                ls='-'
                rbdebsum+=1
                type2 = "DEB"
            elif jname[flag][id].find("R/B")>-1:
                carray="pink"
                ls='-'
                rbrbsum+=1
                type2 = "R/B"
            elif jname[flag][id].find("STARLINK")>-1:
                carray="orange"
                ls='--'
                rbsatsum+=1
                rbstarsum+=1
                type2 = "SL"
            else:
                carray="orange"
                ls='-'
                rbsatsum+=1
                type2 = "SAT"
    
        elif iname[flag][id].find("STARLINK")>-1:
            type1 = "SL"
            if jname[flag][id].find("DEB")>-1:
                carray="blue"
                ls='--'
                debstarsum+=1
                debsatsum+=1
                type2 = "DEB"
            elif jname[flag][id].find("R/B")>-1:
                carray="orange"
                ls='--'
                rbstarsum+=1
                rbsatsum+=1
                type2 = "R/B"
            else:
                carray="red"
                ls='--'
                satsatsum+=1
                satstarsum+=1
                if jname[flag][id].find("STARLINK")>-1:
                    type2 = "SL"
                else:
                    type2 = "SAT"
    
        else:
            type1 = "SAT"
            if jname[flag][id].find("DEB")>-1:
                carray="blue"
                ls='-'
                debsatsum+=1
                type2 = "DEB"
            elif jname[flag][id].find("R/B")>-1:
                carray="orange"
                ls='-'
                rbsatsum+=1
                type2 = "R/B"
            elif jname[flag][id].find("STARLINK")>-1:
                carray="red"
                ls='--'
                satstarsum+=1
                satsatsum+=1
                type2 = "SL"
            else:
                carray="red"
                ls='-'
                satsatsum+=1
                type2 = "SAT"
    
        cclose.append(carray)
        type1close.append(type1)
        type2close.append(type2)
        
        if PLOTSTAR:
            plt.plot(time[flag]/3600,dist[flag]/1e3,ls=ls,color=carray)
        else:
            plt.plot(time[flag]/3600,dist[flag]/1e3,color=carray)
            
        time=np.delete(time,flag)
        dist=np.delete(dist,flag)
        iid=np.delete(iid,flag)
        jid=np.delete(jid,flag)
        v=np.delete(v,flag)
        alt=np.delete(alt,flag)
        iname=np.delete(iname,flag)
        jname=np.delete(jname,flag)

    print("FILE {}: Total detected approaches < {} km is {}".format(int(IFILE),CLOSE/1e3,sum))
    
custom_lines = [Line2D([0], [0], color='green', lw=4),
                Line2D([0], [0], color='blue', lw=4),
                Line2D([0], [0], color='red', lw=4),
                Line2D([0], [0], color='orange', lw=4),
                Line2D([0], [0], color='purple', lw=4),
                Line2D([0], [0], color='pink', lw=4)]
if PLOT_STARLINK:
    plt.legend(custom_lines, ['DEB-DEB: {}'.format(debdebsum), 
                              'SAT-DEB: {}({})'.format(debsatsum,debstarsum), 
                              'SAT-SAT: {}({})'.format(satsatsum,satstarsum),
                              'SAT-R/B: {}({})'.format(rbsatsum,rbstarsum),
                              'DEB-R/B: {}'.format(rbdebsum),
                              'R/B-R/B: {}'.format(rbrbsum)], ncol=3,
               fontsize=14,
               loc='lower center', bbox_to_anchor=(0.0, 1.01, 1., 0.5),
               bbox_transform=ax.transAxes)
else:
    plt.legend(custom_lines, ['DEB-DEB: {}'.format(debdebsum), 
                              'SAT-DEB: {}'.format(debsatsum), 
                              'SAT-SAT: {}'.format(satsatsum),
                              'SAT-R/B: {}'.format(rbsatsum),
                              'DEB-R/B: {}'.format(rbdebsum),
                              'R/B-R/B: {}'.format(rbrbsum)], ncol=3,
               fontsize=14,
               loc='lower center', bbox_to_anchor=(0.0, 1.01, 1., 0.5),
               bbox_transform=ax.transAxes)
plt.xlabel("Time [hrs]", fontsize=16)
plt.ylabel("Conjunction Distance [km]", fontsize=16)
plt.ylim(-0.05,CLOSE/1e3)
plt.xlim(TMIN/3600,TMAX/3600)
ax.tick_params(labelsize=16)
plt.text(0.01,0.025,'Total close encounters < {} km: {}'.format(CLOSE/1e3,sum), 
         fontsize=16, transform=ax.transAxes)
plt.savefig("./out/fullcat_1km/close_approach_tracks_fit{}.pdf".format(PLOT_SUFFIX), bbox_inches="tight")
        
dclose = np.array(dclose)
vclose = np.array(vclose)
altclose = np.array(altclose)
tclose = np.array(tclose)
cclose = np.array(cclose)
type1close = np.array(type1close)
type2close = np.array(type2close)

DF = pd.DataFrame(np.array([tclose,dclose,vclose,altclose,cclose,type1close,type2close]).T)
columnsDF = ['time','dist','dV','alt','type','type1','type2']
DF.columns = columnsDF
DF.to_hdf(OUTPUTFILE,key='data',data_columns=columnsDF,mode='a',append=True)
    
    #========================================================================================
    
if PLOT_ANALYSIS:
    bins=[np.linspace(0,CLOSE/1e3,20), np.linspace(0,16,20)]
    fig, ax = plt.subplots(2,3, figsize=(10,7))
    ax[0,0].hist2d(dclose[cclose=='green']/1e3,vclose[cclose=='green']/1e3,bins=bins)
    ax[0,0].set_title('DEB-DEB', fontsize=16)
    
    ax[0,1].hist2d(dclose[cclose=='blue']/1e3,vclose[cclose=='blue']/1e3,bins=bins)
    ax[0,1].set_title('SAT-DEB', fontsize=16)
    
    ax[0,2].hist2d(dclose[cclose=='red']/1e3,vclose[cclose=='red']/1e3,bins=bins)
    ax[0,2].set_title('SAT-SAT', fontsize=16)
    
    ax[1,0].hist2d(dclose[cclose=='orange']/1e3,vclose[cclose=='orange']/1e3,bins=bins)
    ax[1,0].set_title('SAT-R/B', fontsize=16)
    
    ax[1,1].hist2d(dclose[cclose=='purple']/1e3,vclose[cclose=='purple']/1e3,bins=bins)
    ax[1,1].set_title('DEB-R/B', fontsize=16)
    
    ax[1,2].hist2d(dclose[cclose=='pink']/1e3,vclose[cclose=='pink']/1e3,bins=bins)
    ax[1,2].set_title('R/B-R/B', fontsize=16)
    
    for i in range(2):
        for j in range(3):
            ax[i,j].set_xlabel(r'$\Delta R_{\rm min}$ [km]', fontsize=16)
            ax[i,j].set_ylabel(r'$\Delta v$ at $\Delta R_{\rm min}$ [km/s]', fontsize=16)
            ax[i,j].tick_params(labelsize=14)
            ax[i,j].set_xlim(0.,CLOSE/1e3)
    
    plt.subplots_adjust(wspace=0.35,hspace=0.45)
    if SAVE:
        plt.savefig("./out/fullcat_1km/relative_dist_vel{}.pdf".format(PLOT_SUFFIX), bbox_inches="tight")
    
    #========================================================================================
    
    fig, ax = plt.subplots(figsize=(8,6))
    bins = np.linspace(0.,16.,50)
    h,bins,_ = plt.hist(vclose[cclose=='red']/1e3,bins=bins,
                        label='SAT-SAT, < {} km'.format(CLOSE/1e3));
    plt.xlabel(r'$\Delta v$ at $\Delta R_{\rm min}$ [km/s]', fontsize=16)
    plt.ylabel('Counts', fontsize=16)
    plt.legend(fontsize=16)
    ax.tick_params(labelsize=16)
    if SAVE:
        plt.savefig("./out/fullcat_1km/vclose_sat_hist{}.pdf".format(PLOT_SUFFIX), bbox_inches="tight")
    
    #========================================================================================
    
    fig, ax = plt.subplots(figsize=(8,6))
    plt.hist2d(tclose/3600, altclose/1e3, norm=LogNorm(), bins=30);
    cb = plt.colorbar()
    cb.ax.set_ylabel('No. of close encounters < {} km'.format(int(CLOSE/1e3)), fontsize=16)
    cb.ax.tick_params(labelsize=14)
    plt.xlabel('Time [hr]', fontsize=16)
    plt.ylabel('Altitude [km]', fontsize=16)
    ax.tick_params(labelsize=16)
    if SAVE:
        plt.savefig("./out/fullcat_1km/time_alt_encounters{}.pdf".format(PLOT_SUFFIX), bbox_inches="tight")
    
    #========================================================================================
    
    fig, ax = plt.subplots(figsize=(8,6))
    bins = np.linspace(TMIN/3600,TMAX/3600,30)
    plt.hist(tclose[cclose=='green']/3600, color='green',histtype='step', bins=bins);
    plt.hist(tclose[cclose=='blue']/3600, color='blue',histtype='step',bins=bins);
    plt.hist(tclose[cclose=='red']/3600, color='red',histtype='step',bins=bins);
    plt.hist(tclose[cclose=='orange']/3600, color='orange',histtype='step',bins=bins);
    plt.hist(tclose[cclose=='purple']/3600, color='purple',histtype='step',bins=bins);
    plt.hist(tclose[cclose=='pink']/3600, color='pink',histtype='step',bins=bins);
    
    custom_lines = [Line2D([0], [0], color='green', lw=4),
                    Line2D([0], [0], color='blue', lw=4),
                    Line2D([0], [0], color='red', lw=4),
                    Line2D([0], [0], color='orange', lw=4),
                    Line2D([0], [0], color='purple', lw=4),
                    Line2D([0], [0], color='pink', lw=4)]
    plt.legend(custom_lines, ['DEB-DEB','SAT-DEB','SAT-SAT',
                              'SAT-R/B','DEB-R/B','R/B-R/B'], ncol=3,
               fontsize=14,
               loc='lower center', bbox_to_anchor=(0.0, 1.01, 1., 0.5),
               bbox_transform=ax.transAxes)
        
    plt.xlabel('Time of close encounter < {} km [hr]'.format(int(CLOSE/1e3)),fontsize=16)
    plt.ylabel('Counts',fontsize=16)
    ax.tick_params(labelsize=16)
    plt.xlim(TMIN/(3600),TMAX/(3600))

    if SAVE:
        plt.savefig("./out/fullcat_1km/time_encounters_hist{}.pdf".format(PLOT_SUFFIX), bbox_inches="tight")
    
    #========================================================================================
    
    fig, ax = plt.subplots(figsize=(8,6))
    altbins = np.arange(200,2000,1)
    h,bins = np.histogram(altclose/1e3,bins=altbins)
    binmids = (bins[:-1]+bins[1:])/2
    plt.plot(binmids,MAXSIMTIME/(60*h),color='k')
    plt.xlim(200,2000)
    plt.axhline(MAXSIMTIME/60,ls='--',color='xkcd:grey',label='max sim. time')
    plt.axhline(20,color='xkcd:red',ls='--',label='20 minutes')
    plt.yscale('log')
    plt.legend(fontsize=14)
    plt.xlabel('Altitude [km]', fontsize=16)
    ax.tick_params(labelsize=16)
    ax.set_xticks(np.arange(300,2100,300))
    plt.ylabel('Avg. close encounter rate < {} km [min]'.format(int(CLOSE/1e3)),
              fontsize=16)
    if SAVE:
        plt.savefig("./out/fullcat_1km/frequency{}.pdf".format(PLOT_SUFFIX), bbox_inches="tight")