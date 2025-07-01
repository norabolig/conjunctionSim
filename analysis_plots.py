import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm

MAXSIMTIME = 60*60*24*2  # full sim run time
CLOSE=1e3 # close encounter distance 
DTMAX=0.05
TMIN=0  # min time for plot
TMAX=MAXSIMTIME  # max time for plot
DTSEP=1000
PLOTSTAR=True
PLOT_SUFFIX='_random_starlink_2d_15km'  # plot name suffix
PLOT_STARLINK = True # choose if you want to note number of starlink-specific encounters
SAVE = True  # choose if you want to save plots

OUTPUTFILE = 'out/Starlink_random/encounter_data_random_starlink_2d_15km.hdf'  # FOR IF YOU WANT TO SAVE ENCOUNTER DATA

dat = pd.read_hdf(OUTPUTFILE)

tclose=dat.time.values.astype(float)
altclose=dat.alt.values.astype(float)
dclose=dat.dist.values.astype(float)
vclose=dat.dV.values.astype(float)
cclose=dat.type.values.astype(str)

subset = dat.loc[(dclose>0)&(dclose<CLOSE)&(tclose>TMIN)&(tclose<TMAX)]

#========================================================================================

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
    plt.savefig("./out/Starlink_random/relative_dist_vel{}.pdf".format(PLOT_SUFFIX), bbox_inches="tight")

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
    plt.savefig("./out/Starlink_random/vclose_sat_hist{}.pdf".format(PLOT_SUFFIX), bbox_inches="tight")

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
    plt.savefig("./out/Starlink_random/time_alt_encounters{}.pdf".format(PLOT_SUFFIX), bbox_inches="tight")

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
    plt.savefig("./out/Starlink_random/time_encounters_hist{}.pdf".format(PLOT_SUFFIX), bbox_inches="tight")

#========================================================================================

fig, ax = plt.subplots(figsize=(8,6))
altbins = np.arange(200,2000,1)
h,bins = np.histogram(altclose/1e3,bins=altbins)
binmids = (bins[:-1]+bins[1:])/2
plt.plot(binmids,MAXSIMTIME/(h*60),color='k')
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
    plt.savefig("./out/Starlink_random/frequency{}.pdf".format(PLOT_SUFFIX), bbox_inches="tight")
