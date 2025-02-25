import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import curve_fit

fh=open("track_out.dat","r")

CLOSE=10000
DTMAX=0.05
sum=0

time=[]
dist=[]
iid=[]
jid=[]
v=[]
r=[]
iname=[]
jname=[]
for line in fh:
    tt,dt,vt,rt,it,jt,inme,jnme=line.split(",")
    time.append(float(tt))
    dist.append(float(dt))
    iid.append(int(it))
    jid.append(int(jt))
    v.append(float(vt))
    r.append(float(rt))
    iname.append(inme)
    jname.append(jnme)

fh.close()

time=np.array(time)
dist=np.array(dist)
iid=np.array(iid)
jid=np.array(jid)
v=np.array(v)
r=np.array(r)
iname=np.array(iname)
jname=np.array(jname)

dclose=[]
vclose=[]

plt.figure()
while True:
    if time.size<1: break
    id=iid[0]
    jd=jid[0]
    flag1 = iid==id
    flag2 = jid==jd
    flag = flag1*flag2
    ids=np.nonzero(flag)
    #print(ids)

    if min(dist[flag])<CLOSE:
         sum+=1 
         min_d = min(dist[flag])
         id = np.argmin(dist[flag])
         dclose.append(min_d)
         vclose.append(v[id])

    if iname[id].find("STARLINK")>-1:
        if jname[id].find("STARLINK")>-1:
           carray="red"
        else:
           carray="blue"
    elif jname[id].find("STARLINK")>-1:
        carray="blue"
    else:
        carray="green"

    plt.plot(time[flag]/60,dist[flag],color=carray)
    #if popt is not None: plt.scatter(time[flag],dopt,s=0.1)
    time=np.delete(time,flag)
    dist=np.delete(dist,flag)
    iid=np.delete(iid,flag)
    jid=np.delete(jid,flag)

plt.xlabel("Time in minutes")
plt.ylabel("Conjunction Distance [m]")
plt.savefig("close_approach_tracks.pdf")

print(vclose)
print(dclose)

plt.figure()
plt.scatter(dclose,vclose)

print("Total detected approache < {} m is {}".format(CLOSE,sum))

plt.show()


