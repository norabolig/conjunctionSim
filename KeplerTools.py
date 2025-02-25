#
# By norabolig  (A Boley)
#
# cc-by-sa 4.0
#
# github up soon
#
#

import numpy as np

twopi=np.pi*2.0

def vrad(f,n,a,ecc): return n*a*ecc*np.sin(f)/np.sqrt(1-ecc**2)

def vaz(f,n,a,ecc): return n*a*(1+ecc*np.cos(f))/np.sqrt(1-ecc**2)

def radial(f,a,e): return a*(1-e**2)/(1+e*np.cos(f))

def kepResid(ecc,EA,MA):
  return EA - ecc * np.sin(EA) - MA

def kepResidHyper(ecc,EA,MA):
  return -EA + ecc * np.sinh(EA) - MA

def kepEq(MA,ecc,EA=0,tol=1e-4):
  MA = MA/(twopi)
  MA = (MA-np.floor(MA))*twopi      # ensures between 0 and 2 pi
  A = kepResid(ecc,EA,MA)
  iter=0
  EA=MA*1
  if ecc>0.8: EA=np.pi*1
  while np.abs(A)>tol:
      EA = MA + ecc*np.sin(EA)
      A = kepResid(ecc,EA,MA)
      iter+=1
      if iter > 1000:
         print("Convergence issue for MA {} ecc {} EA {} and iter {}".format(MA,ecc,EA,iter))
         iter=0
  return EA


def kepEqHyper(MA,ecc,EA=0,tol=1e-6):
  A = kepResidHyper(ecc,EA,MA)
  while np.abs(A)>tol:
      dAdE = -1 + ecc * np.cosh(EA)
      EA = EA - A/dAdE
      A = kepResidHyper(ecc,EA,MA)
  return EA

def EAnom(ecc,TA): 
  f=2*np.arctan(np.sqrt( (1-ecc)/(1+ecc) ) * np.tan(TA/2.) )
  if f<0: f+=twopi
  return f


def trueAnom(ecc,EA): 
  f=2*np.arctan(np.sqrt( (1+ecc)/(1-ecc) ) * np.tan(EA/2.) )
  if f<0: f+=twopi
  return f

def trueAnomHyper(ecc,EA): 
  f=2*np.arctan(np.sqrt( (ecc+1)/(ecc-1) ) * np.tanh(EA/2.) )
  if f<0: f+=twopi
  return f

def get_Qs(w,O,inc):
    Q=np.zeros((3,2))
    Q[0][0] =np.cos(w)*np.cos(O)-np.sin(w)*np.sin(O)*np.cos(inc)
    Q[0][1] =-np.sin(w)*np.cos(O)-np.cos(w)*np.sin(O)*np.cos(inc)
    Q[1][0] =np.cos(w)*np.sin(O)+np.sin(w)*np.cos(O)*np.cos(inc) 
    Q[1][1] = -np.sin(w)*np.sin(O)+np.cos(w)*np.cos(O)*np.cos(inc)
    Q[2][0] = np.sin(w)*np.sin(inc)
    Q[2][1] = np.cos(w)*np.sin(inc)
    return Q

def get_dQs(w0,O,inc):
    Q=np.zeros((3,2))
    Q[1][0] =-np.sin(w0)*np.cos(O)*np.sin(inc) 
    Q[1][1] = -np.cos(w0)*np.cos(O)*np.sin(inc)
    Q[2][0] = np.sin(w0)*np.cos(inc)
    Q[2][1] = np.cos(w0)*np.cos(inc)
    return Q

def get_dQOs(w0,O,inc):
    Q=np.zeros((3,2))
    Q[0][0]=-np.cos(w0)*np.sin(O)-np.sin(w0)*np.cos(O)*np.cos(inc)
    Q[0][1]=+np.sin(w0)*np.sin(O)-np.cos(w0)*np.cos(O)*np.cos(inc)
    Q[1][0] =np.cos(w0)*np.cos(O)-np.sin(w0)*np.sin(O)*np.cos(inc) 
    Q[1][1] = -np.sin(w0)*np.cos(O)-np.cos(w0)*np.sin(O)*np.cos(inc)
    Q[2][0] = 0.
    Q[2][1] = 0.
    return Q

def getXYZVVV(nu,a,w0,ecc,O,inc,m0,m1,G=6.6743e-11):
    n = np.sqrt(G*(m0+m1)/a**3)
    r = radial(nu,a,ecc)
    Q = get_Qs(w0,O,inc)
    X=r*np.cos(nu)*Q[0][0]+r*np.sin(nu)*Q[0][1]
    Y=r*np.cos(nu)*Q[1][0]+r*np.sin(nu)*Q[1][1]
    Z=r*np.cos(nu)*Q[2][0]+r*np.sin(nu)*Q[2][1]

    vr = vrad(nu,n,a,ecc)
    vf = vaz(nu,n,a,ecc)
    VXss=vr*np.cos(nu)-vf*np.sin(nu)
    VYss=vr*np.sin(nu)+vf*np.cos(nu)

    VX=VXss*Q[0][0]+VYss*Q[0][1]
    VY=VXss*Q[1][0]+VYss*Q[1][1]
    VZ=VXss*Q[2][0]+VYss*Q[2][1]

    return X,Y,Z,VX,VY,VZ

def dotProduct(x,y):
    a=0.
    for i in range(3): a+=x[i]*y[i]
    return a

def getORBELM(r_vec,v_vec,mu):

    r = np.sqrt(dotProduct(r_vec,r_vec))
    v = np.sqrt(dotProduct(v_vec,v_vec))
    vr = dotProduct(v_vec,r_vec)/r

    vperp = np.sqrt(v**2-vr**2)

    h_vec = np.cross(r_vec,v_vec)
    h = np.sqrt(h_vec[0]**2 + h_vec[1]**2 + h_vec[2]**2)

    inc = np.arccos(h_vec[2]/h)

    N_vec = np.cross(np.array([0,0,1]),h_vec)
    N = np.sqrt(N_vec[0]**2 + N_vec[1]**2 + N_vec[2]**2)
    Omega = np.arccos(N_vec[0]/N)
    if N_vec[1]<0: Omega=2*np.pi-Omega

    ecc_vec = np.cross(v_vec,h_vec)/mu - r_vec/r
    ecc = np.sqrt(ecc_vec[0]**2 + ecc_vec[1]**2 + ecc_vec[2]**2)

    product = 0.
    for i in range(3): product += ecc_vec[i]*N_vec[i]

    omega = np.arccos(dotProduct(ecc_vec,N_vec)/(ecc*N))
    if ecc_vec[2] < 0: omega=2*np.pi-omega

    TA = np.arccos(dotProduct(ecc_vec,r_vec)/(ecc*r))

    if vr < 0: TA = 2*np.pi-TA

    a = 1./(2./r - v*v/mu)

    return a,ecc,omega,inc,Omega,TA





