#
# By norabolig  (A Boley)
#
# cc-by-sa 4.0
#
# github up soon
#
#

import numpy as np
import numba as nb

from numba import njit

twopi=np.pi*2.0
MEarth = 5.97e24

@njit
def vrad(f,n,a,ecc): return n*a*ecc*np.sin(f)/np.sqrt(1-ecc**2)

@njit
def vaz(f,n,a,ecc): return n*a*(1+ecc*np.cos(f))/np.sqrt(1-ecc**2)

@njit
def radial(f,a,e): return a*(1-e**2)/(1+e*np.cos(f))

@njit
def kepResid(ecc,EA,MA):
  return EA - ecc * np.sin(EA) - MA

@njit
def kepResidHyper(ecc,EA,MA):
  return -EA + ecc * np.sinh(EA) - MA

@njit
def kepEq(MA,ecc,EA=0,tol=1e-4):
  MA = MA/(twopi)
  MA = (MA-np.floor(MA))*twopi      # ensures between 0 and 2 pi
  A = kepResid(ecc,EA,MA)
  EA=MA

  if ecc > 0.8: EA = np.pi

  i = 0
  while np.abs(A) > tol and i <= 1000:
      EA = MA + ecc*np.sin(EA)
      A = kepResid(ecc,EA,MA)
      i += 1

  if i > 1000:
      print("Convergence issue for MA {} ecc {} EA {} and iter {}".format(MA,ecc,EA,iter))
      i = 0

  return EA

@njit
def kepEqHyper(MA,ecc,EA=0,tol=1e-6):
  A = kepResidHyper(ecc,EA,MA)
  
  while np.abs(A) > tol:
      dAdE = -1 + ecc * np.cosh(EA)
      EA = EA - A/dAdE
      A = kepResidHyper(ecc,EA,MA)

  return EA

@njit
def EAnom(ecc,TA): 
  f=2*np.arctan(np.sqrt( (1-ecc)/(1+ecc) ) * np.tan(TA/2.) )
  if f<0: f+=twopi
  return f

@njit
def trueAnom(ecc,EA): 
  f=2*np.arctan(np.sqrt( (1+ecc)/(1-ecc) ) * np.tan(EA/2.) )
  if f<0: f+=twopi
  return f

@njit
def meanAnom(ecc,TA):
   EA = EAnom(ecc, TA)
   return EA - ecc*np.sin(EA)

@njit
def trueAnomHyper(ecc,EA): 
  f=2*np.arctan(np.sqrt( (ecc+1)/(ecc-1) ) * np.tanh(EA/2.) )
  if f<0: f+=twopi
  return f

@njit(nb.float64[:, :, :](nb.float64[:], nb.float64[:], nb.float64[:]))
def get_Qs(w, O, inc) -> np.ndarray:
    
    Qs = np.zeros((len(w), 3, 2))

    cosW, sinW = np.cos(w), np.sin(w)
    cosO, sinO = np.cos(O), np.sin(O)
    cosI, sinI = np.cos(inc), np.sin(inc)

    Qs[:, 0, 0] = cosW*cosO - sinW*sinO*cosI
    Qs[:, 0, 1] =-sinW*cosO - cosW*sinO*cosI
    Qs[:, 1, 0] = cosW*sinO + sinW*cosO*cosI 
    Qs[:, 1, 1] =-sinW*sinO + cosW*cosO*cosI
    Qs[:, 2, 0] = sinW*sinI
    Qs[:, 2, 1] = cosW*sinI

    # print(Qs)

    return Qs

# @njit
# def get_dQs(w0,O,inc):
#     Q=np.zeros((3,2))
#     Q[1][0] =-np.sin(w0)*np.cos(O)*np.sin(inc) 
#     Q[1][1] = -np.cos(w0)*np.cos(O)*np.sin(inc)
#     Q[2][0] = np.sin(w0)*np.cos(inc)
#     Q[2][1] = np.cos(w0)*np.cos(inc)
#     return Q

# @njit
# def get_dQOs(w0,O,inc):
#     Q=np.zeros((3,2))
#     Q[0][0]=-np.cos(w0)*np.sin(O)-np.sin(w0)*np.cos(O)*np.cos(inc)
#     Q[0][1]=+np.sin(w0)*np.sin(O)-np.cos(w0)*np.cos(O)*np.cos(inc)
#     Q[1][0] =np.cos(w0)*np.cos(O)-np.sin(w0)*np.sin(O)*np.cos(inc) 
#     Q[1][1] = -np.sin(w0)*np.cos(O)-np.cos(w0)*np.sin(O)*np.cos(inc)
#     Q[2][0] = 0.
#     Q[2][1] = 0.
#     return Q

@njit
def getXYZVVV(nu, a, w0, ecc, O, inc, n, m0=MEarth, m1=0., G=6.6743e-11):
    r = radial(nu,a,ecc)
    Q = get_Qs(w0,O,inc)

    cosNu, sinNu = np.cos(nu), np.sin(nu)

    pos = np.column_stack((r*cosNu*Q[:, 0, 0]+r*sinNu*Q[:, 0, 1], 
                    r*cosNu*Q[:, 1, 0]+r*sinNu*Q[:, 1, 1], 
                    r*cosNu*Q[:, 2, 0]+r*sinNu*Q[:, 2, 1]))

    vr = vrad(nu,n,a,ecc)
    vf = vaz(nu,n,a,ecc)

    VXss=vr*cosNu-vf*sinNu
    VYss=vr*sinNu+vf*cosNu

    vel = np.column_stack((VXss*Q[:, 0, 0]+VYss*Q[:, 0, 1],
                    VXss*Q[:, 1, 0]+VYss*Q[:, 1, 1],
                    VXss*Q[:, 2, 0]+VYss*Q[:, 2, 1]))

    return pos, vel

@njit
def dotProduct(x,y): return np.dot(x, y)

@njit
def getORBELM(r_vec,v_vec,mu):

    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)
    vr = np.dot(v_vec, r_vec)/r

    vperp = np.sqrt(v**2 - vr**2)

    h_vec = np.cross(r_vec,v_vec)
    h = np.linalg.norm(h_vec)

    inc = np.arccos(h_vec[2]/h)

    N_vec = np.cross(np.array([0,0,1]),h_vec)
    N = np.linalg.norm(N_vec)

    Omega = np.arccos(N_vec[0]/N)
    if N_vec[1] < 0: Omega = 2*np.pi - Omega

    ecc_vec = np.cross(v_vec,h_vec)/mu - r_vec/r
    ecc = np.linalg.norm(ecc_vec)

    omega = np.arccos(np.dot(ecc_vec,N_vec)/(ecc*N))
    if ecc_vec[2] < 0: omega = 2*np.pi - omega

    TA = np.arccos(np.dot(ecc_vec,r_vec)/(ecc*r))

    if vr < 0: TA = 2*np.pi - TA

    a = 1./(2./r - v*v/mu)

    return a, ecc, omega, inc, Omega, TA

@njit
def getNorm(arr: np.ndarray) -> np.ndarray:
  norms = np.zeros(len(arr))
  for i in range(len(arr)):
    norms[i] = np.sqrt(arr[i, 0]*arr[i, 0] + arr[i, 1]*arr[i, 1] + arr[i, 2]*arr[i, 2])
  return norms





