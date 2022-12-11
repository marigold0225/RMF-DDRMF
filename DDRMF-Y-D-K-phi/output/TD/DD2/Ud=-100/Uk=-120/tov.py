
from matplotlib import pyplot as plt
import sys
import numpy as np 
import pandas as pd
from scipy import interpolate
import math
# EoS/input file
# Need a header with 'pressure' and 'energy'
# for pressure and energy density columns, respectively.
# The column separator should be a space
# (you can change it in the "sep" option)
df = pd.read_csv('EOS.dat', sep='  ',engine='python')

# print header in the output file
#with open('output_tov.dat', 'a') as f:
#  sys.stdout = f # change standard output to file
#  print('cen_P Ridus Mass Tid k_two yR 2GM/R')

# 'for' with 'step' float values
def range_with_floats(start, stop, step):
  while stop > start:
    yield start
    start += step

# dP/dr TOV equation
def TOV():
  dPdr = -(Rsch/msol)* \
    ((eps + P)*(msol*m + 4*np.pi*(r/fmtokm)**3*P)) \
    /(r*(r - 2*Rsch*m))
  return dPdr

# k2 Love number
def Love():
  k2 =  (8*(C**5)/5)* \
    ((1-2*C)**2)*(2+2*C*(y-1)-y) \
    /(2*C*(6-3*y+3*C*(5*y-8))+4*C*C*C*(13-11*y+C*(3*y-2)+2*C*C*(1+y)) \
      +3*((1-2*C)**2)*(2-y+2*C*(y-1))*math.log(1-2*C))
  return k2

# d(beta)/dr equation
def Tidal():
  e2lambda = 1/(1-2*Rsch*m/r)
  dbetadr = 2*e2lambda*H* \
    (-2*np.pi*(Rsch/(msol*(fmtokm**3)))*(5*eps +9*P +(eps +P)/dPdeps) \
      +3/(r**2) +2*e2lambda*(Rsch*m/(r**2)
      +4*np.pi*(fmtokm**(-3))*(Rsch/msol)*r*P)**2) \
    +2*(beta/r)*e2lambda*(-1 +Rsch*m/r \
      +2*np.pi*(fmtokm**(-3))*(Rsch/msol)*(r**2)*(eps -P))
  return dbetadr

# unit transformations
fmtokm = 1.0e-18 # fm = (fmtokm) km
kmtog = 3.5177e-43 # km^(-1) = (kmtog) g
conv = 5.07e-3 # Mev = (conv) fm^(-1)
conv2 = 1.7825e27 # MeV/fm^3 = (conv2) g/km^3

# constants
msol = (2.0e33*fmtokm)/(conv*kmtog) # solar mass (2.0e33 g -> MeV)
Rsch = 1.476 # Rsch = G*msol (km)

# radius step accuracy
delta_r = 0.001 # km

# surface pressure (usually zero)
P_sup = 0 # MeV/fm^3

# Euler Method
for j in range_with_floats(-26,10,0.25):
# you might need to change the j values depending on your EoS

  # central values/boundary conditions
  P0 = 10**(-j/10) # pressure (MeV/fm^3)
  m0 = 0 # mass (M_sun)
  r0 = 0.001 # radius (km)
  H0 = r0*r0 # eq.16 of arXiv:0711.2420 (km^2)
  beta0 = 2*r0 # beta(r)=H'(r) (km) - arXiv:0711.2420

  # initial values (iteration #0)
  P = P0 # pressure (MeV/fm^3)
  m = m0 # mass (M_sun)
  r = r0 # radius (km)
  H = H0 # linearized perturbation metric (km^2) - arXiv:0711.2420
  beta = beta0 # beta(r) = H'(r) (km) - arXiv:0711.2420
  # Spline Interpolation: energy = energy(pressure)
  tck = interpolate.splrep(df.pressure, df.energy, s=0) # interp. EoS file
  eps = interpolate.splev(P,tck) # energy (MeV/fm^3)
  # Differential equations (TOV)
  dMdr = 4.0*np.pi*(r**2)*eps/(msol*(fmtokm**3))
  dPdr = TOV()

  # iteration #j
  for i in range(1,100000):

    # saves the values from the previous interaction.
    # required to calculate dP/d(eps) for tidal deformability
    r_before = r
    P_before = P
    eps_before = eps
    m_before = m

    if(P <= P_sup): # reached star's surface
      y = r*beta/H \
        -4*np.pi*(r**3)*eps/(msol*m*(fmtokm**3)) # eq.17 arXiv:2010.07448
      C = Rsch*m/r
      k2 = Love() # eq.23 of arXiv:0711.2420
      Lambda = (2/3)*k2/(C**5) # tidal deformability

      # print the result in the output file
      with open('output_tov.dat', 'a') as f:
        sys.stdout = f # change standard output to file
        print('%.3f %.3f %.3f %.3e %.3e %.3e %.3e' \
          % (P0,r,m,Lambda,k2,y,C))
      break # end iteration #j
    else: # obtain values at r = r + delta_r
      P += dPdr*delta_r
      r += delta_r
      m += dMdr*delta_r
      eps=interpolate.splev(P,tck)

      dMdr = 4.0*np.pi*(r**2)*eps/(msol*(fmtokm**3))
      dPdr = TOV()

      dPdeps = (P - P_before)/(eps - eps_before) # dP/d(eps)
      dbetadr = Tidal()
      dHdr = beta
      beta += dbetadr*delta_r
      H += dHdr*delta_r
