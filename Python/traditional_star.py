# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 10:38:19 2016

@author: peterweilnboeck
"""
import math as math
import matplotlib.pyplot as plt
#import numpy as np
from numpy import genfromtxt
from numpy import transpose
#from scipy.interpolate import interp1d
from scipy import interpolate


global ggrav
ggrav = 6.673*10**-8
global clite
clite = 2.99792458*10**10
global g_c2
g_c2 = ggrav/clite**2
global msun
msun = 1.98892*10**33
global pi
pi = math.pi
#global gamma
#gamma = 2.75*10**0
#global eosk
#eosk = 1.98183*10**-6
global min_press
min_press = 1*10**-9
global delta
delta = 0.23
rho_center = 2.0*(10**14)


my_data = genfromtxt('EOS-numbers1.csv', delimiter=',')
my_data1 = transpose(my_data)
global rho_EOS
rho_EOS = my_data1[0]
global p_EOS
p_EOS = my_data1[1]

#rhoFromP = interp1d(p_EOS, rho_EOS)
#pFromRho = interp1d(rho_EOS, p_EOS)
#drho_dp = np.diff(rhoFromP) / np.diff(p_EOS)

rhoFromP = interpolate.InterpolatedUnivariateSpline(p_EOS, rho_EOS, k=1)
pFromRho = interpolate.InterpolatedUnivariateSpline(rho_EOS, p_EOS, k=1)
dRhoDP = rhoFromP.derivative()

def make_a_tov_star(rho_C):
    rad, mass = tov_integrate(rho_C)
    return rad, mass

def tov_integrate( roh_c ):
    press = pFromRho(roh_c) 
    press_now = press
    m_now = 0.0
    count = 0
    m_now = 0.0
    rad = 0.0
    dr = 1.0*10**-20
#    print "P=",press_now
#    print "M=",m_now
#    print "R=",rad
    while press_now >= min_press and count < 1000000:
        count+=1
        press_now, m_now = tov_rk4(press_now, m_now, rad, dr)
        rad = rad + dr
        dpress, dm = tov_RHS(press_now, m_now, rad)
        dr = delta*1/((1/m_now)*dm - (1/press_now)*dpress)
#        if count%100 == 0:
#            print "P=",press_now
#            print "M=",m_now/msun
#            print "R=",rad*(10**-5)
#            print "Count = ",count
#    print "P=",press_now
#    print "M=",m_now/msun
#    print "R=",rad*(10**-5)
#    print "Count = ",count
    return rad*(10**-5),m_now/msun
        
        
def tov_rk4(press_old, m_old, rad, dr):    
    press_temp = press_old
    m_temp = m_old
    press_k1, m_k1 = tov_RHS(press_temp, m_temp, rad)
    press_k1 = press_k1*dr
    m_k1 = m_k1*dr
    
    press_temp = press_old + 0.5*press_k1
    m_temp = m_old + 0.5*m_k1
    press_k2, m_k2 = tov_RHS(press_temp, m_temp, rad + 0.5*dr)
    press_k2 = press_k2*dr
    m_k2 = m_k2*dr
    
    press_temp = press_old + 0.5*press_k2
    m_temp = m_old + 0.5*m_k2
    press_k3, m_k3 = tov_RHS(press_temp, m_temp, rad + 0.5*dr)
    press_k3 = press_k3*dr
    m_k3 = m_k3*dr
    
    press_temp = press_old + press_k3
    m_temp = m_old + m_k3
    press_k4, m_k4 = tov_RHS(press_temp, m_temp, rad + dr)
    press_k4 = press_k4*dr
    m_k4 = m_k4*dr
    
    press_new = press_old + (press_k1 + press_k2 + press_k2 + press_k3 + press_k3 + press_k4) / 6
    m_new = m_old + (m_k1 + m_k2 + m_k2 + m_k3 + m_k3 + m_k4) / 6
    
    return press_new, m_new
    

def tov_RHS(press, mgrav, r):
    press = max(press, min_press)
    rho = rhoFromP(press)
    if(r <= 1*10**-6): 
        beast=(4.0*pi*r*press/clite**2)
    else:
        beast=(mgrav+4.0*pi*r*r*r*press/clite**2)/(r*(r-2.0*ggrav*mgrav/clite**2))

    dpress = -ggrav*( rho+press/clite**2 )*beast
    dm = 4.0*pi*r*r*rho
    return dpress, dm
 
    
rho_center = 1*(10**6)
rhocs = []
radii = []
masses = []
target = open("SimData.txt", 'w')
line = 'rho R(km) M(SM) \n'
target.write(line)
while rho_center <= 5.11*(10**16):
    radius, mass = make_a_tov_star(rho_center)
    print "rho_c=",rho_center," Radius=",radius,"Mass=",mass
    rhocs.append(rho_center)
    radii.append(radius)
    masses.append(mass)
    line = "%f %f %f \n"%(rho_center, radius, mass)
    target.write(line)
#    rho_center = min(rho_center * 1.005, rho_center+1*10**13)
    rho_center = rho_center*1.1
   
target.close()    
plt.plot(radii,masses)
plt.xscale('log')
plt.show