# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 10:38:19 2016

@author: peterweilnboeck
"""
import math as math
#import cmath as math
import matplotlib.pyplot as plt
from numpy import genfromtxt
from numpy import transpose
from scipy import interpolate

global ggrav
ggrav = 6.673*10**-8
G = ggrav
global clite
clite = 2.99792458*10**10
global g_c2
g_c2 = ggrav/clite**2
global msun
msun = 1.98892*10**33
global pi
pi = math.pi
Pi = math.pi
#global gamma
#gamma = 2.75*10**0
#global eosk
#eosk = 1.98183*10**-6
global min_press
min_press = 1*10**-9
global delta
delta = 0.23
rho_center = 2.0*(10**14)
alphaQ = 1.0
Rp = (1.0/1.6)**33
a = 1/2


my_data = genfromtxt('EOS-numbers1.csv', delimiter=',')
my_data1 = transpose(my_data)
global rho_EOS
rho_EOS = my_data1[0]
global p_EOS
p_EOS = my_data1[1]

rhoFromP = interpolate.InterpolatedUnivariateSpline(p_EOS, rho_EOS, k=1)
pFromRho = interpolate.InterpolatedUnivariateSpline(rho_EOS, p_EOS, k=1)
dRhoDP = rhoFromP.derivative()
#drho_dp = np.diff(rhoFromP) / np.diff(p_EOS)

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
    

def tov_RHS(P, M, r):
    P = max(P, min_press)*g_c2
    rho = rhoFromP(P)*g_c2/(clite**2)
    M=M*g_c2
    
    fR = 1 + (16*a*G*Pi*(-3*P + rho))/Rp
    f = (8*G*Pi*(-3*P + rho) + (64*a*G**2*Pi**2*(-3*P + rho)**2)/Rp - 
        (fR/(4.*alphaQ) + 16*G*P*Pi + 8*G*Pi*(-3*P + rho) - 
        (alphaQ*(3*(fR/alphaQ + 8*G*Pi*(-3*P + rho)) + 
        math.sqrt((-32*G*Pi*(P + rho))/alphaQ + (fR/alphaQ + 
        8*G*Pi*(-3*P + rho))**2))**2)/16. + (64*a*G**2*Pi**2*(-3*P + rho)**2)/Rp)/Rp)
    lambda1 = ((math.sqrt(alphaQ)*(3*(fR/alphaQ + 8*G*Pi*(-3*P + rho)) + 
        math.sqrt((-32*G*Pi*(P + rho))/alphaQ + 
        (fR/alphaQ + 8*G*Pi*(-3*P + rho))**2)))/(4.*math.sqrt(2)))  
    sigma1 = (fR/2.0 + math.sqrt(2)*math.sqrt(alphaQ)*
        math.sqrt(-8*G*Pi*(P + rho) + lambda1**2))    
    sigma2 = fR/2. + math.sqrt(2)*math.sqrt(alphaQ)*lambda1
#    print lambda1
#    print fR
#    print sigma1
#    print sigma2
#    print "f=",f
    detSigma = 1/math.sqrt(sigma1 * sigma2 * sigma2 * sigma2)
#    print "detSigma=",detSigma
    tau_rr = detSigma*((f/2 + 8*Pi*G*P))
    tau_tt = detSigma*((f/2 + 8*Pi*G*P)-8*Pi*G*(rho+P))
    Omega = (math.sqrt(-(((240*a*G*P*Pi - 80*a*G*Pi*rho - 5*Rp + 72*alphaQ*G*P*Pi*Rp - 
        24*alphaQ*G*Pi*rho*Rp - alphaQ*Rp*math.sqrt((-32*alphaQ*G*Pi*(P + rho) + 
        (1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))*
        (-96*a*G*P*Pi + 32*a*G*Pi*rho + 2*Rp + 4*math.sqrt(2)*math.sqrt(alphaQ)*Rp*
        math.sqrt(-8*G*Pi*(P + rho) + (alphaQ*((3*(1 - (8*G*Pi*(3*P - rho)*
        (2*a + alphaQ*Rp))/Rp))/alphaQ + math.sqrt((-32*alphaQ*G*Pi*(P + rho) + 
        (1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))**2)/32.)))/Rp**2))/4.)
    dOmagaPDOmega = ((-8*Rp**2*(((-240*a*G*P*Pi + 80*a*G*Pi*rho + 5*Rp - 
        72*alphaQ*G*P*Pi*Rp + 24*alphaQ*G*Pi*rho*Rp + alphaQ*Rp*
        math.sqrt((-32*alphaQ*G*Pi*(P + rho) + (1 - (8*G*Pi*(3*P - rho)*
        (2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))*((-24*a*G*Pi)/Rp + 
        (math.sqrt(alphaQ)*G*Pi*(-16 + ((-9*alphaQ*Rp*(2*a + alphaQ*Rp) + 
        (96*a**2*G*Pi*(3*P - rho) + 6*a*(-1 + 16*alphaQ*G*Pi*(3*P - rho))*Rp + 
        alphaQ*(-5 + 24*alphaQ*G*Pi*(3*P - rho))*Rp**2)/math.sqrt((-32*alphaQ*G*Pi*(P + rho) + 
        (1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))*((3*(1 - 
        (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp))/alphaQ + 
        math.sqrt((-32*alphaQ*G*Pi*(P + rho) + (1 - (8*G*Pi*(3*P - rho)*
        (2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2)))/(alphaQ*Rp**2)))/(2.*math.sqrt(2)*
        math.sqrt(-8*G*Pi*(P + rho) + (alphaQ*((3*(1 - (8*G*Pi*(3*P - rho)*(2*a + 
        alphaQ*Rp))/Rp))/alphaQ + math.sqrt((-32*alphaQ*G*Pi*(P + rho) + 
        (1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))**2)/
        32.))))/(4.*Rp) + (G*Pi*(-12*a*Rp + (-9*alphaQ*Rp*(2*a + alphaQ*Rp) + 
        (96*a**2*G*Pi*(3*P - rho) + 6*a*(-1 + 16*alphaQ*G*Pi*(3*P - rho))*Rp + 
        alphaQ*(-5 + 24*alphaQ*G*Pi*(3*P - rho))*Rp**2)/math.sqrt((-32*
        alphaQ*G*Pi*(P + rho) + (1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/
        Rp)**2)/alphaQ**2))/alphaQ)*(-96*a*G*P*Pi + 32*a*G*Pi*rho + 2*Rp + 
        4*math.sqrt(2)*math.sqrt(alphaQ)*Rp*math.sqrt(-8*G*Pi*(P + rho) + 
        (alphaQ*((3*(1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp))/alphaQ + 
        math.sqrt((-32*alphaQ*G*Pi*(P + rho) + (1 - (8*G*Pi*(3*P - rho)*(2*a + 
        alphaQ*Rp))/Rp)**2)/alphaQ**2))**2)/32.)))/(2.*Rp**3)))/((240*a*G*P*Pi - 
        80*a*G*Pi*rho - 5*Rp + 72*alphaQ*G*P*Pi*Rp - 24*alphaQ*G*Pi*rho*Rp - 
        alphaQ*Rp*math.sqrt((-32*alphaQ*G*Pi*(P + rho) + (1 - (8*G*Pi*(3*P - rho)*
        (2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))*(-96*a*G*P*Pi + 32*a*G*Pi*rho + 
        2*Rp + 4*math.sqrt(2)*math.sqrt(alphaQ)*Rp*math.sqrt(-8*G*Pi*(P + rho) + 
        (alphaQ*((3*(1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp))/alphaQ + 
        math.sqrt((-32*alphaQ*G*Pi*(P + rho) + (1 - (8*G*Pi*(3*P - rho)*(2*a + 
        alphaQ*Rp))/Rp)**2)/alphaQ**2))**2)/32.))))
    S = sigma2**2 / math.sqrt(sigma1*sigma2)
    dSPDS = (4*Rp**2*math.sqrt(-(((240*a*G*P*Pi - 80*a*G*Pi*rho - 5*Rp + 
        72*alphaQ*G*P*Pi*Rp -  24*alphaQ*G*Pi*rho*Rp - alphaQ*
        Rp*math.sqrt((-32*alphaQ*G*Pi*(P + rho) + (1 - (8*G*Pi*(3*P - 
        rho)*(2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))*(-96*a*G*P*Pi + 
        32*a*G*Pi*rho + 2*Rp + 4*math.sqrt(2)*math.sqrt(alphaQ)*Rp*
        math.sqrt(-8*G*Pi*(P + rho) + (alphaQ*((3*(1 - (8*G*Pi*(3*P - rho)*
        (2*a + alphaQ*Rp))/Rp))/alphaQ + math.sqrt((-32*alphaQ*G*Pi*(P + rho) + 
        (1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))**2)/32.)))/
        Rp**2))*((4*G*Pi*(-240*a*G*P*Pi + 80*a*G*Pi*rho + 5*Rp - 72*alphaQ*G*P*Pi*Rp + 
        24*alphaQ*G*Pi*rho*Rp + alphaQ*Rp*math.sqrt((-32*alphaQ*G*Pi*(P + rho) + 
        (1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))*
        (-12*a*Rp + (-9*alphaQ*Rp*(2*a + alphaQ*Rp) + (96*a**2*G*Pi*(3*P - rho) + 
        6*a*(-1 + 16*alphaQ*G*Pi*(3*P - rho))*Rp + alphaQ*(-5 + 24*alphaQ*G*Pi*
        (3*P - rho))*Rp**2)/math.sqrt((-32*alphaQ*G*Pi*(P + rho) + 
        (1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))/
        alphaQ))/(Rp**3*math.sqrt(-(((240*a*G*P*Pi - 80*a*G*Pi*rho - 5*Rp + 
        72*alphaQ*G*P*Pi*Rp - 24*alphaQ*G*Pi*rho*Rp - alphaQ*Rp*
        math.sqrt((-32*alphaQ*G*Pi*(P + rho) + (1 - (8*G*Pi*(3*P - rho)*
        (2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))*(-96*a*G*P*Pi + 32*a*G*Pi*rho + 
        2*Rp + 4*math.sqrt(2)*math.sqrt(alphaQ)*Rp*math.sqrt(-8*G*Pi*(P + rho) + 
        (alphaQ*((3*(1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp))/alphaQ + 
        math.sqrt((-32*alphaQ*G*Pi*(P + rho) + (1 - (8*G*Pi*(3*P - rho)*
        (2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))**2)/32.)))/Rp**2))) - 
        (2*(-240*a*G*P*Pi + 80*a*G*Pi*rho + 5*Rp - 72*alphaQ*G*P*Pi*Rp + 
        24*alphaQ*G*Pi*rho*Rp + alphaQ*Rp*math.sqrt((-32*alphaQ*G*Pi*(P + rho) + 
        (1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp)**2)/
        alphaQ**2))**2*(((-240*a*G*P*Pi + 80*a*G*Pi*rho + 5*Rp - 
        72*alphaQ*G*P*Pi*Rp + 24*alphaQ*G*Pi*rho*Rp + 
        alphaQ*Rp*math.sqrt((-32*alphaQ*G*Pi*(P + rho) + (1 - (8*G*Pi*(3*P - rho)*
        (2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))*((-24*a*G*Pi)/Rp + 
        (math.sqrt(alphaQ)*G*Pi*(-16 + ((-9*alphaQ*Rp*(2*a + alphaQ*Rp) + 
        (96*a**2*G*Pi*(3*P - rho) + 6*a*(-1 + 16*alphaQ*G*Pi*(3*P - rho))*Rp + 
        alphaQ*(-5 + 24*alphaQ*G*Pi*(3*P - rho))*Rp**2)/math.sqrt((-32*alphaQ*G*Pi*(P + rho) + 
        (1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))*((3*(1 - 
        (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp))/alphaQ + 
        math.sqrt((-32*alphaQ*G*Pi*(P + rho) + (1 - (8*G*Pi*(3*P - rho)*
        (2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2)))/(alphaQ*Rp**2)))/(2.*math.sqrt(2)*
        math.sqrt(-8*G*Pi*(P + rho) + (alphaQ*((3*(1 - (8*G*Pi*(3*P - rho)*
        (2*a + alphaQ*Rp))/Rp))/alphaQ + math.sqrt((-32*alphaQ*G*Pi*(P + rho) + 
        (1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp)**2)/
        alphaQ**2))**2)/32.))))/(4.*Rp) + (G*Pi*(-12*a*Rp + 
        (-9*alphaQ*Rp*(2*a + alphaQ*Rp) + (96*a**2*G*Pi*(3*P - rho) + 
        6*a*(-1 + 16*alphaQ*G*Pi*(3*P - rho))*Rp + alphaQ*(-5 + 24*alphaQ*G*Pi*
        (3*P - rho))*Rp**2)/math.sqrt((-32*alphaQ*G*Pi*(P + rho) + 
        (1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))/alphaQ)*
        (-96*a*G*P*Pi + 32*a*G*Pi*rho + 2*Rp + 4*math.sqrt(2)*math.sqrt(alphaQ)*
        Rp*math.sqrt(-8*G*Pi*(P + rho) + (alphaQ*((3*(1 - (8*G*Pi*(3*P - rho)*
        (2*a + alphaQ*Rp))/Rp))/alphaQ + math.sqrt((-32*alphaQ*G*Pi*(P + rho) + 
        (1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))**2)/32.)))/
        (2.*Rp**3)))/(Rp**2*(-(((240*a*G*P*Pi - 80*a*G*Pi*rho - 5*Rp + 
        72*alphaQ*G*P*Pi*Rp - 24*alphaQ*G*Pi*rho*Rp - alphaQ*Rp*
        math.sqrt((-32*alphaQ*G*Pi*(P + rho) + (1 - (8*G*Pi*(3*P - rho)*
        (2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))*(-96*a*G*P*Pi + 32*a*G*Pi*rho + 
        2*Rp + 4*math.sqrt(2)*math.sqrt(alphaQ)*Rp*math.sqrt(-8*G*Pi*(P + rho) + 
        (alphaQ*((3*(1 - (8*G*Pi*(3*P - rho)*(2*a + alphaQ*Rp))/Rp))/alphaQ + 
        math.sqrt((-32*alphaQ*G*Pi*(P + rho) + (1 - (8*G*Pi*(3*P - rho)*(2*a + 
        alphaQ*Rp))/Rp)**2)/alphaQ**2))**2)/32.)))/Rp**2))**1.5)))/(-240*a*G*P*Pi + 
        80*a*G*Pi*rho + 5*Rp - 72*alphaQ*G*P*Pi*Rp + 24*alphaQ*G*Pi*rho*Rp + 
        alphaQ*Rp*math.sqrt((-32*alphaQ*G*Pi*(P + rho) + (1 - (8*G*Pi*(3*P - rho)*
        (2*a + alphaQ*Rp))/Rp)**2)/alphaQ**2))**2
    beta = 2*r* dOmagaPDOmega * (1 - (rho+P)/2 * (3*dOmagaPDOmega/2-(dOmagaPDOmega - dSPDS)))
    alpha = (rho+P)/2 * (dOmagaPDOmega + dSPDS)
#    print "rho",rho
#    print "P",P
#    print "M",M
#    print "taurr",tau_rr
#    print "Omega",Omega
#    print "tautt",tau_tt
#    print "S",S
#    print "r",r
#    print (M - (tau_rr + Omega * tau_tt / S) * (r**3) / 4)
#    print (rho + P)*(M - (tau_rr + Omega * tau_tt / S) * (r**3) / 4)
    dP0 = (rho + P)*(M - (tau_rr + Omega * tau_tt / S) * (r**3) / 4)
    print "dP0=",dP0
    print "alpha=",alpha
    print "beta=",beta
    print "beta*dP0=",beta*dP0
    print dOmagaPDOmega*P+2/r
    print math.sqrt(1-beta*dP0)
    dP = 2 * dP0 /((1-alpha)*(1+math.sqrt(1-beta*dP0)))
    dm = 4.0*pi*r*r*rho
    return dP, dm
 
    
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
#    print "did it"
   
target.close()    
plt.plot(radii,masses)
plt.xscale('log')
plt.show