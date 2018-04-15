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
Rp = (1.0/1.6)**33
alphaQ = 1.0/Rp
a = 1/2
s=1

my_data = genfromtxt('EOS-numbers1.csv', delimiter=',')
my_data1 = transpose(my_data)
global rho_EOS
rho_EOS = my_data1[0]
global p_EOS
p_EOS = my_data1[1]

rhoFromP = interpolate.InterpolatedUnivariateSpline(p_EOS, rho_EOS, k=1)
pFromRho = interpolate.InterpolatedUnivariateSpline(rho_EOS, p_EOS, k=1)
rhoFromP2 = interpolate.InterpolatedUnivariateSpline(p_EOS, rho_EOS, k=2)
dRhoDP = rhoFromP.derivative()
ddRhoDP = rhoFromP2.derivative(2)
#drho_dp = np.diff(rhoFromP) / np.diff(p_EOS)

def make_a_tov_star(rho_C):
    rad, mass = tov_integrate(rho_C)
    return rad, mass

def tov_integrate( roh_c ):
    press = pFromRho(roh_c) 
    press_now = (press)
    m_now = (0.0)
    count = 0
    m_now = 0.0
    rad = 0.0
    dr = 1.0*10**-20
    dm = 0.0
#    print "P=",press_now
#    print "M=",m_now
#    print "R=",rad
    target = open("SingleData.txt", 'w')
    line = 'R(km) P M(SM) dP dM\n'
    target.write(line)
    while press_now >= min_press and count < 1000000:
        count+=1
        press_now, m_now = tov_rk4(press_now, m_now, rad, dr, dm)
        rad = rad + dr
        dpress, dm = tov_RHS(press_now, m_now, rad, dm)
        dr = delta*1/((1/m_now)*dm - (1/press_now)*dpress)
        line = "%.20f %f %f %.20f %.20f\n"%(rad, press_now, m_now,dpress,dm)
        target.write(line)
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
        
        
def tov_rk4(press_old, m_old, rad, dr, dm):    
    press_temp = press_old
    m_temp = m_old
    press_k1, m_k1 = tov_RHS(press_temp, m_temp, rad, dm)
    press_k1 = press_k1*dr
    m_k1 = m_k1*dr
    
    press_temp = press_old + 0.5*press_k1
    m_temp = m_old + 0.5*m_k1
    press_k2, m_k2 = tov_RHS(press_temp, m_temp, rad + 0.5*dr, dm)
    press_k2 = press_k2*dr
    m_k2 = m_k2*dr
    
    press_temp = press_old + 0.5*press_k2
    m_temp = m_old + 0.5*m_k2
    press_k3, m_k3 = tov_RHS(press_temp, m_temp, rad + 0.5*dr, dm)
    press_k3 = press_k3*dr
    m_k3 = m_k3*dr
    
    press_temp = press_old + press_k3
    m_temp = m_old + m_k3
    press_k4, m_k4 = tov_RHS(press_temp, m_temp, rad + dr, dm)
    press_k4 = press_k4*dr
    m_k4 = m_k4*dr
    
    press_new = press_old + (press_k1 + press_k2 + press_k2 + press_k3 + press_k3 + press_k4) / 6
    m_new = m_old + (m_k1 + m_k2 + m_k2 + m_k3 + m_k3 + m_k4) / 6
    
    press_change = (press_k1 + press_k2 + press_k2 + press_k3 + press_k3 + press_k4) / 6
    print "press_change=",press_change
    print "press_old=%.10f"%press_old
    print "press_new=%.10f"%press_new 
    print "press_sho=%.10f"%(press_old + press_change)
    
    return press_new, m_new
    

def tov_RHS(P, M, r, dM):
    P = max(P, min_press)*g_c2/(clite**2)
    rho = rhoFromP(P)*g_c2
    drhoP = dRhoDP(P)*(clite**2)
    ddrhoPP = ddRhoDP(P)*(clite**6)/G
    M=M*g_c2
    
    fR = 1 + (16*a*G*Pi*(-3*P + rho))/Rp
    f = (8*G*Pi*(-3*P + rho) + (64*a*G**2*Pi**2*(-3*P + rho)**2)/Rp - 
        (fR/(4.*alphaQ) + 16*G*P*Pi + 8*G*Pi*(-3*P + rho) - 
        (alphaQ*(3*(fR/alphaQ + 8*G*Pi*(-3*P + rho)) + 
        math.sqrt((-32*G*Pi*(P + rho))/alphaQ + (fR/alphaQ + 
        8*G*Pi*(-3*P + rho))**2))**2)/16. + (64*a*G**2*Pi**2*(-3*P + rho)**2)/Rp)/Rp)
        
    lambda1 = ((math.sqrt(alphaQ)*(3*(fR/alphaQ + 8*G*Pi*(-3*P + rho)) - 
        math.sqrt((-32*G*Pi*(P + rho))/alphaQ + 
        (fR/alphaQ + 8*G*Pi*(-3*P + rho))**2)))/(4.*math.sqrt(2))) 
        
    dlambdaP = (math.sqrt(2)*math.sqrt(alphaQ)*G*Pi*(3*(-3 + drhoP) - 
        ((-2*(1 + drhoP))/alphaQ - (-3 + drhoP)*(fR/alphaQ + 8*G*Pi*(-3*P + rho)))/
        math.sqrt((-32*G*Pi*(P + rho))/alphaQ + (fR/alphaQ + 8*G*Pi*(-3*P + rho))**2)))
        
    sigma1 = (fR/2.0 + math.sqrt(2)*math.sqrt(alphaQ)*
        math.sqrt(-8*G*Pi*(P + rho) + lambda1**2))    
    sigma2 = fR/2. + math.sqrt(2)*math.sqrt(alphaQ)*lambda1
    dsigma1P = ((math.sqrt(2)*math.sqrt(alphaQ)*(dlambdaP*lambda1 - 
        4*(1 + drhoP)*G*Pi))/math.sqrt(lambda1**2 - 8*G*Pi*(P + rho)))
    dsigma2P = math.sqrt(2)*math.sqrt(alphaQ)*dlambdaP  
#    print lambda1
#    print fR
#    print sigma1
#    print sigma2
#    print "f=",f
    detSigma = 1/math.sqrt(sigma1 * sigma2 * sigma2 * sigma2)
#    print "detSigma=",detSigma
    tau_rr = detSigma*((f/2 + 8*Pi*G*P))
    tau_tt = detSigma*((f/2 + 8*Pi*G*P)-8*Pi*G*(rho+P))  
    
    Omega = math.sqrt(sigma1*sigma2)
    dOmegaP = (sigma2*dsigma1P + sigma1*dsigma2P)/(2.*math.sqrt(sigma1*sigma2))
    dOmegaPDOmega = dOmegaP/Omega
    S = sigma2**2 / math.sqrt(sigma1*sigma2)
    dSP = (sigma2**2*(3*dsigma2P*sigma1 - dsigma1P*sigma2))/(2.*(sigma1*sigma2)**1.5)
    dSPDS = S/dSP
    beta = 2*r* dOmegaPDOmega * (1 - (rho+P)/2 * (3*dOmegaPDOmega/2-(dOmegaPDOmega - dSPDS)))
    alpha = (rho+P)/2 * (dOmegaPDOmega + dSPDS)
#    print "rho",rho
#    print "P",P
#    print "M",M
#    print "taurr",tau_rr
    print "Omega",Omega
#    print "tautt",tau_tt
    print "S",S
    print "fR",fR
    print "sigma1",sigma1
    print "sigma2",sigma2
    print "lambda",lambda1
    print "Sqrt(2 * alpha)*lambda",lambda1*math.sqrt(2*alphaQ)
    print "detSigma",1/detSigma
#    print "r",r
#    print (M - (tau_rr + Omega * tau_tt / S) * (r**3) / 4)
#    print (rho + P)*(M - (tau_rr + Omega * tau_tt / S) * (r**3) / 4)
    
    
    dP0 = (rho + P)*(M - (tau_rr + Omega * tau_tt / S) * (r**3) / 4)
#    print "dP0=",dP0
#    print "alpha=",alpha
#    print "beta=",beta
#    print "beta*dP0=",beta*dP0
#    print dOmagaPDOmega*P+2/r
#    print math.sqrt(1-beta*dP0)
    dP = 2 * dP0 /((1-alpha)*(1+math.sqrt(1-beta*dP0))*g_c2/(clite**2))
    print "dP=",dP
#---dP------------------------------------------------------------------------

    ddlambdaPP = (math.sqrt(2)*math.sqrt(alphaQ)*G*Pi*(3*ddrhoPP - 
        (8*G*Pi*((-2*(1 + drhoP))/alphaQ + (-3 + drhoP)*(fR/alphaQ + 
        8*G*Pi*(-3*P + rho)))**2)/((-32*G*Pi*(P + rho))/alphaQ + 
        (fR/alphaQ + 8*G*Pi*(-3*P + rho))**2)**1.5 + (72*alphaQ*G*Pi - 
        48*alphaQ*drhoP*G*Pi + 8*alphaQ*drhoP**2*G*Pi + ddrhoPP*(-2 + fR - 
        24*alphaQ*G*P*Pi + 8*alphaQ*G*Pi*rho))/(alphaQ*math.sqrt((-32*G*Pi*(P + 
        rho))/alphaQ + (fR/alphaQ + 8*G*Pi*(-3*P + rho))**2))))
    ddsigma1PP = (-((math.sqrt(2)*math.sqrt(alphaQ)*(dlambdaP*lambda1 - 
        4*(1 + drhoP)*G*Pi)**2)/(lambda1**2 - 8*G*Pi*(P + rho))**1.5) + 
        (math.sqrt(2)*math.sqrt(alphaQ)*(dlambdaP**2 + ddlambdaPP*lambda1 - 
        4*ddrhoPP*G*Pi))/math.sqrt(lambda1**2 - 8*G*Pi*(P + rho)) + 
        (8*a*ddrhoPP*G*Pi)/Rp)
    ddsigma2PP = math.sqrt(2)*math.sqrt(alphaQ)*ddlambdaPP + (8*a*ddrhoPP*G*Pi)/Rp

    ddOmegaPP = ((2*sigma1*sigma2*(2*dsigma1P*dsigma2P + ddsigma2PP*sigma1 + 
        ddsigma1PP*sigma2) - (dsigma2P*sigma1 + dsigma1P*sigma2)**2)/
        (4.*(sigma1*sigma2)**1.5))
    ddSPP = ((3*dsigma2P**2*sigma1**2 + 6*sigma1*(-(dsigma1P*dsigma2P) + 
        ddsigma2PP*sigma1)*sigma2 + (3*dsigma1P**2 - 
        2*ddsigma1PP*sigma1)*sigma2**2)/(4.*sigma1**2*math.sqrt(sigma1*sigma2)))

    dalpha = ((P + rho)*(-dOmegaPDOmega**2 + ddOmegaPP/Omega - dSPDS**2 + 
        ddSPP/S) + (1 + drhoP)*(dOmegaPDOmega + dSPDS))/2.  
    dbeta = (-(r*(-2*dOmegaPDOmega*dSP**2*Omega**2*(P + rho) + 
        2*ddSPP*dOmegaPDOmega*Omega**2*(P + rho)*S + 
        (dOmegaPDOmega*(1 + drhoP)*(dOmegaPDOmega + 2*dSPDS)*Omega**2 - 
        2*dOmegaP**2*(-2 + dOmegaPDOmega*P + dSPDS*P + (dOmegaPDOmega + dSPDS)*rho) + 
        2*ddOmegaPP*Omega*(-2 + dOmegaPDOmega*P + dSPDS*P + 
        (dOmegaPDOmega + dSPDS)*rho))*S**2))/(2.*Omega**2*S**2))
    
    dtau_ttP = ((2*G*Pi*sigma2**2*(-2*(-3 + drhoP)*(48*a*G*P*Pi - 
        16*a*G*Pi*rho - Rp)*sigma1*sigma2 - (72*a*G*P**2*Pi + 
        8*a*G*Pi*rho**2 - 3*P*Rp - 2*rho*Rp + rho*(-48*a*G*P*Pi + Rp))*
        (3*dsigma2P*sigma1 + dsigma1P*sigma2)))/(Rp*(sigma1*sigma2**3)**1.5))  
    dtau_rrP = ((2*G*Pi*sigma2**2*(2*(2 + ((-3 + drhoP)*(-48*a*G*P*Pi + 
        16*a*G*Pi*rho + Rp))/Rp)*sigma1*sigma2 + ((-8*a*G*Pi*rho**2 + 
        rho*(48*a*G*P*Pi - Rp) + P*(-72*a*G*P*Pi + Rp))*(3*dsigma2P*sigma1 + 
        dsigma1P*sigma2))/Rp))/(sigma1*sigma2**3)**1.5)
    
    Phi = tau_rr + (Omega*tau_tt)/S
    dPhiP = (S*(dtau_ttP*Omega + dtau_rrP*S) + (-(dSP*Omega) + dOmegaP*S)*tau_tt/S**2)    
       
    
    if(r <= 1*10**-6): 
        ddP0 = (dP0*(2 + (3*Phi*r**2)/(-4*M + Phi*r**3) + 
            dM*(2 + 1/(M - (Phi*r**3)/4.)) + 
            dP*((dPhiP*r**3)/(-4*M + Phi*r**3) + (1 + drhoP)/(P + rho))))
        ddP = dP*((dalpha*dP)/(1 - alpha) + ddP0/dP0 + 
            (beta*ddP0*s)/(2.*(math.sqrt(1 - beta*dP0) + s - beta*dP0*s)) + 
            (dbeta*dP*dP0*s)/(2.*(math.sqrt(1 - beta*dP0) + s - beta*dP0*s)))
        beast=r*dP*dOmegaPDOmega
    else:
        A=1-2*M/r
        ddP0 = (dP0*((2*(M - r))/(r*(-2*M + r)) + (3*Phi*r**2)/(-4*M + Phi*r**3) + 
            dM*(2/(-2*M + r) + 1/(M - (Phi*r**3)/4.)) + 
            dP*((dPhiP*r**3)/(-4*M + Phi*r**3) + (1 + drhoP)/(P + rho))))
        ddP = dP*((dalpha*dP)/(1 - alpha) + ddP0/dP0 + 
            (beta*ddP0*s)/(2.*(math.sqrt(1 - beta*dP0) + s - beta*dP0*s)) + 
            (dbeta*dP*dP0*s)/(2.*(math.sqrt(1 - beta*dP0) + s - beta*dP0*s))) 
        beast=((r*(2*M*(6*A*dOmegaP*dP*Omega*S + r*(-3*A*dOmegaP**2*dP**2*S + 
            4*A*(ddP*dOmegaP + ddOmegaPP*dP**2)*Omega*S +  6*Omega**2*S*tau_rr - 
            2*Omega**3*tau_tt)) + r*(-8*A*dOmegaP*dP*Omega*S + 
            r*(3*A*dOmegaP**2*dP**2*S - 4*A*(ddP*dOmegaP + ddOmegaPP*dP**2)*Omega*S - 
            6*Omega**2*S*tau_rr + 2*Omega**3*tau_tt))))/
            (4.*Omega*(2*M - r)*(2*Omega + dOmegaP*dP*r)*S))
    dm = beast/g_c2
    return dP, dm
 
    
rho_center = 1*(10**6)
rhocs = []
radii = []
masses = []
target = open("SimData.txt", 'w')
line = 'rho R(km) M(SM) \n'
target.write(line)
radius, mass = make_a_tov_star(rho_center)
#while rho_center <= 5.11*(10**16):
#    radius, mass = make_a_tov_star(rho_center)
#    print "rho_c=",rho_center," Radius=",radius,"Mass=",mass
#    rhocs.append(rho_center)
#    radii.append(radius)
#    masses.append(mass)
#    line = "%f %f %f \n"%(rho_center, radius, mass)
#    target.write(line)
##    rho_center = min(rho_center * 1.005, rho_center+1*10**13)
#    rho_center = rho_center*1.1
##    print "did it"
   
target.close()    
plt.plot(radii,masses)
plt.xscale('log')
plt.show