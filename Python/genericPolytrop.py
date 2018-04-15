# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 10:38:19 2016

@author: peterweilnboeck
"""
import math as math

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
global gamma
gamma = 2.75*10**0
global eosk
eosk = 1.98183*10**-6
global min_press
min_press = 76
global delta
delta = 0.23

def make_a_tov_star():
    rho_C = 5.0*10**14
    tov_integrate(rho_C)

def tov_integrate( roh_c ):
    press = eosk*(roh_c**gamma)  
    print press
    press_now = press
    m_now = 0.0
    count = 0
    m_now = 0.0
    rad = 0.0
    dr = 1.0*10**-20
    while press_now >= min_press and count < 1000000:
        count+=1
        press_now, m_now = tov_rk4(press_now, m_now, rad, dr)
        rad = rad + dr
        dpress, dm = tov_RHS(press_now, m_now, rad)
        print "P=",press_now
        print "M=",m_now
        print "R=",rad
        print "Count = ",count
        dr = delta*1/((1/m_now)*dm - (1/press_now)*dpress)
        
        
def tov_rk4(press_old, m_old, rad, dr):    
    press_temp = press_old
    m_temp = m_old
    press_k1, m_k1 = tov_RHS(press_temp, m_temp, rad)
    
    press_temp = press_old + 0.5*press_k1
    m_temp = m_old + 0.5*m_k1
    press_k2, m_k2 = tov_RHS(press_temp, m_temp, rad + 0.5*dr)
    
    press_temp = press_old + 0.5*press_k2
    m_temp = m_old + 0.5*m_k2
    press_k3, m_k3 = tov_RHS(press_temp, m_temp, rad + 0.5*dr)
    
    press_temp = press_old + press_k3
    m_temp = m_old + m_k3
    press_k4, m_k4 = tov_RHS(press_temp, m_temp, rad + dr)
    
    press_new = press_old + (press_k1 + press_k2 + press_k2 + press_k3 + press_k3 + press_k4) / 6
    m_new = m_old + (m_k1 + m_k2 + m_k2 + m_k3 + m_k3 + m_k4) / 6
    
    return press_new, m_new
    
def rk4(f, x0, y0, x1, n):
    vx = [0] * (n + 1)
    vy = [0] * (n + 1)
    h = (x1 - x0) / float(n)
    vx[0] = x = x0
    vy[0] = y = y0
    for i in range(1, n + 1):
        k1 = h * f(x, y)
        k2 = h * f(x + 0.5 * h, y + 0.5 * k1)
        k3 = h * f(x + 0.5 * h, y + 0.5 * k2)
        k4 = h * f(x + h, y + k3)
        vx[i] = x = x0 + i * h
        vy[i] = y = y + (k1 + k2 + k2 + k3 + k3 + k4) / 6
    return vx, vy
 
def f(x, y):
    return x * math.sqrt(y)

def tov_RHS(press, mgrav, r):
    press = max(press, min_press)
    rho = (press / eosk)**(1.0 / gamma)
    if(r <= 1*10**-6): 
        beast=(4.0*pi*r*press/clite**2)
    else:
        beast=(mgrav+4.0*pi*r*r*r*press/clite**2)/(r*(r-2.0*ggrav*mgrav/clite**2))

    dpress = -ggrav*( rho+press/clite**2 )*beast
    dm = 4.0*pi*r*r*rho
    return dpress, dm
 
    
make_a_tov_star()
vx, vy = rk4(f, 0, 1, 10, 100)
for x, y in list(zip(vx, vy))[::10]:
    print("%4.1f %10.5f %+12.4e" % (x, y, y - (4 + x * x)**2 / 16))