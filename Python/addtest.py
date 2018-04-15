# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 08:21:39 2016

@author: peterweilnboeck
"""

import decimal
Dec = decimal.Decimal
import mpmath

a = 21815755775858530320384.00
b = 2
print "a=%f"%a
print "b=%f"%b
print "c=%f"%(a+b)

print "-------------------------"
A = Dec('21815755775858530320384.00')
B = Dec('2.00')

print "a=",A
print "b=",B
C=A+B
print "c=",C
#print "-------------------------"
#A1 = Decimal.from_float(a)
#B1 = Decimal.from_float(b)
#
#print "a=",A1
#print "b=",B1
#C1=A1+B1
#print "c=",C1
#print "-------------------------"
#A2 = Decimal('%.2f'%a)
#B2 = Decimal('%.2f'%b)
#
#print "a=",A2
#print "b=",B2
#C2 = A2 + B2
#print "c=",C2
#
#print "-------------------------"
#mpmath.dps = 55
#A3 = mpmath.mpf(a)
#B3 = mpmath.mpf(b)
#
#mpmath.nprint(A3,55)
#mpmath.nprint(B3,50)
#C3 = A3 + B3
#mpmath.nprint(C3,50)
#
#print "-------------------------"
#mpmath.dps = 55
#A4 = mpmath.mpf(21815755775858530320384.00)
#B4 = mpmath.mpf(2.000)
#
#mpmath.nprint(A4,55)
#mpmath.nprint(B4,50)
#C4 = A4 + B4
mpmath.nprint(C4,50)