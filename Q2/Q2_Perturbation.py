# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 20:05:20 2017

@author: amit
"""

#imports
from scipy.integrate import nquad
from scipy.special import hermite
import numpy as np
from numpy import sqrt, exp, sin, cos, tan, pi, inf
from math import factorial
import matplotlib.pyplot as plt
#Global constants
mass = 1
hbar = 1

A = 1 / 2000

def HermitePoly(x, n):
    if n == 0:
        return 1
    if n == 1:
        return 2 * x
    if n == 2:
        return 4 * pow(x,2) - 2
    if n == 3:
        return 8 * pow(x,3) - 12 * x
    if n == 4:
        return 16 * pow(x,4) - 48 * pow(x,2) + 12
    if n == 5:
        return 32 * pow(x,5) - 160 * pow(x,3) + 120 * x
    if n == 6:
        return 64 * pow(x,6) - 480 * pow(x,4) + 720 * pow(x,2) - 120
    if n == 7:
        return 128 * pow(x,7) - 1344 * pow(x,5) + 3360 * pow(x,3) - 1680 * x
    if n == 8:
        return 256 * pow(x,8) - 3584 * pow(x,6) + 13440 * pow(x,4) - 13440 * pow(x,2) + 1680
    if n == 9:
        return 512 * pow(x,9) - 9216 * pow(x,7) + 48384 * pow(x,5) - 80640 * pow(x,3) + 30240
    if n == 10:
        return 1024 * pow(x,10) - 23040 * pow(x,8) + 161280 * pow(x,6) - 403200 * pow(x, 4) + 302400 * pow(x,2) - 30240  

def HOS(x, n, omega):
    alpha = mass * omega / hbar
    return ( 1/sqrt( pow(2,n) * factorial(n) ) ) * pow(alpha/pi, 1/4) * exp( -alpha * pow(x, 2) / 2) * HermitePoly(sqrt(alpha) * x, n)

def psi(x1, x2, omega):
    return pow(1 / 2, 1 / 2) * ( HOS(x1, 0, omega) * HOS(x2, 1, omega) - HOS(x2, 0, omega) * HOS(x1, 1, omega) )

def first_order_correction(x1, x2, omega):
    return psi(x1, x2, omega) * pow( abs(x1 - x2) + A, -1 ) * psi(x1, x2, omega)

omega_list = []
correced_E_list = []

omega = 1
while(omega <= 5):
    
    E_c = nquad( first_order_correction, [ [ -inf, inf ], [ -inf, inf ] ], args= [omega] )
    print(omega, '\t', E_c[0])
    omega_list.append(omega)
    correced_E_list.append(E_c[0])
    omega+= 0.5

plt.plot(omega_list, correced_E_list)
plt.show()