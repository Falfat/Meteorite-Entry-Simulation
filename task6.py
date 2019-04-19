# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 13:29:23 2018

@author: Gareth
"""


import Read_AltitudeDensityTable as rd
import scipy.optimize as so
import meteor as m
import numpy as np
import matplotlib.pyplot as plt

alt, energy = rd.readcurvedata('data/ChelyabinskEnergyAltitude.csv')
maxe = energy.max()
maxalt = alt[energy.argmax()]


m1 = m.meteor(v=19.16E3, m=1200E3, theta=(18/180) * np.pi, r=30, sigma0 = 3E6) # may be useful

# fire iterate now returns peak enery loss position and maximum energy loss.

def minimisefunc(x):
    r = x[0]
    sigma0 = x[1]
    maxenergydiff, maxenergylossalt, _, _, _ = m1.fireiterate(r, sigma0)
    val = np.sqrt((maxe - maxenergydiff)**2 + (maxalt - maxenergylossalt)**2)
    return val

# need intial guess for r and sigma0
#answer = so.minimize(minimisefunc, [15, 2E6])

def minimiseplot(r=12.1, sigma0=3e6):
    alt, energy = rd.readcurvedata('data/ChelyabinskEnergyAltitude.csv')
    maxe = energy.max()
    maxalt = alt[energy.argmax()]
    z, energydiff, v, m, r = m1.fireiterate(r, sigma0)

    plt.figure()

    plt.plot(energydiff, z, label = r'Mathilde simulation $r \ = \ 8.5m, \ \sigma_{0} = 4 \times 10^{6}kg$')
    plt.xlabel(r'Energy release $\frac{kt \ TNT}{km}$')
    plt.ylabel(r'Altitude $km$')
    plt.title('Comparison of energy release curve with Chelyabinsk impact')
    plt.plot(energy, alt, label = 'Chelyabinsk empirical data')
    plt.legend(loc = 'best')
    plt.xlim(left = 0)


minimiseplot(8.5, 4e6)

#popt, pcov = so.curve_fit(m1.fireiterate, alt, energy)
