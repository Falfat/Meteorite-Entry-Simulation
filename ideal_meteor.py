# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 16:50:16 2018

@author: ming-rui
"""
import numpy as np
import scipy.integrate as sp
import MeteorPackage.meteor as mt
import matplotlib.pyplot as plt
import sympy as sy
from sympy import Function, dsolve, Eq, sin, cos, exp, symbols

class ideal_meteor(mt.meteor):
    
    def __init__(self, v, m, theta, r, pm=3300, z0=100E3, g=0, Rp=6400E30, sigma0=np.inf, pa_tab=False, planet='Earth'):
        self.planet = planet
        self.v = v
        self.m = m
        self.theta = theta * np.pi/180
        self.z = z0
        self.x = 0
        self.r = r
        self.z0 = z0
        self.g = g
        self.Rp = Rp  # planetary radius
        self.Cd = 1
        self.Ch = 0
        self.Q = 10**7
        self.Cl = 0
        self.a = 0
        self.H = 8000
        self.sigma0 = sigma0
        self.pm = pm
        self.A = np.pi * self.r**2
        self.p0 = 1.2  # atmospheric pressure at z0
        self.pa_tab = pa_tab

    def func(self, t, y):
        """
        rewrite differentiable function pre ablation for scipy solver.
        """
        v = y[0]
        m = y[1]
        theta = y[2]
        z = y[3]
        x = y[4]
        r = y[5]
        pa = self.p0 *np.exp(-z / self.H) # added to try and fix the issue with graph

        dv = ((-self.Cd * pa * self.A * v**2) / (2*m))
        dm = 0
        dtheta = 0
        dz = -v * np.sin(theta)
        dx = v * np.cos(theta)
        dr = 0

        return [dv, dm, dtheta, dz, dx, dr]  

    def compare(self, plot=False):
        res = super(ideal_meteor, self).analytical()
        time = res[2][0]
        data_num = res[2]
        data_ana = res[:2]
        if plot:
            fig, ax = plt.subplots(1, 1, figsize=(8, 6))
            ax.plot(data_num[1][0], data_num[1][3],'r', label = '$numerical$')
            ax.plot(data_ana[0], data_ana[1],'b--', label = '$analytical$')
            ax.set_xlabel('$v$', fontsize=14)
            ax.set_ylabel('$z$', fontsize=14)
            ax.set_title('Compare the numerical and analytical solutions', fontsize=16)
            ax.legend(loc='best', fontsize=12)
            plt.show()

        return data_num[1][0], data_num[1][3], data_ana[0], data_ana[1], time

    def ideal_plot(self):
        super(ideal_meteor, self).plot()



tt = ideal_meteor(v=19.2E3, m=1200E3, theta=18.3, r=50, z0=100E3, pa_tab=False)
x11, y11, x21, y21, t1 = tt.compare()

tt = ideal_meteor(v=19.2E3, m=1200E3, theta=18.3, r=20, z0=100E3, pa_tab=False)
x12, y12, x22, y22, t2 = tt.compare()

tt = ideal_meteor(v=19.2E3, m=1200E3, theta=18.3, r=10, z0=100E3, pa_tab=False)
x13, y13, x23, y23, t3 = tt.compare()

tt = ideal_meteor(v=19.2E3, m=1200E3, theta=18.3, r=2, z0=100E3, pa_tab=False)
x14, y14, x24, y24, t4 = tt.compare()

tt = ideal_meteor(v=19.2E3, m=1200E3, theta=18.3, r=100, z0=100E3, pa_tab=False)
x10, y10, x20, y20, t0 = tt.compare()



fig, ax = plt.subplots(1, 1, figsize=(8, 6))

ax.plot(x11/1E3, y11/1E3,'r', label='$numerical-ivp$', ms=2)
ax.text(10, 45, '$r=50$')
ax.plot(x21/1E3, y21/1E3,'b--', label='$analytical$',linewidth=2.5)

ax.plot(x12/1E3, y12/1E3,'r', ms=2)
ax.text(10, 29, '$r=20$')
ax.plot(x22/1E3, y22/1E3,'b--', linewidth=2.5)

ax.plot(x13/1E3, y13/1E3,'r', ms=2)
ax.text(10, 18, '$r=10$')
ax.plot(x23/1E3, y23/1E3,'b--', linewidth=2.5)

ax.plot(x14/1E3, y14/1E3,'r', ms=2)
ax.text(16.5, 8, '$r=2$')
ax.plot(x24/1E3, y24/1E3,'b--', linewidth=2.5)

ax.plot(x10/1E3, y10/1E3,'r', ms=2)
ax.text(10, 57, '$r=100$')
ax.plot(x20/1E3, y20/1E3,'b--', linewidth=2.5)


ax.set_xlabel('$v(km/s)$', fontsize=14)
ax.set_ylabel('$z(km)$', fontsize=14)
ax.set_title('Compare the numerical and analytical solutions', fontsize=16)
ax.legend(loc='best', fontsize=12)
ax.grid()
plt.show()



#absolutely error
fig1, ax1 = plt.subplots(1, 1, figsize=(8, 6))
ax1.plot(t1[10:-1], abs((x11-x21)[10:-1]),'r', label='$v-error$', linewidth=2, ms=2)
ax1.plot(t1[10:-1], abs((y11-y21)[10:-1]),'b', label='$z-error$', linewidth=2, ms=2)
ax1.set_yscale('log')
ax1.set_xlabel('$t$', fontsize=14)
ax1.set_ylabel('$error$', fontsize=14)
ax1.set_title('$Absolute Errors$', fontsize=16)
ax1.legend(loc='best', fontsize=12)
ax1.grid()
plt.show()

#relative error
fig2, ax2 = plt.subplots(1, 1, figsize=(8, 6))
ax2.plot(t1, abs((x11-x21))*100/(x21),'r', label = '$v-error$', linewidth=2, ms=2)
ax2.plot(t1, abs((y11-y21))*100/(y21),'b', label = '$z-error$', linewidth=2, ms=2)
ax2.set_xlabel('$t$', fontsize=14)
ax2.set_ylabel('$error(percent)$', fontsize=14)
ax2.set_title('$Relative Errors$', fontsize=16)
ax2.legend(loc='best', fontsize=12)
ax2.grid()
plt.show()

#v-error-altitude
fig3, ax3 = plt.subplots(1, 1, figsize=(8, 6))
ax3.plot((x11-x21), y11/1E3,'r', label = '$v-error$', linewidth=2, ms=2)
#ax2.plot(abs((y11-y21)), y11,'bx', label = '$z-error$', linewidth=2, ms=2)
ax3.set_xlabel('$error(m/s)$', fontsize=14)
ax3.set_ylabel('$Altitude(km)$', fontsize=14)
ax3.set_title('$Error$$-$$Altitude$', fontsize=16)
ax3.legend(loc='best', fontsize=12)
ax3.grid()
plt.show()

#z-error-altitude
fig3, ax3 = plt.subplots(1, 1, figsize=(8, 6))
ax3.plot((y11-y21), y11/1E3,'b', label = '$z-error$', linewidth=2, ms=2)
#ax2.plot(abs((y11-y21)), y11,'bx', label = '$z-error$', linewidth=2, ms=2)
ax3.set_xlabel('$error(m)$', fontsize=14)
ax3.set_ylabel('$Altitude(km)$', fontsize=14)
ax3.set_title('$Error$$-$$Altitude$', fontsize=16)
ax3.legend(loc='best', fontsize=12)
ax3.grid()
plt.show()

