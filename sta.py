#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 20:16:40 2018
@author: ming-rui
"""

import numpy as np
import matplotlib.pyplot as plt
import meteor as mt

def statistic(v, theta, r, sigma_v, sigma_theta, sigma_r, PLOT=True):
    '''
    Calculate the uncertainty in altitude and magnitude of peak energy release 
    --------
    Parameters:
    mean value and sigma of v, theta, r
    --------
    Return:
    Distribution of the peak energy loss and altitude of the peak point
    mean value and sigma of both distributions
    
    '''
    mu_v = v#19E3
    
    sigma_v = sigma_v#0.05*mu_v
    sampleNo = 1000
    np.random.seed(0)
    v = np.random.normal(mu_v, sigma_v, sampleNo)#Distribution v
    
    mu_theta = theta/180*np.pi#18.3/180 * np.pi
    
    sigma_theta = sigma_theta#0.1*mu_theta
    sampleNo = 1000
    np.random.seed(0)
    theta = np.random.normal(mu_theta, sigma_theta, sampleNo)#Distribution theta
    
    mu_r = r#8.5
    
    sigma_r = sigma_r#0.1*mu_r
    sampleNo = 1000
    np.random.seed(0)
    r = np.random.normal(mu_r, sigma_r, sampleNo)#Distribution  r
    
    #Distribution mass
    mass = (4/3) * np.pi * r**3 * 3300
    #mu_mass = (4/3) * np.pi * mu_r**3 * 3300
    
    #######################################################
    
    data_set = zip(v, theta, mass, r)
    max_energyloss = []
    loss_altitude = []
    ssum = 0
    for v, theta, mass, r in data_set:
        sta = mt.meteor(v=v, m=mass, theta=theta, r=r, sigma0=4e6, pa_tab=False, planet="Earth")
        ediff, maxenergyloss, elossaltitude = sta.energyrelease()
        max_energyloss.append(maxenergyloss)
        loss_altitude.append(elossaltitude)
        ssum += 1
        print('boom', ssum)
    
    if PLOT == True:
        max_energyloss = np.array(max_energyloss)
        mu_me = max_energyloss.mean()
        print('mu_me', mu_me)
        sigama_me = max_energyloss.std()**2
        print('sigam_me', sigama_me)
        
        loss_altitude = np.array(loss_altitude)
        mu_la = loss_altitude.mean()
        print('mu_la', mu_la)
        sigama_la = loss_altitude.std()**2
        print('sigam_la', sigama_la)
        
        n, bins, patch = plt.hist(x=max_energyloss, bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85, density=True)
        plt.grid(axis='y', alpha=0.75)
        plt.xlabel('$energy$ $loss (kt)$')
        plt.ylabel('$Probability$')
        plt.title('$Distribution$ $of$ $Max$ $Energy$ $Loss$')
        plt.xlim(xmax=200)
        plt.show()
        
        #weights = np.ones_like(loss_altitude)/float(len(loss_altitude)) #if bins='auto' not work, use these weights
        n1, bins1, patch1 = plt.hist(loss_altitude,bins='auto', color='#0504aa', alpha=0.7, rwidth=0.8, density=True)
        plt.grid(axis='y', alpha=0.75)
        plt.xlabel('$Altitude$ $of$ $max$ $engergy$ $loss$ $point$ $(km)$')
        plt.ylabel('$Probability$')
        plt.title('$Distribution$ $of$ $Altitude$')
        plt.xlim(xmin=25, xmax=35)
        plt.show()
        
    return max_energyloss, mu_me, sigama_me, loss_altitude, mu_la, sigama_la






