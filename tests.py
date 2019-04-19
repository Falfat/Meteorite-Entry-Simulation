# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 15:24:17 2018

@author: Gareth
"""

import meteor as m
import numpy as np


def test_notablated(): #checks if doesnt ablate
    m1 = m.meteor(v=19.16E3, m=1200E3, theta=(18/180) * np.pi, r=30, sigma0 = 64E6)
    time, yfinal, ablated, t_events = m1.fire()
    assert(ablated is False)
    assert(len(t_events) > 0)

def test_ablated(): # checks if does ablate
    m1 = m.meteor(v=19.16E3, m=1200E3, theta=(18/180) * np.pi, r=30, sigma0 = 64E2)
    time, yfinal, ablated, t_events = m1.fire()
    assert(ablated is True)
    assert(len(t_events) > 0)

def test_energypatch():
    m1 = m.meteor(v=19.16E3, m=1200E3, theta=(18/180) * np.pi, r=30, sigma0 = 64E6)
    time, energydiff, _ = m1.energypatch()

    assert(len(time) == len(energydiff))


