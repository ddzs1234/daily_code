#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 13:54:45 2020

@author: ashley
"""
# not finished yet

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


Bh=13.6
yita=10**-10
me=9.1*10**-27

left=0

func=lambda t: left-(0.26*(me/t)**(2/3)*np.exp(-Bh/t)/yita)

# t_initial_guess=0
# t_solve=fsolve(func,t_initial_guess)

# print('solution : ',t_solve)
