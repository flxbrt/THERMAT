# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 20:56:30 2024

@author: felix
"""



mu = fluid['combustion']['frozen']['mu']['array']
rho = fluid['combustion']['frozen']['rho static']['array']
vel = fluid['combustion']['frozen']['velocity']['array']

axial = geometry['1D']['axial']['array']
delta_x = axial[1] - axial[0]

delta_vel = vel[1:]-vel[:-1]
dvel_dx = delta_vel/delta_x

K = mu[1:]/rho[1:]/vel[1:]**2*dvel_dx