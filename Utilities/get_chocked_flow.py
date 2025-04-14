# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 18:58:50 2024

@author: felix
"""



import cantera as ct
import numpy as np

def set_combustion():
    comb = ct.Solution('gri30_WARR.yaml')
    return comb

def obtain_combustion_properties(comb, p):
    def heat_capacity_ratio(comb):
        return comb.cp_mass / comb.cv_mass
    
    ref_t = 273
    ROF = 6

    comb.Y = {'CH4': 1, 'O2': ROF}
    comb.TP = ref_t, p
    comb.equilibrate('HP')
    
    return comb.T, comb.mean_molecular_weight, heat_capacity_ratio(comb)

def calc_chocked_mass_flow(d_th, gamma, M, T, p_cc):
        # theta = np.sqrt(gamma*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))))
        theta = np.sqrt(gamma)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))
        R_ideal = 8314
        m_dot = np.pi*d_th**2/4*p_cc/np.sqrt(T*R_ideal/M)*theta
        return m_dot

comb = set_combustion()
T, M, gamma = obtain_combustion_properties(comb, p=19e5)
m_dot = calc_chocked_mass_flow(7.3e-3, gamma, M, T, p_cc=19e5)

# m_dot_nitrogen = calc_chocked_mass_flow(4e-3, 1.4, 28, 270, 300e5)

