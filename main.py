# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 08:55:09 2024

@author: felix
"""


'''
To Dos:
    x hc implementieren
    x counter flow
    x druckverlust
    x fin efficiency
    - hot gas radiation hinzufügen
    - expand equilibrium flow --> physik anchvollziehen und anpassen
    - equilibrium bartz --> welche terme müssen da angepasst werden?
    - überprüfung geometrieannahmen und rechnungen
    x plot function
    - TODOs abarbeiten
    - korrelationen für rauigkeit hinzufügen (druckverlust coolant + wärmeübergang bei hot has und coolant)
    - geeingete korrelation für ethanol finden --> Max
    
    - area_wetted_hg ist nicht kleiner, sondern gleich groß wie area_wetted_inner (oder wie die variable heißt) --> erst ab bestimmtem index
    - druckverlust korrelation korrigieren
'''


import os 
os.chdir(r'C:\Users\felix\Desktop\Therm')

import therm_functions as tfc



# config = {'ID': 'Hopper 223kg 40% thrust',
#           'chamber': {
#                 'd_th': {'Value': 34.5e-3, 'Description': 'Throat Diameter [m]'},
#                 'l_star': {'Value': 1, 'Description': 'Characteristic Chamber Length including convergent section [m]'},
#                 'alpha': {'Value': 30, 'Description': 'Half angle Convergent section [°]'},
#                 'beta': {'Value': 15, 'Description': 'Half angle Divergent section [°]'},
#                 'epsilon_c': {'Value': 8, 'Description': 'Contraction Ratio [-]'},
#                 'epsilon_e': {'Value': 4.1, 'Description': 'Expansion Ratio [-]'}},
#           'cooling': {
#                 'N': {'Value': 93, 'Description': 'Number of channels [-]'},
#                 'h': {'Value': 10e-3, 'Description': 'height of cooling channel [m]'},
#                 'tc': {'Value': 1e-3, 'Description': 'wall thickness between cooling channels in circumferential direction [m]'},
#                 'tl': {'Value': 1e-3, 'Description': 'thickness of liner between chamber and cooling channel [m]'},
#                 'lam': {'Value': 25, 'Description': 'heat conductivity of wall material [W/m/K]'},
#                 'p_inlet': {'Value': 30e5, 'Description': 'inlet pressure into cooling channel [Pa]'},
#                 'T_wall': {'Value': 700, 'Description': ' in case "fix" is chosen as thermal type, this wall temperature is used to calculate the heat flux into the wall [K]'},
#                 'T_inlet': {'Value': 293, 'Description': 'inlet temperature into cooling channel [K]'},
#                 'fin_efficiency': {'Value': True, 'Description': 'Fin efficiency'}},
#           'nodes': {'Value': 100, 'Description': 'Number of axially discretized equidistant elements [-]'},
#           'load point': {
#                 'p_c': {'Value': 9.46e5, 'Description': 'total chamber pressure [Pa]'},
#                 # 'p_c': {'Value': 9.7e5, 'Description': 'total chamber pressure [Pa]'},
#                 'ROF': {'Value': 1.1, 'Description': 'oxidizer to fuel ratio [-]'},
#                 'm_dot': {'Value': 0.55, 'Description': 'total mass flow in combustion chamber [kg/s]'},
#                 'T_oxidizer': {'Value': 100, 'Description': 'Oxidizer Injection Temperature [K]'},
#                 'T_fuel': {'Value': 293, 'Description': 'Fuel Injection Temperature [K]'},
#                 # 'm_dot': {'Value': 0.544, 'Description': 'total mass flow in combustion chamber [kg/s]'},
#                 'fuel_cp': {'Value': 'C2H6O', 'Description': 'chemical formula of fuel as required by coolprop library'},
#                 'fuel_ct': {'Value': 'C2H5OH', 'Description': 'chemical formula of fuel as required by cantera library'}},
#           'settings': {
#                 'direction': {'Value': 'counter', 'Description': 'Flow direction of fuel in regenerative cooling channels with respect to hot gas flow direction; either "counter" or "co"'},
#                 'thermal type': {'Value': 'fix', 'Description': 'Either "fix" --> fixed wall temperature (needs to be specified) or "regen" --> solves regenerative cooling problem'},
#                 'flow type': {'Types': ['frozen'], 'Description': 'Either "frozen" or "equilibrium" or "both" --> affects the chemical composition along the axial direction AND transport properties (cp, lambda, mu etc.)'}
#               }}

config = {'ID': 'Hopper 223kg full thrust',
          'chamber': {
                'd_th': {'Value': 37e-3, 'Description': 'Throat Diameter [m]'},
                'l_star': {'Value': 1, 'Description': 'Characteristic Chamber Length including convergent section [m]'},
                'alpha': {'Value': 30, 'Description': 'Half angle Convergent section [°]'},
                'beta': {'Value': 15, 'Description': 'Half angle Divergent section [°]'},
                'epsilon_c': {'Value': 8, 'Description': 'Contraction Ratio [-]'},
                'epsilon_e': {'Value': 3, 'Description': 'Expansion Ratio [-]'}},
          'cooling': {
                'N': {'Value': 50, 'Description': 'Number of channels [-]'},
                'h': {'Value': 0.8e-3, 'Description': 'height of cooling channel [m]'},
                'tc': {'Value': 1.5e-3, 'Description': 'wall thickness between cooling channels in circumferential direction along arc [m]'},
                'tl': {'Value': 1e-3, 'Description': 'thickness of liner between chamber and cooling channel [m]'},
                'lam': {'Value': 25, 'Description': 'heat conductivity of wall material [W/m/K]'},
                'p_inlet': {'Value': 40e5, 'Description': 'inlet pressure into cooling channel [Pa]'},
                'T_wall': {'Value': 700, 'Description': ' in case "fix" is chosen as thermal type, this wall temperature is used to calculate the heat flux into the wall [K]'},
                'T_inlet': {'Value': 293, 'Description': 'inlet temperature into cooling channel [K]'},
                'fin_efficiency': {'Value': True, 'Description': 'Fin efficiency'}},
          'nodes': {'Value': 100, 'Description': 'Number of axially discretized equidistant elements [-]'},
          'load point': {
                'p_c': {'Value': 20e5, 'Description': 'total chamber pressure [Pa]'},
                # 'p_c': {'Value': 9.7e5, 'Description': 'total chamber pressure [Pa]'},
                'ROF': {'Value': 3, 'Description': 'oxidizer to fuel ratio [-]'},
                'm_dot': {'Value': 1.45, 'Description': 'total mass flow in combustion chamber [kg/s]'},
                'T_oxidizer': {'Value': 293, 'Description': 'Oxidizer Injection Temperature [K]'},
                'T_fuel': {'Value': 293, 'Description': 'Fuel Injection Temperature [K]'},
                # 'm_dot': {'Value': 0.544, 'Description': 'total mass flow in combustion chamber [kg/s]'},
                'oxidizer': {'Value': 'N2O', 'Description': 'chemical formula of oxidizer'},
                'fuel_cp': {'Value': 'C2H6O', 'Description': 'chemical formula of fuel as required by coolprop library'},
                'fuel_ct': {'Value': 'C2H5OH', 'Description': 'chemical formula of fuel as required by cantera library'}},
          'settings': {
                'direction': {'Value': 'counter', 'Description': 'Flow direction of fuel in regenerative cooling channels with respect to hot gas flow direction; either "counter" or "co"'},
                'thermal type': {'Value': 'fix', 'Description': 'Either "fix" --> fixed wall temperature (needs to be specified) or "regen" --> solves regenerative cooling problem'},
                'flow type': {'Types': ['frozen'], 'Description': 'Either "frozen" or "equilibrium" or "both" --> affects the chemical composition along the axial direction AND transport properties (cp, lambda, mu etc.)'}
              }}
#%%
geometry = tfc.calc_chamber_geometry(config)

geometry = tfc.discretize_chamber_geometry(config, geometry)

geometry = tfc.calc_discretized_cooling_geometry(config, geometry)

# tfc.plot_geometry(geometry, config)

#%%

fluid = tfc.set_thermal_problem(config)

fluid = tfc.calc_flow(config, geometry, fluid)


#%%
import time
t_start = time.time()
fluid = tfc.thermal_analysis(config, geometry, fluid)
t_end = time.time()
print(f'Process took {(t_end-t_start):.2f} [s]')

# tfc.plot_flow(geometry, fluid)

# tfc.plot_constant_wall(geometry, fluid)

# tfc.plot_regen_cooling(geometry, fluid)


#%%

# fc.plot_constant_wall(geometry, fluid)

# fc.plot_regen_cooling(geometry, fluid)

print(f'Heat absorbed: {tfc.obtain_heat_absorbed(fluid)/1e3} [kW]')

print(f'Heat Capacity {tfc.global_heat_capacity(config, p_outlet=28e5, T_outlet=["Margin", 0])[0]/1e3} [kW]')#, T_outlet=["Value", 450]))

# print(tfc.obtain_heat_of_combustion(config), tfc.obtain_heat_absorbed(fluid)/1e3/tfc.obtain_heat_of_combustion(config)*100)
#%% senity checks

# import matplotlib.pyplot as plt

# Q_dot_hg = fluid['cooling']['Q_dot_hg']['array']

# Q_dot_con = fluid['cooling']['Q_dot_conduction']['array']

# Q_dot_c = fluid['cooling']['Q_dot_c']['array']


# plt.plot(Q_dot_hg)
# plt.plot(Q_dot_con)
# plt.plot(Q_dot_c)

# print(sum(Q_dot_hg))


# m_dot = 1.08
# ROF = 1.4
# m_dot_fuel = m_dot*1/(1+ROF)

# deltah_array = fluid['cooling']['deltah']['array']

# print(sum(deltah_array)*m_dot_fuel)



# h2 = psi('H', 'P', 22.5e5, 'T', 420, 'C2H6O')
# h1 = psi('H', 'P', 40e5, 'T', 293, 'C2H6O')

# deltah_cumulative = h2 - h1

# print(deltah_cumulative*m_dot_fuel)