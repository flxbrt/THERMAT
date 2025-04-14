# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 13:06:10 2024

@author: felix
"""



import numpy as np
import therm_functions as tfc

config = {'ID': 'comparison_RPA',
           'chamber': {
                 'd_th': {'Value': 50e-3, 'Description': 'Throat Diameter [m]'},
                 'l_star': {'Value': 1, 'Description': 'Characteristic Chamber Length including convergent section [m]'},
                 'alpha': {'Value': 30, 'Description': 'Half angle Convergent section [°]'},
                 'beta': {'Value': 15, 'Description': 'Half angle Divergent section [°]'},
                 'epsilon_c': {'Value': 5, 'Description': 'Contraction Ratio [-]'},
                 'epsilon_e': {'Value': 5, 'Description': 'Expansion Ratio [-]'}},
          'cooling': {
                'N': {'Value': 50, 'Description': 'Number of channels [-]'},
                
                'h': {'Value': [1e-3], 'Description': 'height of cooling channel [m]'},
                
                # 'h': {'Value': np.array([2e-3]), 'Description': 'height of cooling channel [m]'},
                'tc': {'Value': 2e-3, 'Description': 'wall thickness between cooling channels in circumferential direction [m]'},
                'tl': {'Value': 2e-3, 'Description': 'thickness of liner between chamber and cooling channel [m]'},
                'lam': {'Value': 100, 'Description': 'heat conductivity of wall material [W/m/K]'},
                # https://iopscience.iop.org/article/10.1088/1742-6596/1382/1/012175/pdf
                'p_inlet': {'Value': 30e5, 'Description': 'inlet pressure into cooling channel [Pa]'},
                'T_wall': {'Value': 700, 'Description': ' in case "fix" is chosen as thermal type, this wall temperature is used to calculate the heat flux into the wall [K]'},
                'T_inlet': {'Value': 293, 'Description': 'inlet temperature into cooling channel [K]'},
                'fin_efficiency': {'Value': True, 'Description': 'Fin efficiency'}},
          'nodes': {'Value': 53, 'Description': 'Number of axially discretized equidistant elements [-]'},
          'load point': {
                'p_c': {'Value': 20e5, 'Description': 'total chamber pressure [Pa]'},
                # 'p_c': {'Value': 9.7e5, 'Description': 'total chamber pressure [Pa]'},
                'ROF': {'Value': 1, 'Description': 'oxidizer to fuel ratio [-]'},
                ### massenstrom eig. 1.63kg/s
                ### failed aber nur nicht, wenn ich auf 1.52kg/s runter gehe
                'm_dot': {'Value': 2.55, 'Description': 'total mass flow in combustion chamber [kg/s]'},
                'T_oxidizer': {'Value': 100, 'Description': 'Oxidizer Injection Temperature [K]'},
                'T_fuel': {'Value': 300, 'Description': 'Fuel Injection Temperature [K]'},
                # 'm_dot': {'Value': 0.544, 'Description': 'total mass flow in combustion chamber [kg/s]'},
                'oxidizer': {'Value': 'O2', 'Description': 'chemical formula of oxidizer'},
                'fuel_cp': {'Value': 'C2H6O', 'Description': 'chemical formula of fuel as required by coolprop library'},
                'fuel_ct': {'Value': 'C2H5OH', 'Description': 'chemical formula of fuel as required by cantera library'}},
          'settings': {
                'direction': {'Value': 'co', 'Description': 'Flow direction of fuel in regenerative cooling channels with respect to hot gas flow direction; either "counter" or "co"'},
                'thermal type': {'Value': 'regen', 'Description': 'Either "fix" --> fixed wall temperature (needs to be specified) or "regen" --> solves regenerative cooling problem'},
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

fluid = tfc.thermal_analysis(config, geometry, fluid)

# tfc.plot_flow(geometry, fluid)

# tfc.plot_constant_wall(geometry, fluid)

tfc.plot_regen_cooling(geometry, fluid)
