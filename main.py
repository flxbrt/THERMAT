"""
Title: main.py
Project: THERMAT - Thermal Analysis Tool for Rocket Engine Thrust Chambers
Author: @flxbrt
Version: 2.0

Description:
    Entry point script for THERMAT, which models thermal transport in rocket engine thrust chambers 
    and supports the design of regenerative cooling channels.

Usage:
    python main.py

Dependencies:
    - core.therm_functions
"""

import os
# Get current file directory path
path = os.path.dirname(__file__)

#%% Imports and configuration

import core.therm_functions as tfc

# Define full system configuration for the chamber, cooling, and load point
# Includes geometry, fluid properties, and model settings
# config = {
#     'ID': 'Hopper 223kg full thrust',
#     'chamber': {
#         'd_th': {'Value': 37e-3, 'Description': 'Throat Diameter [m]'},
#         'l_star': {'Value': 1, 'Description': 'Characteristic Chamber Length [m]'},
#         'alpha': {'Value': 30, 'Description': 'Half angle Convergent section [°]'},
#         'beta': {'Value': 15, 'Description': 'Half angle Divergent section [°]'},
#         'epsilon_c': {'Value': 8, 'Description': 'Contraction Ratio [-]'},
#         'epsilon_e': {'Value': 3, 'Description': 'Expansion Ratio [-]'}
#     },
#     'cooling': {
#         'N': {'Value': 50, 'Description': 'Number of channels [-]'},
#         'h': {'Value': 0.8e-3, 'Description': 'Channel height [m]'},
#         'tc': {'Value': 1.5e-3, 'Description': 'Wall thickness circumferentially [m]'},
#         'tl': {'Value': 1e-3, 'Description': 'Wall thickness liner [m]'},
#         'lam': {'Value': 25, 'Description': 'Wall material thermal conductivity [W/m/K]'},
#         'p_inlet': {'Value': 40e5, 'Description': 'Coolant inlet pressure [Pa]'},
#         'T_wall': {'Value': 700, 'Description': 'Fixed wall temperature if used [K]'},
#         'T_inlet': {'Value': 293, 'Description': 'Coolant inlet temperature [K]'},
#         'fin_efficiency': {'Value': True, 'Description': 'Whether to consider fin efficiency'}
#     },
#     'nodes': {'Value': 100, 'Description': 'Number of axial nodes [-]'},
#     'load point': {
#         'p_c': {'Value': 20e5, 'Description': 'Chamber pressure [Pa]'},
#         'ROF': {'Value': 3, 'Description': 'Oxidizer-to-fuel ratio [-]'},
#         'm_dot': {'Value': 1.45, 'Description': 'Mass flow rate [kg/s]'},
#         'T_oxidizer': {'Value': 293, 'Description': 'Oxidizer injection temperature [K]'},
#         'T_fuel': {'Value': 293, 'Description': 'Fuel injection temperature [K]'},
#         'oxidizer': {'Value': 'N2O', 'Description': 'Oxidizer chemical formula'},
#         'fuel_cp': {'Value': 'C2H6O', 'Description': 'Fuel formula (CoolProp)'},
#         'fuel_ct': {'Value': 'C2H5OH', 'Description': 'Fuel formula (Cantera)'}
#     },
#     'settings': {
#         'direction': {'Value': 'counter', 'Description': 'Flow direction of coolant (co/counter)'},
#         'thermal type': {'Value': 'fix', 'Description': '"fix" uses constant wall temperature'},
#         'flow type': {
#             'Types': ['frozen'],
#             'Description': 'Frozen/equilibrium chemistry assumptions for flow'
#         }
#     }
# }

config = {
    'ID': 'Pintle High Load Point',
    'chamber': {
        'd_th': {'Value': 0.03563152896630277, 'Description': 'Throat Diameter [m]'},
        'l_star': {'Value': 2, 'Description': 'Characteristic Chamber Length [m]'},
        'alpha': {'Value': 45, 'Description': 'Half angle Convergent section [°]'},
        'beta': {'Value': 15, 'Description': 'Half angle Divergent section [°]'},
        'epsilon_c': {'Value': 16.33, 'Description': 'Contraction Ratio [-]'},
        'epsilon_e': {'Value': 2.87, 'Description': 'Expansion Ratio [-]'}
    },
    'cooling': {
        'N': {'Value': 50, 'Description': 'Number of channels [-]'},
        'h': {'Value': 0.8e-3, 'Description': 'Channel height [m]'},
        'tc': {'Value': 1.5e-3, 'Description': 'Wall thickness circumferentially [m]'},
        'tl': {'Value': 1e-3, 'Description': 'Wall thickness liner [m]'},
        'lam': {'Value': 15, 'Description': 'Wall material thermal conductivity [W/m/K]'},
        'p_inlet': {'Value': 40e5, 'Description': 'Coolant inlet pressure [Pa]'},
        'T_wall': {'Value': 300, 'Description': 'Fixed wall temperature if used [K]'},
        'T_inlet': {'Value': 293, 'Description': 'Coolant inlet temperature [K]'},
        'fin_efficiency': {'Value': True, 'Description': 'Whether to consider fin efficiency'}
    },
    'nodes': {'Value': 100, 'Description': 'Number of axial nodes [-]'},
    'load point': {
        'p_c': {'Value': 20e5, 'Description': 'Chamber pressure [Pa]'},
        'ROF': {'Value': 1.17, 'Description': 'Oxidizer-to-fuel ratio [-]'},
        'm_dot': {'Value': 1.2551501159552974, 'Description': 'Mass flow rate [kg/s]'},
        'T_oxidizer': {'Value': 100, 'Description': 'Oxidizer injection temperature [K]'},
        'T_fuel': {'Value': 400, 'Description': 'Fuel injection temperature [K]'},
        'oxidizer': {'Value': 'O2', 'Description': 'Oxidizer chemical formula'},
        'fuel_cp': {'Value': 'C2H6O', 'Description': 'Fuel formula (CoolProp)'},
        'fuel_ct': {'Value': 'C2H5OH', 'Description': 'Fuel formula (Cantera)'}
    },
    'settings': {
        'direction': {'Value': 'counter', 'Description': 'Flow direction of coolant (co/counter)'},
        'thermal type': {'Value': 'fix', 'Description': '"fix" uses constant wall temperature'},
        'flow type': {
            'Types': ['frozen'],
            'Description': 'Frozen/equilibrium chemistry assumptions for flow'
        }
    }
}

# config = {
#     'ID': 'Pintle Low Load Point',
#     'chamber': {
#         'd_th': {'Value': 0.03563152896630277, 'Description': 'Throat Diameter [m]'},
#         'l_star': {'Value': 2, 'Description': 'Characteristic Chamber Length [m]'},
#         'alpha': {'Value': 45, 'Description': 'Half angle Convergent section [°]'},
#         'beta': {'Value': 15, 'Description': 'Half angle Divergent section [°]'},
#         'epsilon_c': {'Value': 15.44, 'Description': 'Contraction Ratio [-]'},
#         'epsilon_e': {'Value': 3, 'Description': 'Expansion Ratio [-]'}
#     },
#     'cooling': {
#         'N': {'Value': 50, 'Description': 'Number of channels [-]'},
#         'h': {'Value': 0.8e-3, 'Description': 'Channel height [m]'},
#         'tc': {'Value': 1.5e-3, 'Description': 'Wall thickness circumferentially [m]'},
#         'tl': {'Value': 1e-3, 'Description': 'Wall thickness liner [m]'},
#         'lam': {'Value': 25, 'Description': 'Wall material thermal conductivity [W/m/K]'},
#         'p_inlet': {'Value': 40e5, 'Description': 'Coolant inlet pressure [Pa]'},
#         'T_wall': {'Value': 700, 'Description': 'Fixed wall temperature if used [K]'},
#         'T_inlet': {'Value': 293, 'Description': 'Coolant inlet temperature [K]'},
#         'fin_efficiency': {'Value': True, 'Description': 'Whether to consider fin efficiency'}
#     },
#     'nodes': {'Value': 100, 'Description': 'Number of axial nodes [-]'},
#     'load point': {
#         'p_c': {'Value': 925779.4962934465, 'Description': 'Chamber pressure [Pa]'},
#         'ROF': {'Value': 0.92, 'Description': 'Oxidizer-to-fuel ratio [-]'},
#         'm_dot': {'Value': 0.6222634347397111, 'Description': 'Mass flow rate [kg/s]'},
#         'T_oxidizer': {'Value': 100, 'Description': 'Oxidizer injection temperature [K]'},
#         'T_fuel': {'Value': 400, 'Description': 'Fuel injection temperature [K]'},
#         'oxidizer': {'Value': 'O2', 'Description': 'Oxidizer chemical formula'},
#         'fuel_cp': {'Value': 'C2H6O', 'Description': 'Fuel formula (CoolProp)'},
#         'fuel_ct': {'Value': 'C2H5OH', 'Description': 'Fuel formula (Cantera)'}
#     },
#     'settings': {
#         'direction': {'Value': 'counter', 'Description': 'Flow direction of coolant (co/counter)'},
#         'thermal type': {'Value': 'fix', 'Description': '"fix" uses constant wall temperature'},
#         'flow type': {
#             'Types': ['frozen'],
#             'Description': 'Frozen/equilibrium chemistry assumptions for flow'
#         }
#     }
# }

#%% Geometry and discretization

# Calculate chamber geometry from throat, L*, and area ratios
geometry = tfc.calc_chamber_geometry(config)

# Discretize chamber geometry axially
geometry = tfc.discretize_chamber_geometry(config, geometry)

# Discretize cooling channel geometry over chamber nodes
geometry = tfc.calc_discretized_cooling_geometry(config, geometry)

# Optional visualization of chamber and cooling layout
tfc.plot_geometry(geometry, config)

#%% Flow setup and calculation

# Set up fluid properties and flow configuration
fluid = tfc.set_thermal_problem(config)

# Solve isentropic axial flow properties over chamber
fluid = tfc.calc_flow(config, geometry, fluid)

#%% Thermal simulation

import time
t_start = time.time()

# Perform thermal analysis (either constant wall or regenerative cooling)
fluid = tfc.thermal_analysis(config, geometry, fluid)

t_end = time.time()
print(f'Process took {(t_end-t_start):.2f} [s]')

#%% Results and plotting

# Show temperature, heat flux and other distributions along chamber
tfc.plot_flow(geometry, fluid)
tfc.plot_constant_wall(geometry, fluid)

# Uncomment to visualize regenerative cooling specifically
# tfc.plot_regen_cooling(geometry, fluid)

#%% Heat balance outputs

# Show absorbed heat [kW]
print(f'Heat absorbed: {tfc.obtain_heat_absorbed(fluid)/1e3} [kW]')

# Show heat capacity of coolant [kW]
print(f'Heat Capacity {tfc.global_heat_capacity(config, p_outlet=28e5, T_outlet=["Margin", 0])[0]/1e3} [kW]')