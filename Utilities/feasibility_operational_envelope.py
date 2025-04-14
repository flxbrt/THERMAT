# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 19:42:42 2024

@author: felix
"""



import os 
os.chdir(r'C:\Users\felix\Desktop\script_directory')

import therm_functions as tfc

import numpy as np
import cantera as ct

def calc_chocked_mass_flow(d_th, gamma, M, T, p_cc):
        # theta = np.sqrt(gamma*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))))
        theta = np.sqrt(gamma)*((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))
        R_ideal = 8314
        m_dot = np.pi*d_th**2/4*p_cc/np.sqrt(T*R_ideal/M)*theta
        return m_dot

def obtain_combustion_properties(comb):
    def heat_capacity_ratio(comb):
        return comb.cp_mass / comb.cv_mass
    return comb.T, comb.mean_molecular_weight, heat_capacity_ratio(comb)


def calc_effective_exhaust_velocity(p_cc, m_dot, d_th, epsilon, p_ratio, T, M, gamma): # !!! ,p_ex=1e5):
    # p_ratio --> pressure expansion ratio
    
    A_exit = d_th**2/4*np.pi*epsilon
    
    p_ex = p_ratio*p_cc
    
    R_ideal = 8314
    v_ex = np.sqrt((2*gamma*R_ideal)/(gamma-1)*(T)/(M)*(1-(p_ex/p_cc)**((gamma-1)/gamma)))

    p_amb = 1e5
    
    v_eff = v_ex + A_exit/m_dot*(p_ex-p_amb)
    
    return v_eff

p_chamber = np.linspace(9e5, 20e5, 12)
rof = np.linspace(0.7, 1.6, 30)
pinjector_loss = np.linspace(2e5, 8e5, 12)
margin = 20

heat_capacity = np.zeros((len(rof), len(p_chamber)))
heat_capacity_margin = np.zeros((len(rof), len(p_chamber)))
heat_absorbed = np.zeros((len(rof), len(p_chamber)))

isp_eff = np.zeros((len(rof), len(p_chamber)))

#%%


# index_htc_reduction = 0 # htc = heat transfer coefficient


# 0%, 10%, 25%, 50%, 100%


fuel_ct = 'C2H5OH'
fuel_cp = 'C2H6O'
ox = 'O2'

T_fuel = 400
T_ox = 100
T_wall = 700

d_throat = 0.03758523249440576
exp_ratio = 2.92447450715475

for rof_index, rof_val in enumerate(rof):

    for index, p_cc in enumerate(p_chamber):
        
        comb = tfc.set_combustion(p_cc, fuel_ct, fuel_cp, ox, T_fuel, T_ox, rof_val)
        T, M, gamma = obtain_combustion_properties(comb)
        m_dot = calc_chocked_mass_flow(d_throat, gamma, M, T, p_cc)
    
        config = {'ID': f'{p_cc}',
              'chamber': {
                    'd_th': {'Value': d_throat, 'Description': 'Throat Diameter [m]'},
                    'l_star': {'Value': 1, 'Description': 'Characteristic Chamber Length including convergent section [m]'},
                    'alpha': {'Value': 30, 'Description': 'Half angle Convergent section [°]'},
                    'beta': {'Value': 15, 'Description': 'Half angle Divergent section [°]'},
                    'epsilon_c': {'Value': 8, 'Description': 'Contraction Ratio [-]'},
                    'epsilon_e': {'Value': exp_ratio, 'Description': 'Expansion Ratio [-]'}},
              'cooling': {
                    'N': {'Value': 35, 'Description': 'Number of channels [-]'},
                    'h': {'Value': [10e-3], 'Description': 'height of cooling channel [m]'},
                    'tc': {'Value': 1e-3, 'Description': 'wall thickness between cooling channels in circumferential direction along arc [m]'},
                    'tl': {'Value': 1e-3, 'Description': 'thickness of liner between chamber and cooling channel [m]'},
                    'lam': {'Value': 25, 'Description': 'heat conductivity of wall material [W/m/K]'},
                    'p_inlet': {'Value': 40e5, 'Description': 'inlet pressure into cooling channel [Pa]'},
                    'T_wall': {'Value': T_wall, 'Description': ' in case "fix" is chosen as thermal type, this wall temperature is used to calculate the heat flux into the wall [K]'},
                    'T_inlet': {'Value': 293, 'Description': 'inlet temperature into cooling channel [K]'},
                    'fin_efficiency': {'Value': True, 'Description': 'Fin efficiency'}},
              'nodes': {'Value': 100, 'Description': 'Number of axially discretized equidistant elements [-]'},
              'load point': {
                    'p_c': {'Value': p_cc, 'Description': 'total chamber pressure [Pa]'},
                    # 'p_c': {'Value': 9.7e5, 'Description': 'total chamber pressure [Pa]'},
                    'ROF': {'Value': rof_val, 'Description': 'oxidizer to fuel ratio [-]'},
                    'm_dot': {'Value': m_dot, 'Description': 'total mass flow in combustion chamber [kg/s]'},
                    'T_oxidizer': {'Value': T_ox, 'Description': 'Oxidizer Injection Temperature [K]'},
                    'T_fuel': {'Value': T_fuel, 'Description': 'Fuel Injection Temperature [K]'},
                    # 'm_dot': {'Value': 0.544, 'Description': 'total mass flow in combustion chamber [kg/s]'},
                    'oxidizer': {'Value': ox, 'Description': 'chemical formula of oxidizer'},
                    'fuel_cp': {'Value': fuel_cp, 'Description': 'chemical formula of fuel as required by coolprop library'},
                    'fuel_ct': {'Value': fuel_ct, 'Description': 'chemical formula of fuel as required by cantera library'}},
              'settings': {
                    'direction': {'Value': 'counter', 'Description': 'Flow direction of fuel in regenerative cooling channels with respect to hot gas flow direction; either "counter" or "co"'},
                    'thermal type': {'Value': 'fix', 'Description': 'Either "fix" --> fixed wall temperature (needs to be specified) or "regen" --> solves regenerative cooling problem'},
                    'flow type': {'Types': ['frozen'], 'Description': 'Either "frozen" or "equilibrium" or "both" --> affects the chemical composition along the axial direction AND transport properties (cp, lambda, mu etc.)'}
                  }}
    
        geometry = tfc.calc_chamber_geometry(config)
        
        geometry = tfc.discretize_chamber_geometry(config, geometry)
        
        geometry = tfc.calc_discretized_cooling_geometry(config, geometry)
    
        # tfc.plot_geometry(geometry, config)
    
        fluid = tfc.set_thermal_problem(config)
        
        fluid = tfc.calc_flow(config, geometry, fluid)
    
    
        # import time
        # t_start = time.time()
        fluid = tfc.thermal_analysis(config, geometry, fluid)
        # t_end = time.time()
        # print(f'Process took {(t_end-t_start):.2f} [s]')
    
        # tfc.plot_flow(geometry, fluid)
        
        # tfc.plot_constant_wall(geometry, fluid)
        
        # tfc.plot_regen_cooling(geometry, fluid)
        
        heat_absorbed[rof_index, index] = tfc.obtain_heat_absorbed(fluid)
        
        heat_capacity[rof_index, index], _ = tfc.global_heat_capacity(config, p_outlet=p_chamber[index]+pinjector_loss[index], T_outlet=["Margin", 0])
        
        heat_capacity_margin[rof_index, index], _ = tfc.global_heat_capacity(config, p_outlet=p_chamber[index]+pinjector_loss[index], T_outlet=["Margin", margin])
        
        print(f'Outer loop {rof_index+1}')
        print(f'Inner loop {index+1}')



        # Isp abschätzung für jeweiligen rof wert
        p_ratio = 0.0635285993684143
        isp_eff[rof_index, index] = calc_effective_exhaust_velocity(p_cc, m_dot, d_throat, exp_ratio, p_ratio, T, M, gamma)










        
#%%

# import matplotlib.pyplot as plt


# plt.plot(p_chamber/1e5, heat_absorbed[0,:]/1e3, label=f'Heat absorbed [kW] at {T_wall=} [K]')
# plt.plot(p_chamber/1e5, heat_capacity[0,:]/1e3, label='Heat capacity [kW]')
# plt.plot(p_chamber/1e5, heat_capacity_margin[0,:]/1e3, label=f'Heat capacity [kW] at {margin}% margin')
# plt.xlabel('Chamber pressure [bar]')
# plt.ylabel('Heat flow [kW]')
# plt.legend()
# plt.grid()

#%%

# import matplotlib.pyplot as plt
# plt.plot(p_chamber/1e5, heat_capacity[0,:]/1e3, label='Heat capacity [kW]', linestyle = 'dashed')
# plt.plot(p_chamber/1e5, heat_capacity_margin[0,:]/1e3, label=f'Heat capacity [kW] at {margin}% margin', linestyle = 'dashed')

# reduction = [0, 10, 25, 50]

# for ii in range(len(heat_absorbed[:,0])-1):
    

#     plt.plot(p_chamber/1e5, heat_absorbed[ii,:]/1e3, label=f'Heat absorbed [kW] at {T_wall=} [K] with hg reduction of {reduction[ii]}%')

# plt.xlabel('Chamber pressure [bar]')
# plt.ylabel('Heat flow [kW]')
# plt.legend()
# plt.grid()

#%%

# import matplotlib.pyplot as plt

# colors = ['r', 'b', 'g', 'c', 'm', 'y', 'k', 'orange']

# # rof = [0.7, 0.9, 1.1, 1.3]

# for rof_index, rof_val in enumerate(rof):
#     # rof_index = rof_index * 2
#     plt.plot(p_chamber/1e5, heat_capacity[rof_index,:]/1e3, label='Heat capacity [kW]', linestyle = 'dashed', color=colors[rof_index])
#     plt.plot(p_chamber/1e5, heat_capacity_margin[rof_index,:]/1e3, label=f'Heat capacity [kW] at {margin}% margin', linestyle = '-.', color=colors[rof_index])
#     plt.plot(p_chamber/1e5, heat_absorbed[rof_index,:]/1e3, label=f'Heat absorbed [kW] at {T_wall=} [K] with {rof_val=:.1f}', color=colors[rof_index])

# plt.xlabel('Chamber pressure [bar]')
# plt.ylabel('Heat flow [kW]')
# plt.legend()
# plt.grid()


############################


# import matplotlib.pyplot as plt


# for rof_index, rof_val in enumerate(rof):
    
#     plt.plot(p_chamber/1e5, isp_eff[rof_index,:], label=f'Isp effective [m/s] at {rof_val=:.1f}')

# plt.xlabel('Chamber pressure [bar]')
# plt.ylabel('Effective Exhaust Velocity [m/s]')
# plt.legend()
# plt.grid()


############################


feasible = np.zeros(np.shape(heat_absorbed))
feasible_margin = np.zeros(np.shape(heat_absorbed))


for rof_index, rof_val in enumerate(rof):

    for index, p_cc in enumerate(p_chamber):
        
        if heat_capacity[rof_index, index] > heat_absorbed[rof_index, index]:
            feasible[rof_index, index] = 1
        else:
            feasible[rof_index, index] = 0
            
        if heat_capacity_margin[rof_index, index] > heat_absorbed[rof_index, index]:
            feasible_margin[rof_index, index] = 1
        else:
            feasible_margin[rof_index, index] = 0
            
rof_feasible = np.zeros(len(p_chamber))
isp_feasible = np.zeros(len(p_chamber))

for index, p_cc in enumerate(p_chamber):
    for rof_index in range(len(rof)):
        if feasible[rof_index, index] == 1:
            if rof_index == len(rof)-1:
                rof_feasible[index] = rof[-1]
                isp_feasible[index] = isp_eff[-1, index]
        else: # this conditions implies, that the rof is just not feasiible anymore --> requires on rof dicretisation step lower
            if rof_index == 0:
                print(f'chamber pressure of {p_cc/1e5}bar at rof={rof[rof_index]} not feasible')
                rof_feasible[index] = 0
                isp_feasible[index] = 0
            else:
                rof_feasible[index] = rof[rof_index-1]
                isp_feasible[index] = isp_eff[rof_index-1, index]
            break
            
rof_feasible_margin = np.zeros(len(p_chamber))
isp_feasible_margin = np.zeros(len(p_chamber))

for index, p_cc in enumerate(p_chamber):
    for rof_index in range(len(rof)):
        if feasible_margin[rof_index, index] == 1:
            if rof_index == len(rof)-1:
                rof_feasible_margin[index] = rof[-1]
                isp_feasible_margin[index] = isp_eff[-1, index]
        else:
            if rof_index == 0:
                print(f'chamber pressure of {p_cc/1e5}bar at rof={rof[rof_index]} not feasible')
                rof_feasible_margin[index] = 0
                isp_feasible_margin[index] = 0
            else:
                rof_feasible_margin[index] = rof[rof_index-1]
                isp_feasible_margin[index] = isp_eff[rof_index-1, index]
            break

import matplotlib.pyplot as plt

plt.figure(figsize=(12,8))

plt.subplot(121)
plt.plot(p_chamber/1e5, rof_feasible, label='Feasible ROF')
plt.plot(p_chamber/1e5, rof_feasible_margin, label=f'Feasible ROF at {margin}% margin')
plt.xlabel('Chamber pressure [bar]')
plt.ylabel('ROF [-]')
plt.legend()
plt.grid()

plt.subplot(122)
plt.plot(p_chamber/1e5, isp_feasible, label='Feasible Isp eff')
plt.plot(p_chamber/1e5, isp_feasible_margin, label=f'Feasible Isp eff at {margin}% margin')
plt.xlabel('Chamber pressure [bar]')
plt.ylabel('Effective Exhaust Velocity [m/s]')
plt.legend()
plt.grid()
