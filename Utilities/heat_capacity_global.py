# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 14:17:28 2024

@author: felix
"""










'''

für gedrosselten zustand ist p in anders!!!

'''

from CoolProp.CoolProp import PropsSI as psi
import matplotlib.pyplot as plt
import numpy as np


fuel = 'C2H6O'
p_in = 40e5
T_in = 293



p_array = np.linspace(10e5, 40e5, 100)
T_array = np.linspace(300, 500, 100)

p_array, T_array = np.meshgrid(np.linspace(10e5, 40e5, 500), np.linspace(300, 500, 500))


heat_capacity = np.zeros((len(T_array[:,0]), len(p_array[:,0])))


h_inlet = psi('H', 'P', p_in,'T', T_in, fuel)


for p_index, p in enumerate(p_array[0,:]):
    for T_index, T in enumerate(T_array[:,0]):
        h_outlet = psi('H', 'P', p,'T', T, fuel)
        
        heat_capacity[T_index, p_index] = h_outlet - h_inlet


#%%

fig, ax = plt.subplots()
c = ax.pcolormesh(T_array, p_array/1e5, heat_capacity/1e6/0.8, cmap='Blues', vmax=1)#, vmin = 0, vmax = 700)
# d = ax.pcolormesh(temperature, pressure/1e5, path, cmap='Greys', vmin = 0, vmax = 1)
# e = ax.pcolormesh(temperature, pressure/1e5, speed_of_sound, cmap='Greys', vmin = 0, vmax = 1000)

# c = ax.pcolormesh(temperature, pressure/1e5, specific_heat, cmap='Blues', vmin = 2000, vmax = 1e4)

fig.colorbar(c, ax=ax)
# fig.colorbar(d, ax=ax)
# fig.colorbar(e, ax=ax)


ax.set_xlabel('Temperature [K]')
ax.set_ylabel('Pressure [bar]')
# fig.suptitle(fluid_name)
# plt.hlines(y=63.5, xmin = min(temperature[:,0]), xmax = max(temperature[:,0]), color='r')
# plt.vlines(x=513.9, ymin = min(pressure[0,:]/1e5), ymax = max(pressure[0,:]/1e5), color='r')