# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 08:56:53 2024

@author: felix
"""



#%% import libraries

import copy
import numpy as np
import casadi as ca
import cantera as ct
from scipy import optimize
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI as psi

#%% geometry

def calc_chamber_geometry(config):
    d_th = config['chamber']['d_th']['Value']
    l_star = config['chamber']['l_star']['Value']
    alpha = config['chamber']['alpha']['Value']
    beta = config['chamber']['beta']['Value']
    epsilon_c = config['chamber']['epsilon_c']['Value']
    epsilon_e = config['chamber']['epsilon_e']['Value']
    
    if epsilon_c == None:
       epsilon_c = 8*(d_th*100)**(-0.6)+1.25
    A_th = d_th**2/4*np.pi
    V = l_star*A_th # chamber volume including conical section
    A_c = A_th*epsilon_c # cross sectional chamber area
    d_c = np.sqrt(A_c*4/np.pi) # chamber diameter
    A_e = A_th * epsilon_e # cross sectional exit area
    d_e = np.sqrt(A_e*4/np.pi) # exit diameter
    
    l_con = (d_c - d_th)/2/np.tan(alpha/180*np.pi) # length convergent nozzle section
    V_con = l_con*np.pi/3*((d_th/2)**2+d_th*d_c/4+(d_c/2)**2) # volume convergent nozzle section
    
    V_c = V - V_con # volume combustion chamber
    l_c = V_c/(d_c/2)**2/np.pi # length combustion chamber
    
    l_e = (d_e-d_th)/2/np.tan(beta/180*np.pi) # length divergent nozzle section
    
    geometry = {'0D': {
                'epsilon_c': {'Value': epsilon_c, 'Description': 'Contraction Ratio [-]'},
                'd_c': {'Value': d_c, 'Description': 'diameter cylindrical chamber section [m]'},
                'd_e': {'Value': d_e, 'Description': 'diameter at nozzle exit [m]'},
                'd_th': {'Value': d_th, 'Description': 'diameter at throat [m]'},
                'l_con': {'Value': l_con, 'Description': 'length of convergent nozzle section [m]'},
                'l_e': {'Value': l_e, 'Description': 'length of divergent nozzle section [m]'},
                'l_c': {'Value': l_c, 'Description': 'length of cylindrical chamber section [m]'}}}
    
    return geometry 

def discretize_chamber_geometry(config, geometry, axial=None, radial=None):
    
    if geometry != None:
        
        d_th = geometry['0D']['d_th']['Value']
        l_c = geometry['0D']['l_c']['Value']
        d_c = geometry['0D']['d_c']['Value']
        l_con = geometry['0D']['l_con']['Value']
        l_e = geometry['0D']['l_e']['Value']
        d_e = geometry['0D']['d_e']['Value']
        # alpha = config['chamber']['alpha']['Value']
        # beta = config['chamber']['beta']['Value']
        steps = config['nodes']['Value']
        
        # calc radius from diameter  
        r_c = d_c/2
        r_th = d_th/2
        r_e = d_e/2
        
        # obtain total chamber length
        l_tot = l_c + l_con + l_e
        
        # discretize in equidistant steps
        axial = np.linspace(0, l_tot, steps)
        # radius solution array
        radial = np.zeros(steps)
        
        # calculate chamber section lengths according to axial discretisation
        # goal: maintain throat diameter --> throat axial position is therefore shifted downstream or upstream
        def find_nearest(array, value):
            array = np.asarray(array)
            idx = (np.abs(array - value)).argmin()
            return idx, array[idx]
        
        l_c_index, l_c_disc = find_nearest(axial, l_c)
        l_con_index, l_con_disc = find_nearest(axial, l_c+l_con) # l_con_disc --> end of convergent section --> throat should be there
        
        l_con_disc = l_con_disc - l_c_disc
        
        l_e_disc = l_tot - l_c_disc - l_con_disc
        
        # calculate radius for each axial position    
        for l_index, l in enumerate(axial):
    
            if l<=l_c_disc:
                radial[l_index] = r_c
                # l_c_index = l_index
            elif l>l_c_disc and l<=(l_con_disc+l_c_disc):
                radial[l_index] = r_c - (r_c-r_th)/l_con_disc*(l-l_c_disc)
            else:
                radial[l_index] = r_th + (r_e-r_th)/l_e_disc*(l-l_c_disc-l_con_disc)
    
    # calculate wetted combustion chamber surface area at each node
    area_wetted = np.zeros(len(axial)-1)
    for l_index, l in enumerate(axial):         
        
        if l_index == len(axial)-1:
            break
        delta_l = axial[l_index+1] - axial[l_index]
        if geometry != None:
            if l_index < l_c_index: # if elif not required since elif condition always applies --> this also makes following else condition unnecessary
                area_wetted[l_index] = 2*np.pi*radial[l_index]*delta_l
            elif l_index >= l_c_index:
                delta_r = radial[l_index+1]-radial[l_index]
                len_wetted = np.sqrt(delta_r**2 + delta_l**2)
                area_wetted[l_index] = np.pi*(radial[l_index+1] + radial[l_index])*len_wetted
        else:
            delta_r = radial[l_index+1]-radial[l_index]
            len_wetted = np.sqrt(delta_r**2 + delta_l**2)
            area_wetted[l_index] = np.pi*(radial[l_index+1] + radial[l_index])*len_wetted
    
    # calculate chamber cross sectional area at each node
    area_cross = np.zeros(len(axial))
    for r_index, r in enumerate(radial):  
        area_cross[r_index] = r**2*np.pi
        
    if geometry == None:
        geometry = {}
        d_th = config['chamber']['d_th']['Value']
        geometry = {'0D': {'d_th': {'Value': d_th, 'Description': 'diameter at throat [m]'}}}

    geometry['1D'] = {}
    geometry['2D'] = {}
    
    geometry['1D']['axial'] = {'array': axial, 'Description': 'axial position of discretized nodes [m]'}
    geometry['1D']['radial'] = {'array': radial, 'Description': 'radius of chamber according to each axial position [m]'}

    geometry['2D']['area_wetted_hg'] = {'array': area_wetted, 'Description': 'wetted combustion chamber surface area between two nodes [m^2]'}
    geometry['2D']['area_cross_hg'] = {'array': area_cross, 'Description': 'cross sectional chamber in flow direction at each node [m^2]'}

    return geometry

def plot_geometry(geometry, config):
    
    axial = geometry['1D']['axial']['array']
    radial = geometry['1D']['radial']['array']
    area_cross_hg = geometry['2D']['area_cross_hg']['array']
    area_wetted_hg = geometry['2D']['area_wetted_hg']['array']
    
    plt.figure(figsize=(18,8))

    plt.subplot(241)
    plt.plot(axial*1000, radial*1000, color='black')
    plt.plot(axial*1000, -radial*1000, color='black')
    plt.grid()
    plt.xlabel('Chamber X-Coordinate [mm]')
    plt.ylabel('Chamber Radius [mm]')
    plt.axis('equal')
    
    plt.subplot(242)
    plt.plot(axial*1000, area_cross_hg)#, color='black')
    plt.grid()
    plt.xlabel('Chamber X-Coordinate [mm]')
    plt.ylabel('Cross Sectional Area Hot Gas Side [m^2]')
    
    plt.subplot(243)
    plt.plot(axial[:-1]*1000, area_wetted_hg)#, color='black')
    plt.grid()
    plt.xlabel('Chamber X-Coordinate [mm]')
    plt.ylabel('Wetted Area Hot Gas Side [m^2]')
    
    if config['settings']['thermal type']['Value'] == 'regen':
    
        area_wetted_inner_c = geometry['2D']['area_wetted_inner_c']['array']
        area_cross_c = geometry['2D']['area_cross_c']['array']  
        dia_hydraulic = geometry['1D']['dia_hydraulic']['array']
        width = geometry['1D']['width']['array']
        height = geometry['1D']['height']['array']
    
        plt.subplot(244)
        plt.plot(axial[:-1]*1000, area_wetted_inner_c)#, color='black')
        plt.grid()
        plt.xlabel('Chamber X-Coordinate [mm]')
        plt.ylabel('Wetted Area Coolant Side [m^2]')
        
        plt.subplot(245)
        plt.plot(axial*1000, area_cross_c)#, color='black')
        plt.grid()
        plt.xlabel('Chamber X-Coordinate [mm]')
        plt.ylabel('Cross Sectional Area Coolant Side [m^2]')
        
        plt.subplot(246)
        plt.plot(axial*1000, dia_hydraulic)#, color='black')
        plt.grid()
        plt.xlabel('Chamber X-Coordinate [mm]')
        plt.ylabel('Hydraulic diameter [m]')
        
        plt.subplot(247)
        plt.plot(axial*1000, width)#, color='black')
        plt.grid()
        plt.xlabel('Chamber X-Coordinate [mm]')
        plt.ylabel('Channel width [m]')
        
        plt.subplot(248)
        plt.plot(axial*1000, height)#, color='black')
        plt.grid()
        plt.xlabel('Chamber X-Coordinate [mm]')
        plt.ylabel('Channel height [m]')
    
# def plot_geometry(geometry):
    
#     axial = geometry['1D']['axial']['array']
#     radial = geometry['1D']['radial']['array']
    
#     plt.plot(axial*1000, radial*1000, color='black')
#     plt.plot(axial*1000, -radial*1000, color='black')
#     plt.grid()
#     plt.xlabel('Chamber X-Coordinate [mm]')
#     plt.ylabel('Chamber Radius [mm]')
#     plt.axis('equal')

def calc_discretized_cooling_geometry(config, geometry):
    
    axial = geometry['1D']['axial']['array']
    radial = geometry['1D']['radial']['array']
    
    N_channels = config['cooling']['N']['Value']
    h = config['cooling']['h']['Value']
    tl = config['cooling']['tl']['Value'] # thickness between hot gas side and coolant side
    tc = config['cooling']['tc']['Value'] # channel seperator width
    
    if len(h)==1:
        h = np.ones(len(axial))*h[0]
    
    # feasibility check, whether wall thickness between cooling channels is narrow enough
    # only approximation of inner radius, since thickness tl actually measured perpendicular to coolant flow direction; here its in radial direction no matter of the axial position
    r_inner = radial + tl
    circum_inner = 2*r_inner*np.pi
    
    if tc*N_channels >= min(circum_inner):
        raise ValueError('Check cooling channel geometry data. \nReduce channel number N or separator thickness tc.')
    
    # solution array
    # A_inner is wetted surface area of individual cooling channel
    A_inner = np.zeros(len(axial)-1)    
    length = np.zeros(len(axial)-1) 
    
    for index in range(len(A_inner)):  
        deltar = r_inner[index+1]-r_inner[index]
        deltax = axial[index+1]-axial[index]
        leng = np.sqrt(deltax**2+deltar**2)
        A_web = leng*tc
        if deltar == 0:
            # phi = 0
            #############################################################
            ### radial mit r_inner ausgetauscht
            #############################################################
            # A_inner[index] = radial[index]*2*np.pi*leng/N_channels - A_web
            A_inner[index] = r_inner[index]*2*np.pi*leng/N_channels ###################################### - A_web
        else:
            phi = abs(np.arctan(deltar/deltax))
            # v1
            # s1 = radial[index]/np.sin(phi)
            # s2 = radial[index+1]/np.sin(phi)
            # A_inner[index] = abs(s1*radial[index] - s2*radial[index+1])*np.pi/N_channels - A_web
            # v2
            # s1 = r_inner[index]/np.sin(phi)
            # s2 = r_inner[index+1]/np.sin(phi)
            # A_inner[index] = abs(s1*r_inner[index] - s2*r_inner[index+1])*np.pi/N_channels - A_web
            # v3
            #https://studyflix.de/mathematik/oberflaeche-kegel-2612
            r1 = radial[index] + tl/np.cos(phi)
            r2 = radial[index+1] + tl/np.cos(phi)
            s1 = r1/np.sin(phi)
            s2 = r2/np.sin(phi)
            A_inner[index] = abs(s1*r1 - s2*r2)*np.pi/N_channels ###################################### - A_web
        
    A_cross_per_channel = np.zeros(len(axial))    
    
    for l_index, l in enumerate(axial):
        if l_index != len(axial)-1:
            deltar = r_inner[l_index+1]-r_inner[l_index]
            deltax = axial[l_index+1]-axial[l_index]
            leng = np.sqrt(deltax**2+deltar**2)
            if deltar == 0:
                phi = 0
            else:
                phi = np.arctan(deltar/deltax)
            length[l_index] = leng
        # else:
            # pass
        # v1
        # also works for cylindrical part since then it holds that s2 = ro and s1 = ri
        # s1 = (radial[l_index])/np.cos(phi) + tl
        # s2 = s1 + h[l_index]
        # ri = radial[l_index] + tl*np.cos(phi)
        # ro = radial[l_index] + (tl+h[l_index])*np.cos(phi)
        # # https://studyflix.de/mathematik/oberflaeche-kegel-2612
        # A = np.pi*(s2*ro - s1*ri)
        # A_cross_per_channel[l_index] = A/N_channels - h[l_index]*tc
        # v2
        # also works for cylindrical part since then it holds that s2 = ro and s1 = ri
        ri = radial[l_index] + tl/np.cos(phi)
        ro = radial[l_index] + tl/np.cos(phi)+h[l_index]*np.cos(phi)
        s1 = ri/np.cos(phi) # not = (radial[l_index])/np.cos(phi) + tl
        s2 = ro / np.cos(phi) # s1 + h[l_index]
        # https://studyflix.de/mathematik/oberflaeche-kegel-2612
        A = np.pi*(s2*ro - s1*ri)
        A_cross_per_channel[l_index] = A/N_channels - h[l_index]*tc
        
    # calculate hydraulic diameter
    # circum_outer = 2*r_outer*np.pi
    # circum_per_channel = 2*h+(circum_outer+circum_inner)/N_channels - 2*tc
    # dia_hydraulic = 4*A_cross_per_channel/circum_per_channel
    
    circum_per_channel = np.zeros(len(axial))   
    for l_index, l in enumerate(axial):
        if l_index != len(axial)-1:
            deltar = r_inner[l_index+1]-r_inner[l_index]
            deltax = axial[l_index+1]-axial[l_index]
            if deltar == 0:
                phi = 0
            else:
                phi = np.arctan(deltar/deltax)
        else:
            pass
        # v1
        # Ri = (radial[l_index])/np.cos(phi) + tl
        # Ro = Ri + h[l_index]
        # circum_outer = Ro*2*np.pi      
        # circum_inner = Ri*2*np.pi
        # circum_per_channel[l_index] = (circum_outer+circum_inner)/N_channels - 2*tc + 2*h[l_index]
        # v2
        Ri = radial[l_index] + tl/np.cos(phi)
        Ro = Ri + h[l_index]*np.cos(phi)
        circum_outer = Ro*2*np.pi      
        circum_inner = Ri*2*np.pi
        circum_per_channel[l_index] = (circum_outer+circum_inner)/N_channels - 2*tc + 2*h[l_index]

    dia_hydraulic = 4*A_cross_per_channel/circum_per_channel
    
    ###########################################################################
    # calcualte channel width at inner radius of cooling channel --> arc width
    # hopefully bug free version below
    # circum_inner = 2*r_inner*np.pi
    # width = circum_inner/N_channels - tc
    
    circum_inner = np.zeros(len(axial))
    r_inner_equi = np.zeros(len(axial))
    r_hg_equi = np.zeros(len(axial))
    for l_index, l in enumerate(axial):
        if l_index != len(axial)-1:
            deltar = r_inner[l_index+1]-r_inner[l_index]
            deltax = axial[l_index+1]-axial[l_index]
            if deltar == 0:
                phi = 0
            else:
                phi = np.arctan(deltar/deltax)
        r_hg_equi[l_index] = radial[l_index]/np.cos(phi)
        r_inner_equi[l_index] = r_hg_equi[l_index] + tl
        circum_inner[l_index] = (radial[l_index] + tl/np.cos(phi))*2*np.pi
    
    width = circum_inner/N_channels - tc
    
    geometry['1D']['r_hg_equi'] = {'array': r_hg_equi, 'Description': 'equivalent inner radius of combustion chamber [m]'}
    geometry['1D']['r_inner_c_equi'] = {'array': r_inner_equi, 'Description': 'equivalent inner boundary radius of cooling channel [m]'}
    
    ###########################################################################
    
    geometry['2D']['area_wetted_inner_c'] = {'array': A_inner, 'Description': 'Wetted cooling channel area per channel [m^2]'}
    # geometry['2D']['area_outer_c'] = {'array': A_outer, 'Description': '[m^2]'}
    geometry['2D']['area_cross_c'] = {'array': A_cross_per_channel, 'Description': '[m^2]'}
    
    geometry['1D']['dia_hydraulic'] = {'array': dia_hydraulic, 'Description': 'Hydraulic diameter of cooling channel [m]'}
    geometry['1D']['r_inner_c'] = {'array': r_inner, 'Description': 'inner boundary radius of cooling channel [m]'}
    
    geometry['1D']['wetted_length_c'] = {'array': length, 'Description': 'wetted length of cooling channel in axial direction between two adjacent nodes [m]'}
    geometry['1D']['width'] = {'array': width, 'Description': 'channel width at inner radius of channel [m]'}
    geometry['1D']['height'] = {'array': h, 'Description': 'channel height [m]'}

    return geometry

#%% thermal analyis
    
# def set_combustion(p, fuel_ct, fuel_cp, ox, T_fuel, T_ox, ROF):
    
#     comb = ct.Solution('gri30_WARR.yaml')
#     # comb = ct.Solution('gri30_WARR.yaml')
#     comb.Y = {fuel_ct: 1, ox: ROF}
#     comb.TP = 273, p
#     comb.equilibrate('HP')
#     return comb

def set_combustion(p, fuel_ct, fuel_cp, ox, T_fuel, T_ox, ROF):
    def get_delta_h(injec_t, ref_t, p, fluid):
        
        vap_t = psi('T','P',p,'Q',1,fluid)
        
        # ref_t: temperature with which cantera object is initialized
        # vap_t: vaporization temperature of fluid at given pressure
        # injec_t: propellant temperature when injected / tank propellant temperature
        # assumption: propellants are always liquid --> vaporisation temperature is always higher than injection temperature
        # 3 cases need to be distinguished
        # 1: injec_t < vap_t < ref_t
        if ref_t > vap_t:
            h_low = psi('H', 'P', p,'T', injec_t, fluid)
            h_high = psi('H', 'P', p,'T', ref_t, fluid)
            deltah = h_high - h_low
            # print(1)
        # 2: injec_t < ref_t < vap_t
        if ref_t < vap_t and ref_t > injec_t:
            deltah = psi('H', 'P', p,'T', ref_t, fluid) - psi('H', 'P', p,'T', injec_t, fluid) + psi('H','P',p,'Q',1,fluid) - psi('H','P',p,'Q',0,fluid)
            # print(2)
        # 3: ref_t < injec_t < vap_t
        if ref_t < injec_t:
            # negative sign for the first paranthesis, since deltah is subtracted and therefore negative frist term yields less enthalpy being subtracted from cantera object
            deltah =  - (psi('H', 'P', p,'T', injec_t, fluid) - psi('H', 'P', p,'T', ref_t, fluid)) + psi('H','P',p,'Q',1,fluid) - psi('H','P',p,'Q',0,fluid)
            # print(3)
        return deltah
    
    ref_t = 273
    
    comb = ct.Solution('gri30_WARR.yaml')
    # comb = ct.Solution('gri30_WARR.yaml')
    comb.Y = {fuel_ct: 1, ox: ROF}
    comb.TP = 273, p
    comb.equilibrate('HP')
    
    delta_h_ox = get_delta_h(injec_t=T_ox, ref_t=ref_t, p=p, fluid=ox)
    delta_h_fuel = get_delta_h(injec_t=T_fuel, ref_t=ref_t, p=p, fluid=fuel_cp)
    
    deltah = ROF/(1+ROF)*delta_h_ox + 1/(1+ROF)*delta_h_fuel
    
    comb.HP = comb.h-deltah, p
    
    return comb

def set_thermal_problem(config):
    p_c = config['load point']['p_c']['Value']
    fuel_ct = config['load point']['fuel_ct']['Value']
    fuel_cp = config['load point']['fuel_cp']['Value']
    ox = config['load point']['oxidizer']['Value']
    T_fuel = config['load point']['T_fuel']['Value']
    T_ox = config['load point']['T_oxidizer']['Value']

    ROF = config['load point']['ROF']['Value']

    comb = set_combustion(p_c, fuel_ct, fuel_cp, ox, T_fuel, T_ox, ROF)
    comb = ct.Quantity(comb, constant='HP')
    # comb_init = ct.Quantity(comb, constant='HP')
    
    fluid = {'combustion':{
                  'ct object': {'object': comb, 'Description': 'cantera combustion object containing the equilibrium hot gas properties at the injector for an infinite area combustor'}},
            'cooling':{}}
    
    for flow_type in config['settings']['flow type']['Types']:
        if flow_type not in ['frozen', 'equilibrium', 'finite rate']:
            raise NameError(f'Check flow type: {flow_type}')
        fluid['combustion'][flow_type] = {}
    
    # if config['settings']['flow type']['Value'] == 'equilibrium':
    #     fluid['combustion']['equilibrium'] = {}
    # elif config['settings']['flow type']['Value'] == 'frozen':
    #     fluid['combustion']['frozen'] = {}
    # elif config['settings']['flow type']['Value'] == 'both':
    #     fluid['combustion']['equilibrium'] = {}
    #     fluid['combustion']['frozen'] = {}
    # else:
    #     raise NameError(f'Check flow type: {config["settings"]["flow type"]["Value"]}')
    
    d_th = config['chamber']['d_th']['Value']
    m_dot = config['load point']['m_dot']['Value']
    c_star = p_c*d_th**2/4*np.pi/m_dot
    config['load point']['c_star'] = {'Value': c_star, 'Description': 'characteristic velocity [m/s]'}
    
    return fluid

def expand_gas_equilibrium(comb, m_dot, area, section):
    p0 = comb.P
    T0 = comb.T
    s0 = comb.s
    h0 = comb.h
    
    err = 1
    
    dummy = copy.copy(comb)
    
    gamma = dummy.cp / dummy.cv
    
    if section == 'subsonic':
        p_up = p0
        p_low = p0*(2/(gamma+1))**(gamma/(gamma-1))
        while abs(err/area) >= 1e-3:
            p_mid = (p_up + p_low) / 2
            # set the state using (p,s0)
            dummy.SP = s0, p_mid
            dummy.equilibrate('SP')
            
            v = np.sqrt(2.0*(h0 - dummy.h))      # h + V^2/2 = h0
            A_calc = m_dot/(dummy.density*v)    # rho*v*A = constant
            delta_A = A_calc - area
            err = delta_A
            gamma = dummy.cp / dummy.cv
            sonic_vel = np.sqrt(gamma * ct.gas_constant * dummy.T / dummy.mean_molecular_weight)
            Ma = v/sonic_vel
            if delta_A > 0:
                p_up = p_mid
            elif delta_A < 0:
                p_low = p_mid
    elif section == 'throat':
        p_mid = p0*(2/(gamma+1))**(gamma/(gamma-1))
        dummy.SP = s0, p_mid
        sonic_vel = np.sqrt(gamma * ct.gas_constant * dummy.T / dummy.mean_molecular_weight)
        v = sonic_vel
        Ma = 1
    elif section == 'supersonic':
            p_up = p0/1.6
            p_low = 0
            while abs(err/area) >= 1e-6:
                p_mid = (p_up + p_low) / 2
                # set the state using (p,s0)
                dummy.SP = s0, p_mid
                dummy.equilibrate('SP')
                
                v = np.sqrt(2.0*(h0 - dummy.h))      # h + V^2/2 = h0
                A_calc = m_dot/(dummy.density*v)    # rho*v*A = constant
                delta_A = A_calc - area
                err = delta_A
                gamma = dummy.cp / dummy.cv
                sonic_vel = np.sqrt(gamma * ct.gas_constant * dummy.T / dummy.mean_molecular_weight)
                Ma = v/sonic_vel
                if delta_A > 0:
                    p_low = p_mid
                elif delta_A < 0:
                    p_up = p_mid
    
    T = T0*(1+(gamma-1)/2*Ma**2)**(-1)
    p = p_mid # p0*(1+(gamma-1)/2*Ma**2)**(-gamma/(gamma-1))
    
    rho = dummy.density
    
    return p, T, sonic_vel, v, Ma, rho, gamma

def expand_gas_frozen(comb, m_dot, area, section, p_prev = None):
    p0 = comb.P
    # T0 = comb.T
    s0 = comb.s
    h0 = comb.h
    
    err = 1
    
    dummy = copy.copy(comb)
    
    counter = 0
    
    if section == 'subsonic':
        p_up = p0
        p_low = p0/2        
        while abs(err)/area >= 1e-4:
            p_mid = (p_up + p_low) / 2
            # set the state using (p,s0)
            dummy.SP = s0, p_mid
            v = np.sqrt(2.0*(h0 - dummy.h))      # h + V^2/2 = h0
            A_calc = m_dot/(dummy.density*v)    # rho*v*A = constant
            delta_A = A_calc - area
            err = delta_A
            gamma = dummy.cp / dummy.cv
            
            if delta_A > 0:
                p_up = p_mid
            elif delta_A < 0:
                p_low = p_mid
            # print(err, p_up/1e5, p_low/1e5)
            
            counter += 1
            if counter > 10000:
                print('Failed on subsonic node')
                break
            
    elif section == 'throat':
        ## option 1
        # p = p_prev
        # while  abs(err)/area >= 1e-4:
        #     dummy.SP = s0, p
        #     v = np.sqrt(2.0*(h0 - dummy.h))      # h + V^2/2 = h0
        #     A_calc = m_dot/(dummy.density*v)    # rho*v*A = constant
        #     delta_A = A_calc - area
        #     err = delta_A
        #     p -= 0.1e5
            # print(abs(err)/area)
        
        ## option 2            
        p_up = p0/1.6
        p_low = p0/2.0
        while  abs(err)/area >= 1e-2:
            p_mid = (p_up + p_low) / 2
            # set the state using (p,s0)
            dummy.SP = s0, p_mid
            v = np.sqrt(2.0*(h0 - dummy.h))      # h + V^2/2 = h0
            A_calc = m_dot/(dummy.density*v)    # rho*v*A = constant
            delta_A = A_calc - area
            err = delta_A
            gamma = dummy.cp / dummy.cv
            if delta_A > 0:
                # break
                p_up = p_mid
            elif delta_A < 0:
                p_low = p_mid
            # print(err/area)
                
        ## option 3
        # while abs(err) >= 1e-10:
        #     gamma = dummy.cp / dummy.cv
        #     p_mid = p0*(2/(gamma+1))**(gamma/(gamma-1))
        #     dummy.SP = s0, p_mid
        #     sonic_vel = np.sqrt(gamma * ct.gas_constant * dummy.T / dummy.mean_molecular_weight)
        #     err = p_mid - p0*(2/(gamma+1))**(gamma/(gamma-1))    
        # v = sonic_vel

            counter += 1
            if counter > 10000:
                print('Failed at throat node')
                break

    elif section == 'supersonic':
            p_up = p0/1.6
            p_low = 0
            while abs(err)/area >= 1e-4:
                
            # while abs(err/area) >= 1e-3:
                p_mid = (p_up + p_low) / 2
                # set the state using (p,s0)
                dummy.SP = s0, p_mid
                v = np.sqrt(2.0*(h0 - dummy.h))      # h + V^2/2 = h0
                A_calc = m_dot/(dummy.density*v)    # rho*v*A = constant
                delta_A = A_calc - area
                err = delta_A
                gamma = dummy.cp / dummy.cv
                # sonic_vel = np.sqrt(gamma * ct.gas_constant * dummy.T / dummy.mean_molecular_weight)
                # Ma = v/sonic_vel
                if delta_A > 0:
                    p_low = p_mid
                elif delta_A < 0:
                    p_up = p_mid
    
                counter += 1
                if counter > 10000:
                    print('Failed on supersonic node')
                    break
    
    sonic_vel = np.sqrt(gamma * ct.gas_constant * dummy.T / dummy.mean_molecular_weight)
    Ma = v/sonic_vel
    
    cp = dummy.cp
    
    T = dummy.T # T0*(1+(gamma-1)/2*Ma**2)**(-1)
    p = p_mid # p0*(1+(gamma-1)/2*Ma**2)**(-gamma/(gamma-1))
    
    rho = dummy.density
    
    mu = dummy.viscosity 
    lam = dummy.thermal_conductivity
    
    return p, T, sonic_vel, v, Ma, rho, gamma, cp, mu, lam

def calc_flow(config, geometry, fluid):
    
    area_cross = geometry['2D']['area_cross_hg']['array']
    comb = fluid['combustion']['ct object']['object']
    m_dot = config['load point']['m_dot']['Value']
    
    throat_index = np.where(area_cross == min(area_cross))[0][0]
    
    axial = geometry['1D']['axial']['array']
    radial = geometry['1D']['radial']['array']
    
    # if config['settings']['flow type']['Value'] == 'frozen' or config['settings']['flow type']['Value'] == 'both':
    
        # frozen
        # solution arrays
    T_array = np.zeros(len(area_cross))
    p_array = np.zeros(len(area_cross))
    rho_array = np.zeros(len(area_cross))
    a_array = np.zeros(len(area_cross))
    v_array = np.zeros(len(area_cross))
    Ma_array = np.zeros(len(area_cross))
    gamma_array = np.zeros(len(area_cross))
    cp_array = np.zeros(len(area_cross))
    viscosity_array = np.zeros(len(area_cross))
    conductivity_array = np.zeros(len(area_cross))
    Pr_array = np.zeros(len(area_cross))
    Re_array = np.zeros(len(area_cross))

    p_prev = None

    for area_index, area in enumerate(area_cross):
        
        if area_index < throat_index:
            section = 'subsonic'
        elif area_index == throat_index:
            section = 'throat'
            p_prev = p_array[area_index-1]
            # break
        elif area_index > throat_index:
            section = 'supersonic'
            # break
        
        p, T, a, v, Ma, rho, gamma, cp, mu, lam = expand_gas_frozen(comb, m_dot, area, section, p_prev)
        
        T_array[area_index] = T
        p_array[area_index] = p
        a_array[area_index] = a
        v_array[area_index] = v
        Ma_array[area_index] = Ma
        rho_array[area_index] = rho
        gamma_array[area_index] = gamma
        cp_array[area_index] = cp
        viscosity_array[area_index] = mu
        conductivity_array[area_index] = lam
        Pr_array[area_index] = Prandtl(mu, cp, lam)
    
        Re_array[area_index] = rho*v*2*radial[area_index]/mu
    
    # wetted_length = geometry['1D']['wetted_length_c']['array']
    
    delta_l = np.sqrt(np.diff(axial)**2 + np.diff(radial)**2)
    delta_v = np.diff(v_array)    
    Klam_array = mu/rho/v**2*delta_v/delta_l
    
    fluid['combustion']['frozen']['T static'] = {'array': T_array, 'Description': 'Static hot gas temperature [K]'}
    fluid['combustion']['frozen']['p static'] = {'array': p_array, 'Description': 'Static hot gas pressure [p]'}
    fluid['combustion']['frozen']['sonic velocity static'] = {'array': a_array, 'Description': 'Static hot gas sonic velocity [m/s]]'}
    fluid['combustion']['frozen']['velocity'] = {'array': v_array, 'Description': 'velocity [m/s]]'}
    fluid['combustion']['frozen']['Ma static'] = {'array': Ma_array, 'Description': 'Static Mach number [-]'}
    fluid['combustion']['frozen']['rho static'] = {'array': rho_array, 'Description': 'Static hot gas density [kg/m^3]'}
    fluid['combustion']['frozen']['gamma'] = {'array': gamma_array, 'Description': 'Specific heat ratio [-]'}
    fluid['combustion']['frozen']['cp'] = {'array': cp_array, 'Description': 'Specific heat at constant pressure [J/kg/K]'}
    fluid['combustion']['frozen']['mu'] = {'array': viscosity_array, 'Description': 'Dynamic viscosity [Pa s]'}
    fluid['combustion']['frozen']['lam'] = {'array': conductivity_array, 'Description': 'Thermal conductivity [W/m/K]'}
    fluid['combustion']['frozen']['Pr'] = {'array': Pr_array, 'Description': 'Prandtl number [-]'}
    fluid['combustion']['frozen']['K lam'] = {'array': Klam_array, 'Description': 'Relaminarisation Factor [-]'}
    fluid['combustion']['frozen']['Re'] = {'array': Re_array, 'Description': 'Reynolds Number [-]'}

    print("Solved frozen flow")
    
    # if config['settings']['flow type']['Value'] == 'equilibrium' or config['settings']['flow type']['Value'] == 'both':

    #     # equilibrium
    #     # solution arrays
    #     T_array = np.zeros(len(area_cross))
    #     p_array = np.zeros(len(area_cross))
    #     rho_array = np.zeros(len(area_cross))
    #     a_array = np.zeros(len(area_cross))
    #     v_array = np.zeros(len(area_cross))
    #     Ma_array = np.zeros(len(area_cross))
    #     gamma_array = np.zeros(len(area_cross))
        
    #     for area_index, area in enumerate(area_cross):
    #         if area_index < throat_index:
    #             section = 'subsonic'
    #         elif area_index == throat_index:
    #             section = 'throat'
    #         elif area_index > throat_index:
    #             section = 'supersonic'
            
    #         p, T, a, v, Ma, rho, gamma = expand_gas_equilibrium(comb, m_dot, area, section)
            
    #         print(f"Solved {area_index+1} of {len(area_cross)} nodes for equilibrium flow")
            
    #         T_array[area_index] = T
    #         p_array[area_index] = p
    #         a_array[area_index] = a
    #         v_array[area_index] = v
    #         Ma_array[area_index] = Ma
    #         rho_array[area_index] = rho
    #         gamma_array[area_index] = gamma
        
    #     fluid['combustion']['equilibrium']['T static'] = {'array': T_array, 'Description': 'Static hot gas temperature [K]'}
    #     fluid['combustion']['equilibrium']['p static'] = {'array': p_array, 'Description': 'Static hot gas pressure [p]'}
    #     fluid['combustion']['equilibrium']['sonic velocity static'] = {'array': a_array, 'Description': 'Static hot gas sonic velocity [m/s]]'}
    #     fluid['combustion']['equilibrium']['velocity'] = {'array': v_array, 'Description': 'velocity [m/s]]'}
    #     fluid['combustion']['equilibrium']['Ma static'] = {'array': Ma_array, 'Description': 'Static Mach number [-]'}
    #     fluid['combustion']['equilibrium']['rho static'] = {'array': rho_array, 'Description': 'Static hot gas density [kg/m^3]'}
    #     fluid['combustion']['equilibrium']['gamma'] = {'array': gamma_array, 'Description': 'Specific heat ratio [-]'}

    return fluid

def Prandtl(mu, cp, lam):
    Pr = (mu * cp) / lam
    return Pr

def adiabatic_wall(T_array, T0, Pr):
    # Thesis Johanna Kärner p-19 + Diss Kirchberger p- 23
    # TODO: check assumption: turbulent flow
    recovery = Pr**(1/3)
    T_aw = T_array + recovery*(T0-T_array)
    return T_aw

def solve_cooling_constant_wall(config, geometry, fluid):
    d_th = geometry['0D']['d_th']['Value']
    c_star = config['load point']['c_star']['Value']
    T0 = fluid['combustion']['ct object']['object'].T
    # p_chamber = config['load point']['p_c']['Value']
    if 'Value' in config['cooling']['T_wall'].keys():
        T_wall = config['cooling']['T_wall']['Value']
        typ = 'const'
    else:
        T_wall = config['cooling']['T_wall']['Array']
        typ = 'variable'
   
    p_hg_array = fluid['combustion']['frozen']['p static']['array']
    mu_hg_array = fluid['combustion']['frozen']['mu']['array']
    Cp_hg_array = fluid['combustion']['frozen']['cp']['array']
    Pr_hg_array = fluid['combustion']['frozen']['Pr']['array']
    area_cross_hg = geometry['2D']['area_cross_hg']['array']
    gamma_hg_array = fluid['combustion']['frozen']['gamma']['array']
    Ma_array = fluid['combustion']['frozen']['Ma static']['array']
    area_wetted_hg = geometry['2D']['area_wetted_hg']['array']
    # Pr_hg_array = np.flip(fluid['combustion']['frozen']['Pr']['array'])
    T_hg_static = fluid['combustion']['frozen']['T static']['array']
    T_aw_array = adiabatic_wall(T_hg_static, T0, Pr_hg_array)
    
    C = 0.026
    omega = 0.6
    A_star = d_th**2/4*np.pi
    
    # solution arrays
    hg_array = np.zeros(len(area_wetted_hg))
    Q_dot_array = np.zeros(len(area_wetted_hg))
    q_dot_array = np.zeros(len(area_wetted_hg))

    gamma = gamma_hg_array[0]
    mu = mu_hg_array[0]
    cp = Cp_hg_array[0]
    Pr = Pr_hg_array[0]
    
    
    for index in range(len(area_wetted_hg)):
        if typ == 'variable':
            T_w = T_wall[index]
        else:
            T_w = T_wall
            
        p_hg = p_hg_array[0]
        
        sigma = 1 / (((((1 / 2) * (T_w / T0) * (1 + (((gamma - 1) / 2) * Ma_array[index]**2))) + (1 / 2))**(0.8 - (omega / 5))) * ((1 + (((gamma - 1) / 2) * Ma_array[index]**2))**(omega / 5)))
        # assumption r_curvature * d_th --> this term therefore does not appear here
        # hg = (C / (d_th**0.2)) * ((mu**0.2 * cp)/(Pr**0.6)) * ((p_hg_array[index] / c_star)**0.8) * ((A_star / area_cross_hg[index])**0.9) * sigma
        fac = (d_th/(1.5*d_th))**0.1
        hg = (C / (d_th**0.2)) * ((mu**0.2 * cp)/(Pr**0.6)) * ((p_hg / c_star)**0.8) * ((A_star / area_cross_hg[index])**0.9) * sigma # * fac
        
        # hg *= 2/3
        
        hg_array[index] = hg
        
        Q_dot_array[index] = hg*area_wetted_hg[index]*(T_aw_array[index] - T_w)
        q_dot_array[index] = hg*(T_aw_array[index] - T_w)
    
    fluid['combustion']['frozen']['hg'] = {'array': hg_array, 'Description': 'heat transfer coefficient hot gas side [W/m^2/K]'}
    fluid['combustion']['frozen']['T_aw'] = {'array': T_aw_array, 'Description': 'adiabatic wall temperature [K]'}
    fluid['cooling']['Q_dot'] = {'array': Q_dot_array, 'Description': 'Heat flow PERPENDICULAR to chamber wall [W]'}
    fluid['cooling']['q_dot_hg'] = {'array': q_dot_array, 'Description': 'Heat flux with respect to hot gas side wall PERPENDICULAR to chamber wall [W/m^2]'}

    return fluid
  
def critical_heat_flux(T_bulk, p, vel, fuel):
    try:
        T_sat = psi('T', 'P', p,'Q', 1, fuel)
        valid = True
        if p<10e5 or p>48e5:
            valid = False
        if vel < 3 or vel > 23:
            valid = False
        
        vel = vel * 3.28084 # conversion m/s in ft/s
        T_sat = 9/5*(T_sat-273.15)+32 # conversion K in °F
        T_bulk = 9/5*(T_bulk-273.15)+32
        p = p/1e5*14.5038 # conversion from Pa in psi
        
        CHF = 0.1003 + 0.05264*np.sqrt(vel*(T_sat-T_bulk))
        CHF = CHF*1.17-8.56e-4*p
        
        CHF = CHF * 1635346.27# conversion from BTU/in^2/s into W/m^2
    except ValueError:
        CHF = 0
        valid = False
    return CHF, valid
    
  
def thermal_analysis(config, geometry, fluid):
    
    # fixed wall temperature heat flux
    if config['settings']['thermal type']['Value'] == 'fix':
        
        fluid = solve_cooling_constant_wall(config, geometry, fluid)
        
    # regnerative cooling --> wall temperature resolves this
    elif config['settings']['thermal type']['Value'] == 'regen':
        
        d_th = geometry['0D']['d_th']['Value']
        
        c_star = config['load point']['c_star']['Value']
        T0 = fluid['combustion']['ct object']['object'].T
        # p_chamber = config['load point']['p_c']['Value']
        
        T_hg_static = fluid['combustion']['frozen']['T static']['array']
        
        fluid_cp = config['load point']['fuel_cp']['Value']
        m_dot = config['load point']['m_dot']['Value']
        ROF = config['load point']['ROF']['Value']
        m_dot_fuel = m_dot / (1+ROF)
        # m_dot_fuel_channel = m_dot_fuel / config['cooling']['N']['Value']
        N_channels = config['cooling']['N']['Value']
        
        eta_f_flag = config['cooling']['fin_efficiency']['Value']
        
        height = geometry['1D']['height']['array']
        tc = config['cooling']['tc']['Value']
        tl = config['cooling']['tl']['Value']
        
        lam_wall = config['cooling']['lam']['Value']
        
        fuel_cp = config['load point']['fuel_cp']['Value']
        
        ###
        radius_hg = geometry['1D']['radial']['array']
        radius_hg_equi = geometry['1D']['r_hg_equi']['array']
        radius_inner_c = geometry['1D']['r_inner_c']['array']
        wetted_length_c = geometry['1D']['wetted_length_c']['array']
        dia_hydraulic = geometry['1D']['dia_hydraulic']['array']
        width = geometry['1D']['width']['array']
        
        area_cross_hg = geometry['2D']['area_cross_hg']['array']
        area_wetted_hg = geometry['2D']['area_wetted_hg']['array']
        area_wetted_inner_c = geometry['2D']['area_wetted_inner_c']['array']
        area_cross_c = geometry['2D']['area_cross_c']['array']
        
        # frozen flow
        mu_hg = fluid['combustion']['frozen']['mu']['array']
        gamma_hg = fluid['combustion']['frozen']['gamma']['array']
        cp_hg = fluid['combustion']['frozen']['cp']['array']
        Pr_hg = fluid['combustion']['frozen']['Pr']['array']
        Ma_hg = fluid['combustion']['frozen']['Ma static']['array']
        p_hg_static = fluid['combustion']['frozen']['p static']['array']
        
        # T_aw_array = np.flip(adiabatic_wall(T_hg_static, T0, Pr_hg))
        T_aw_array = adiabatic_wall(T_hg_static, T0, Pr_hg)
        
        mu = mu_hg[0]
        cp = cp_hg[0]
        Pr = Pr_hg[0]
        p_hg = p_hg_static[0]
    
        if config['settings']['direction']['Value'] == 'counter':
            height = np.flip(geometry['1D']['height']['array'])
            radius_hg = np.flip(geometry['1D']['radial']['array'])
            radius_hg_equi = np.flip(geometry['1D']['r_hg_equi']['array'])
            radius_inner_c = np.flip(geometry['1D']['r_inner_c']['array'])
            wetted_length_c = np.flip(geometry['1D']['wetted_length_c']['array'])
            dia_hydraulic = np.flip(geometry['1D']['dia_hydraulic']['array'])
            width = np.flip(geometry['1D']['width']['array'])
            
            area_cross_hg = np.flip(geometry['2D']['area_cross_hg']['array'])
            area_wetted_hg = np.flip(geometry['2D']['area_wetted_hg']['array'])
            area_wetted_inner_c = np.flip(geometry['2D']['area_wetted_inner_c']['array'])
            area_cross_c = np.flip(geometry['2D']['area_cross_c']['array'])
            
            # frozen flow
            mu_hg = np.flip(fluid['combustion']['frozen']['mu']['array'])
            gamma_hg = np.flip(fluid['combustion']['frozen']['gamma']['array'])
            cp_hg = np.flip(fluid['combustion']['frozen']['cp']['array'])
            Pr_hg = np.flip(fluid['combustion']['frozen']['Pr']['array'])
            Ma_hg = np.flip(fluid['combustion']['frozen']['Ma static']['array'])
            p_hg_static = np.flip(fluid['combustion']['frozen']['p static']['array'])
            
            T_aw_array = np.flip(adiabatic_wall(T_hg_static, T0, Pr_hg))
        
        # elif config['settings']['direction']['Value'] == 'co':
        #     radius_hg = geometry['1D']['radial']['array']
        #     radius_inner_c = geometry['1D']['r_inner_c']['array']
        #     wetted_length_c = geometry['1D']['wetted_length_c']['array']
        #     dia_hydraulic = geometry['1D']['dia_hydraulic']['array']
        #     width = geometry['1D']['width']['array']
            
        #     area_cross_hg = geometry['2D']['area_cross_hg']['array']
        #     area_wetted_hg = geometry['2D']['area_wetted_hg']['array']
        #     area_wetted_inner_c = geometry['2D']['area_wetted_inner_c']['array']
        #     area_cross_c = geometry['2D']['area_cross_c']['array']
    
        #     # frozen flow
        #     mu_hg = fluid['combustion']['frozen']['mu']['array']
        #     gamma_hg = fluid['combustion']['frozen']['gamma']['array']
        #     cp_hg = fluid['combustion']['frozen']['cp']['array']
        #     Pr_hg = fluid['combustion']['frozen']['Pr']['array']
        #     Ma_hg = fluid['combustion']['frozen']['Ma static']['array']
        #     p_hg_static = fluid['combustion']['frozen']['p static']['array']
            
        #     T_aw_array = adiabatic_wall(T_hg_static, T0, Pr_hg)
        
        # solution arrays
        p_c_array = np.ones(len(radius_hg))*config['cooling']['p_inlet']['Value']
        T_c_array = np.ones(len(radius_hg))*config['cooling']['T_inlet']['Value']
        h_coolant_array = np.ones(len(radius_hg))*psi('H', 'P', p_c_array[0], 'T', T_c_array[0], fluid_cp) # enthalpy
        rho_c_array = np.ones(len(radius_hg))*psi('D', 'P', p_c_array[0], 'T', T_c_array[0], fluid_cp)
        vel_c_array = np.zeros(len(radius_hg)-1)
        
        CHF_array = np.zeros(len(radius_hg)-1)
        CHF_status = np.zeros(len(radius_hg)-1)
        
        hg_array = np.zeros(len(radius_hg)-1)
        hc_array = np.zeros(len(radius_hg)-1)
        Q_dot_array = np.zeros(len(radius_hg)-1)
        T_w_hg_array = np.zeros(len(radius_hg)-1)
        T_w_c_array = np.zeros(len(radius_hg)-1)
        q_dot_hg_array = np.zeros(len(radius_hg)-1)
        eta_f_array = np.zeros(len(radius_hg)-1)
        
        Q_dot_hg_array = np.zeros(len(radius_hg)-1)
        Q_dot_con_array = np.zeros(len(radius_hg)-1)
        Q_dot_c_array = np.zeros(len(radius_hg)-1)
        deltah_array = np.zeros(len(radius_hg)-1)
    
        hg_wo_sigma = 0
        sigma = 0
        
        for index in range(len(area_wetted_hg)):
            print(index+1)
            
            # initialization for casadi solver --> resolves some convergence issues
            init = [hg_array[index-1], hg_wo_sigma, sigma, T_w_hg_array[index-1], 
                    T_w_c_array[index-1], Q_dot_array[index-1]]
            
            h = height[index]
            r_hg = radius_hg[index]
            r_hg_equi = radius_hg_equi[index]
            r_inner_c = radius_inner_c[index]
            node_length = wetted_length_c[index]
            dia_hyd = dia_hydraulic[index]
            w = width[index]
            A_cross_hg = area_cross_hg[index]
            A_hotgas = area_wetted_hg[index]
            A_cooling = area_wetted_inner_c[index]
            A_cross_cooling = area_cross_c[index]
            T_aw = T_aw_array[index]
            T_c = T_c_array[index]
            p_c = p_c_array[index]
            h_coolant = h_coolant_array[index]
            # mu = mu_hg[index]
            gamma = gamma_hg[index]
            # cp = cp_hg[index]
            # Pr = Pr_hg[index]
            Ma = Ma_hg[index]
            # p_cc = p_hg_static[index]
            # print(p_hg/1e5)

            # hg, hg_wo_sigma, sigma, hc, T_w1, T_w2, Q_dot, q_dot_hg, T_c_next, p_c_next, eta_f, h_coolant, Q_dot_hg, Q_dot_conduction, Q_dot_c, rho, vel, deltah = solve_regen_cooling_node(fluid_cp, 
            #             m_dot_fuel, d_th, r_hg, r_inner_c, node_length, dia_hyd, h, tc, w, A_cross_hg, 
            #             A_hotgas, A_cooling, A_cross_cooling, lam_wall, T_aw, p_cc, T_c, 
            #             p_c, h_coolant, c_star, T0, mu, gamma, cp, Pr, Ma, init, eta_f_flag)
            
            hg, hg_wo_sigma, sigma, hc, T_w1, T_w2, Q_dot, q_dot_hg, T_c_next, p_c_next, eta_f, h_coolant, Q_dot_hg, Q_dot_conduction, Q_dot_c, rho, vel, deltah = solve_regen_cooling_node(fluid_cp, 
                        m_dot_fuel, d_th, r_hg, r_hg_equi, r_inner_c, node_length, dia_hyd, h, tc, tl, w, A_cross_hg, 
                        A_hotgas, A_cooling, A_cross_cooling, lam_wall, T_aw, p_hg, T_c, 
                        p_c, h_coolant, c_star, T0, mu, gamma, cp, Pr, Ma, init, N_channels, eta_f_flag)
        
            hg_array[index] = hg
            hc_array[index] = hc
            h_coolant_array[index+1] = h_coolant
            T_w_hg_array[index] = T_w1
            T_w_c_array[index] = T_w2
            Q_dot_array[index] = Q_dot
            T_c_array[index+1] = T_c_next
            p_c_array[index+1] = p_c_next
            q_dot_hg_array[index] = q_dot_hg
            eta_f_array[index] = eta_f
            Q_dot_hg_array[index] = Q_dot_hg
            Q_dot_con_array[index] = Q_dot_conduction
            Q_dot_c_array[index] = Q_dot_c
            rho_c_array[index+1] = rho
            vel_c_array[index] = vel
            deltah_array[index] = deltah
            
            CHF_array[index], CHF_status[index] = critical_heat_flux(T_c_array[index], p_c_array[index], vel, fuel_cp)
            
        if config['settings']['direction']['Value'] == 'counter':
            hg_array = np.flip(hg_array)
            hc_array = np.flip(hc_array)
            h_coolant_array = np.flip(h_coolant_array)
            T_w_hg_array = np.flip(T_w_hg_array)
            T_w_c_array = np.flip(T_w_c_array)
            Q_dot_array = np.flip(Q_dot_array)
            T_c_array = np.flip(T_c_array)
            p_c_array = np.flip(p_c_array)
            q_dot_hg_array = np.flip(q_dot_hg_array)
            T_aw_array = np.flip(T_aw_array)
            eta_f_array = np.flip(eta_f_array)
            Q_dot_hg_array = np.flip(Q_dot_hg_array)
            Q_dot_con_array = np.flip(Q_dot_con_array)
            Q_dot_c_array = np.flip(Q_dot_c_array)
            rho_c_array = np.flip(rho_c_array)
            vel_c_array = np.flip(vel_c_array)
            deltah_array = np.flip(deltah_array)
            CHF_array = np.flip(CHF_array)
        
        fluid['combustion']['frozen']['hg'] = {'array': hg_array, 'Description': 'heat transfer coefficient hot gas side [W/m^2/K]'}
        fluid['combustion']['frozen']['T_aw'] = {'array': T_aw_array, 'Description': 'adiabatic wall temperature [K]'}
        fluid['cooling']['hc'] = {'array': hc_array, 'Description': 'heat transfer coefficient coolant side [W/m^2/K]'}
        fluid['cooling']['T_w_hg'] = {'array': T_w_hg_array, 'Description': 'Wall temperature on hot gas side [K]'}
        fluid['cooling']['T_w_c'] = {'array': T_w_c_array, 'Description': 'Wall temperature on coolant side [K]'}
        fluid['cooling']['Q_dot'] = {'array': Q_dot_array, 'Description': 'Heat flow PERPENDICULAR to chamber wall [W]'}
        fluid['cooling']['T_c'] = {'array': T_c_array, 'Description': 'Coolant temperature  [K]'}
        fluid['cooling']['p_c'] = {'array': p_c_array, 'Description': 'Coolant pressure  [Pa]'}
        fluid['cooling']['q_dot_hg'] = {'array': q_dot_hg_array, 'Description': 'Heat flux with respect to hot gas side wall PERPENDICULAR to chamber wall [W/m^2]'}
        fluid['cooling']['fin efficiency'] = {'array': eta_f_array, 'Description': 'Fin efficiency [-]'}
        fluid['cooling']['Q_dot_hg'] = {'array': Q_dot_hg_array, 'Description': 'Heat flow PERPENDICULAR to chamber wall [W]'}
        fluid['cooling']['Q_dot_conduction'] = {'array': Q_dot_con_array, 'Description': 'Heat flow PERPENDICULAR to chamber wall [W]'}
        fluid['cooling']['Q_dot_c'] = {'array': Q_dot_c_array, 'Description': 'Heat flow PERPENDICULAR to chamber wall [W]'}
        fluid['cooling']['vel_c'] = {'array': vel_c_array, 'Description': 'Velocity in cooling channel [m/s]'}
        fluid['cooling']['rho_c'] = {'array': rho_c_array, 'Description': 'Density in cooling channel [kg/m^3]'}
        fluid['cooling']['deltah'] = {'array': deltah_array, 'Description': 'Enthalpy difference across cooling node [J/kg]'}
        fluid['cooling']['CHF_array'] = {'array': CHF_array, 'Description': 'Critical Heat Flux [W/m^2]'}
        fluid['cooling']['CHF_status'] = {'array': CHF_status, 'Description': 'Critical Heat Flux Correlation Validity [-]'}

    else:
        raise ValueError('Check thermal type')
    
    return fluid

def Reynolds(rho, vel, L, mu):
    Re = rho*vel*L/mu
    return Re

def heat_transfer_cooling(m_dot, A_cross_cooling, dia_hyd, h_c, p_c, N_channels, substance):
    
    m_dot = m_dot/N_channels
    
    rho = psi('D','P',p_c,'H',h_c,substance) # density
    mu = psi('V','P',p_c,'H',h_c,substance) # viscosity
    lam = psi('conductivity','P',p_c,'H',h_c,substance) # thermal conductivity
    cp = psi('C','P',p_c,'H',h_c,substance) # constant pressure specific heat
    vel = m_dot/A_cross_cooling/rho
    Re = Reynolds(rho, vel, dia_hyd, mu)
    Pr = Prandtl(mu, cp, lam)
    # Nu = 0.332*Re**0.5*Pr**0.33 # previous correlation
    # Nu = 0.024*Re**0.8*Pr**0.37
    # Nu = 0.023*Re**0.8*Pr**0.3 # suggested by chat gpt --> dittus boelter with 0.3 for cooling (0.4 for heating)
    Nu = 0.023*Re**0.8*Pr**0.4 # RPA
    # Nu = 0.033*Re**0.8*Pr**0.4 # RPA for hydrogen
    
    # temporary for methane
    # Nu = 0.0185*Re**0.8*Pr**0.4
    
    hc = Nu*lam/dia_hyd
    
    return hc, lam

def fin_efficiency(hc, h, tc, w, lam, eta_f_flag):
    # Christoph Kirchberger
    # term = np.sqrt(2*hc*tc/lam) 
    # eta_f = np.tanh(h/tc*term)/term  
    # if eta_f_flag:
    #     hc = hc*(w+eta_f*2*h)/(w+tc)
    # return hc, eta_f
    # according to RPA --> significant idfference to kirchberger
    term = np.sqrt(2*hc*tc/lam)*h/tc 
    eta_f = w/(w+tc)+2*h/(w+tc)*np.tanh(term)/term
    if eta_f_flag:
        hc = hc*eta_f
    return hc, eta_f

def pressure_loss(m_dot, A_cross_cooling, dia_hyd, node_length, h_c, p_c, N_channels, substance):
    def colebrook_white(f, epsilon, dia_hyd, Re):
        # https://en.wikipedia.org/wiki/Darcy%E2%80%93Weisbach_equation
        # https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.newton.html
        # print(args)
        # epsilon = args[0]
        # dia_hyd = args[1]
        # Re = args[2]
        return 1/np.sqrt(f) + 2*np.log10(epsilon/3.7/dia_hyd + 2.51/Re/np.sqrt(f))
    
    m_dot = m_dot/N_channels
    
    rho = psi('D','P',p_c,'H',h_c,substance) # density
    mu = psi('V','P',p_c,'H',h_c,substance) # viscosity
    vel = m_dot/A_cross_cooling/rho
    Re = Reynolds(rho, vel, dia_hyd, mu)
    
    # TODO: check for realistic epsilon value
    epsilon = 0
    # args = (epsilon, dia_hyd, Re_coolant)
    
    f = optimize.newton(colebrook_white, 1e-3, args=(epsilon, dia_hyd, Re,))
    delta_p = f*node_length*rho/2*vel**2/dia_hyd
    return delta_p, rho, vel

def solve_regen_cooling_node(substance, m_dot_coolant, d_th, r_hg, r_hg_equi, r_inner_c, node_length, dia_hyd, h, tc, tl, w, A_cross_hg, A_hotgas, A_cooling, A_cross_cooling, lam_wall, T_aw, p_cc, T_c, p_c, h_coolant, c_star, T0, mu_hg, gamma_hg, Cp_hg, Pr_hg, Ma, init, N_channels, eta_f_flag):
    
    hc, lam_coolant = heat_transfer_cooling(m_dot_coolant, A_cross_cooling, dia_hyd, h_coolant, p_c, N_channels, substance)
    # hc, eta_f = fin_efficiency(hc, h, tc, w, lam_coolant, eta_f_flag)
    
    # print(hc)
    hc, eta_f = fin_efficiency(hc, h, tc, w, lam_wall, eta_f_flag)
    # print(eta_f)
    
    # hc = 11000
    
    # hc *= 2
    
    # 2) solve implicit problem with casadi
    # constants
    C = 0.026
    omega = 0.6
    A_star = d_th**2/4*np.pi
    
    opti = ca.Opti()
    s = opti.variable(6)

    hg = s[0]
    hg_wo_sigma = s[1]  
    sigma = s[2]
    # hc = s[]
    T_w1 = s[3]             
    T_w2 = s[4]
    Q_dot = s[5]  
    
    opti.subject_to(hg - hg_wo_sigma*sigma == 0)
    opti.subject_to(hg_wo_sigma - (C / (d_th**0.2)) * ((mu_hg**0.2 * Cp_hg)/(Pr_hg**0.6)) * ((p_cc / c_star)**0.8) * ((A_star / A_cross_hg)**0.9) == 0)
    opti.subject_to(sigma - (1 / (((((1 / 2) * (T_w1 / T0) * (1 + (((gamma_hg - 1) / 2) * Ma**2))) + (1 / 2))**(0.8 - (omega / 5))) * ((1 + (((gamma_hg - 1) / 2) * Ma**2))**(omega / 5)))) == 0)
    
    opti.subject_to(Q_dot - A_hotgas*hg*(T_aw-T_w1) == 0)
    opti.subject_to(Q_dot - A_cooling*N_channels*hc*(T_w2-T_c) == 0)   
    
    # v1
    # opti.subject_to(Q_dot - (T_w1-T_w2)/(np.log(r_inner_c/r_hg)/(lam_wall*2*ca.pi*node_length)) == 0)
    # 1D cartesian coordinates as a test
    # opti.subject_to(Q_dot - (T_w1-T_w2)/((r_inner_c-r_hg)/(lam_wall*2*ca.pi*node_length*r_hg)) == 0)
    # v2
    opti.subject_to(Q_dot - (T_w1-T_w2)/(np.log((r_hg_equi+tl)/r_hg_equi)/(lam_wall*2*ca.pi*node_length)) == 0)

    opts = {'ipopt.print_level':0, 'print_time':0}
    opti.solver('ipopt', opts)

    opti.set_initial(s, init)

    sol = opti.solve()
    params = opti.value(s)
    hg = params[0]
    hg_wo_sigma = params[1]
    sigma = params[2]
    # hc = params[3]
    T_w1 = params[3]
    T_w2 = params[4]
    Q_dot = params[5]
    
    q_dot_hg = hg*(T_aw-T_w1)
    
    # 3) obtain pressure losses
    # TODO: rework delta p correlation
    delta_p, rho, vel = pressure_loss(m_dot_coolant, A_cross_cooling, dia_hyd, node_length, h_coolant, p_c, N_channels, substance)
    # print(delta_p/1e5)
    # delta_p = 0
    p_c_next = p_c - delta_p
    
    # 4) calcualte new fluid temperature
    # T_in = T_c
    # h1 = psi('H','P',p_c,'T',T_in,substance)
    h1 = h_coolant
    h2 = h1 + Q_dot/m_dot_coolant
    T_out = psi('T','P',p_c_next,'H',h2,substance)
    T_c_next = T_out
    
    # senity check
    Q_dot_hg = A_hotgas*hg*(T_aw-T_w1)
    Q_dot_conduction = (T_w1-T_w2)/(ca.log(r_inner_c/r_hg)/(lam_wall*2*ca.pi*node_length))
    Q_dot_c = A_cooling*N_channels*hc*(T_w2-T_c)
    deltah = h2 - h1
    
    return hg, hg_wo_sigma, sigma, hc, T_w1, T_w2, Q_dot, q_dot_hg, T_c_next, p_c_next, eta_f, h2, Q_dot_hg, Q_dot_conduction, Q_dot_c, rho, vel, deltah

def plot_flow(geometry, fluid):
    axial = geometry['1D']['axial']['array']
    
    p_static = fluid['combustion']['frozen']['p static']['array']
    T_static = fluid['combustion']['frozen']['T static']['array']
    T_aw = fluid['combustion']['frozen']['T_aw']['array']
    vel = fluid['combustion']['frozen']['velocity']['array']
    a_static = fluid['combustion']['frozen']['sonic velocity static']['array']
    Ma_static = fluid['combustion']['frozen']['Ma static']['array']
    Ma_static_calc = vel/a_static
    gamma = fluid['combustion']['frozen']['gamma']['array']
    cp = fluid['combustion']['frozen']['cp']['array']
    lam = fluid['combustion']['frozen']['lam']['array']
    mu = fluid['combustion']['frozen']['mu']['array']
    Pr = fluid['combustion']['frozen']['Pr']['array']
    rho = fluid['combustion']['frozen']['rho static']['array']
    Re = fluid['combustion']['frozen']['Re']['array']
    
    plt.figure(figsize=(18,8))

    plt.subplot(231)
    plt.plot(axial*1000, p_static/1e5)#, label='Hot Gas side')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Static Pressure [bar]')
    plt.grid()
    
    plt.subplot(232)
    plt.plot(axial*1000, T_static, label='Static')
    plt.plot(axial*1000, T_aw, label='Adiabatic Wall')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Temperature [K]')
    plt.legend()
    plt.grid()
    
    plt.subplot(233)
    plt.plot(axial*1000, vel, label='Flow Velocity')
    plt.plot(axial*1000, a_static, label='Static Speed of Sound')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Velocity [m/s]')
    plt.legend()
    plt.grid()
    
    plt.subplot(234)
    plt.plot(axial*1000, Ma_static, label='Tool')
    plt.plot(axial*1000, Ma_static_calc, label='Manually Calculated')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Mach Number [-]')
    plt.legend()
    plt.grid()
    
    plt.subplot(235)
    plt.plot(axial*1000, rho)#, label='Flow Velocity')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Density [kg/m^3]')
    plt.grid()
    
    plt.plot()
    
    plt.figure(figsize=(18,8))
    
    plt.subplot(231)
    plt.plot(axial*1000, cp)#, label='Hot Gas side')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Specific Heat [J/kg/K]')
    plt.title('Cp')
    plt.grid()
    
    plt.subplot(232)
    plt.plot(axial*1000, lam)#, label='Hot Gas side')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Heat Conductivity [W/m/K]')
    plt.title('Lambda')
    plt.grid()
    
    plt.subplot(233)
    plt.plot(axial*1000, mu)#, label='Hot Gas side')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Viscosity [Pa s]')
    plt.title('Mu')
    plt.grid()
    
    plt.subplot(234)
    plt.plot(axial*1000, gamma)#, label='Hot Gas side')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Heat Capacity Ratio [-]')
    plt.title('Gamma')
    plt.grid()
    
    plt.subplot(235)
    plt.plot(axial*1000, Pr)#, label='Hot Gas side')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Prandtl Number [-]')
    plt.title('Pr')
    plt.grid()
    
    plt.subplot(236)
    plt.plot(axial*1000, Re)#, label='Hot Gas side')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Reynolds Number [-]')
    plt.title('Pr')
    plt.grid()
    
    plt.plot()

def plot_regen_cooling(geometry, fluid):
    
    axial = geometry['1D']['axial']['array']
    
    hg = fluid['combustion']['frozen']['hg']['array']
    T_aw = fluid['combustion']['frozen']['T_aw']['array']
    
    T_c = fluid['cooling']['T_c']['array']
    hc = fluid['cooling']['hc']['array']
    p_c = fluid['cooling']['p_c']['array']
    q_dot_hg = fluid['cooling']['q_dot_hg']['array']
    T_w_c = fluid['cooling']['T_w_c']['array']
    T_w_hg = fluid['cooling']['T_w_hg']['array']
    Q_dot = fluid['cooling']['Q_dot']['array']
    CHF = fluid['cooling']['CHF_array']['array']
    
    # Q_dot_hg = fluid['cooling']['Q_dot_hg']['array']
    # Q_dot_con = fluid['cooling']['Q_dot_conduction']['array']
    # Q_dot_c = fluid['cooling']['Q_dot_c']['array']
    
    # eta_f_array = fluid['cooling']['fin efficiency']['array']
    
    vel_c_array = fluid['cooling']['vel_c']['array']
    rho_c_array = fluid['cooling']['rho_c']['array']

    
    plt.figure(figsize=(24,8))

    plt.subplot(251)
    plt.plot(axial[:-1]*1000, T_w_hg, label='Hot Gas side')
    plt.plot(axial[:-1]*1000, T_w_c, label='Coolant side')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Wall Temperature [K]')
    plt.legend()
    plt.grid()
    
    plt.subplot(252)
    plt.plot(axial*1000, T_c, label='Coolant')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Temperature [K]')
    plt.legend()
    plt.grid()
    
    plt.subplot(253)
    plt.plot(axial*1000, p_c/1e5, label='Coolant')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Pressure [bar]')
    plt.legend()
    plt.grid()
    
    plt.subplot(254)
    # plt.plot(axial[:-1]*1000, hg, label='Hot gas side')
    plt.plot(axial[:-1]*1000, hc, label='Coolant side')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Heat transfer coefficient [W/m^2/K]')
    # plt.legend()
    plt.grid()
    
    plt.subplot(255)
    plt.plot(axial[:-1]*1000, hg, label='Hot gas side')
    # plt.plot(axial[:-1]*1000, hc, label='Coolant side')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Heat transfer coefficient [W/m^2/K]')
    # plt.legend()
    plt.grid()
    
    # plt.subplot(245)
    # plt.plot(axial*1000, T_aw, label='Adiabatic Wall')
    # # plt.plot(timesteps, pref[0,:-1]/1e5)
    # plt.xlabel('Distance from faceplate [mm]')
    # plt.ylabel('Temperature [K]')
    # plt.legend()
    # plt.grid()
    
    plt.subplot(256)
    plt.plot(axial[:-1]*1000, q_dot_hg/1e6, label='Predicted Heat Flux')
    plt.plot(axial[:-1]*1000, CHF/1e6, label='CHF')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Heat flux [MW/m^2]')
    plt.legend()
    plt.grid()
    
    plt.subplot(257)
    plt.plot(axial[:-1]*1000, vel_c_array)
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Velocity [m/s]')
    plt.grid()
    
    plt.subplot(258)
    plt.plot(axial[:]*1000, rho_c_array, label ='hot gas')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Density [kg/m^3]')
    plt.grid()
    plt.legend()
    
    plt.subplot(259)
    plt.plot(axial[:-1]*1000, np.cumsum(Q_dot)/1e3)
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Cummulative Heat Flow until respective axial distance [kW]')
    plt.grid()
    
    # plt.subplot(247)
    # plt.plot(axial[:-1]*1000, eta_f_array)
    # # plt.plot(timesteps, pref[0,:-1]/1e5)
    # plt.xlabel('Distance from faceplate [mm]')
    # plt.ylabel('Fin efficiency [-]')
    # plt.grid()
    
    # plt.subplot(248)
    # plt.plot(axial[:-1]*1000, Q_dot_hg, label ='hot gas')
    # plt.plot(axial[:-1]*1000, Q_dot_con, label ='conduction')
    # plt.plot(axial[:-1]*1000, Q_dot_c, label ='coolant')
    # # plt.plot(timesteps, pref[0,:-1]/1e5)
    # plt.xlabel('Distance from faceplate [mm]')
    # plt.ylabel('Heat flow [W]')
    # plt.grid()
    # plt.legend()


def obtain_heat_of_combustion(config):
    # source = https://www.energie-lexikon.info/ethanol.html
    heizwert = 26.8 # [MJ/kg]
    
    ROF = config['load point']['ROF']['Value']
    m_dot = config['load point']['m_dot']['Value']
    m_dot_ethanol = 1/(1+ROF)*m_dot
    
    return heizwert*m_dot_ethanol*1000 # [kW]
    
    
# theoretically the following function should result in the correct energy release
# however since cantera uses HP to come up with the chemical equilibrium, the enthalpies before and after the reaction are equivalent --> results is thereofre 0
# def obtain_heat_of_combustion(config, fluid, water='gas'):
#     # does not account for liquid injection
#     # ROF_stoichiometric = 2.09
#     fuel_ct = config['load point']['fuel_ct']['Value']
#     ROF = config['load point']['ROF']['Value']
#     comb_before_reaction = ct.Solution('gri30_WARR.yaml')
#     comb_before_reaction.Y = {fuel_ct: 1, 'O2': ROF}
#     reactants = comb_before_reaction.species_names
#     mole_fractions_reactants = comb_before_reaction.X
    
#     comb = fluid['combustion']['ct object']['object']
#     products = comb.species_names
#     mole_fractions_products = comb.X
    
#     t_ref = 298.15
#     deltaH_reaction = 0
#     for species in reactants:
#         index = reactants.index(species)
#         # outputs kJ/mole when considering /1e6
#         # Si unit is J/kmole
#         deltaH_reaction -= mole_fractions_reactants[index]*comb.species(species).thermo.h(t_ref)/1e6
    
#     for species in products:
#         index = reactants.index(species)
#         # outputs kJ/mole when considering /1e6
#         # Si unit is J/kmole
#         deltaH_reaction += mole_fractions_reactants[index]*comb.species(species).thermo.h(t_ref)/1e6
    
#     # convert into kJ/kg
#     molar_mass_ethanol = 46
#     deltaH_reaction = deltaH_reaction*molar_mass_ethanol*1000
    
#     m_dot = config['load point']['m_dot']['Value']
#     m_dot_ethanol = 1/(1+ROF)*m_dot
    
#     # in kJ
#     heat_released = m_dot_ethanol*deltaH_reaction
    
#     return heat_released    
    

def obtain_heat_absorbed(fluid):
    Q_dot = sum(fluid['cooling']['Q_dot']['array'])
    return Q_dot

def global_heat_capacity(config, p_outlet, T_outlet=[None, None]):
    fuel = config['load point']['fuel_cp']['Value']
    
    m_dot = config['load point']['m_dot']['Value']
    ROF = config['load point']['ROF']['Value']
    m_dot_fuel = 1/(1+ROF)*m_dot
    
    # channel inlet conditions
    T_inlet = config['cooling']['T_inlet']['Value']
    p_inlet = config['cooling']['p_inlet']['Value']
    
    # channel outlet conditions
    if T_outlet[0] == None:
        T_outlet = psi('T', 'P', p_outlet,'Q', 1, fuel) - 1
        
    elif T_outlet[0].lower() == 'margin':
        margin = T_outlet[1]/100
        T_outlet = psi('T', 'P', p_outlet,'Q', 1, fuel) - 1
        h_outlet = psi('H', 'P', p_outlet,'T', T_outlet, fuel)
        h_inlet = psi('H', 'P', p_inlet,'T', T_inlet, fuel)
        
        deltah = h_outlet - h_inlet
        h_outlet_margin = h_inlet + deltah*(1-margin)
        T_outlet = psi('T', 'P', p_inlet,'H', h_outlet_margin, fuel)
        
    elif T_outlet[0].lower() == 'value':
        T_outlet = T_outlet[1]
    
    h_outlet = psi('H', 'P', p_outlet,'T', T_outlet, fuel)
    h_inlet = psi('H', 'P', p_inlet,'T', T_inlet, fuel)
    
    deltah = h_outlet - h_inlet
    
    Q = deltah*m_dot_fuel
    
    return Q, T_outlet
    

def plot_constant_wall(geometry, fluid):
    
    axial = geometry['1D']['axial']['array']
    
    hg = fluid['combustion']['frozen']['hg']['array']
    T_aw = fluid['combustion']['frozen']['T_aw']['array']
    
    q_dot_hg = fluid['cooling']['q_dot_hg']['array']
    Q_dot = fluid['cooling']['Q_dot']['array']
    
    Klam = fluid['combustion']['frozen']['K lam']['array']
    
    plt.figure(figsize=(12,8))

    plt.subplot(231)
    plt.plot(axial[:-1]*1000, hg, label='Hot gas side')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Heat transfer coefficient [W/m^2/K]')
    plt.legend()
    plt.grid()
    
    plt.subplot(232)
    plt.plot(axial*1000, T_aw, label='Adiabatic Wall')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Temperature [K]')
    plt.legend()
    plt.grid()
    
    plt.subplot(233)
    plt.plot(axial[:-1]*1000, q_dot_hg/1e6)
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Heat flux [MW/m^2]')
    plt.grid()
    
    plt.subplot(234)
    plt.plot(axial[:-1]*1000, Q_dot)
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Heat flow per element [W]')
    plt.grid()
    
    plt.subplot(235)
    plt.plot(axial[:-1]*1000, np.cumsum(Q_dot)/1e3)
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Cummulative Heat Flow until respective axial distance [kW]')
    plt.grid()
    
    plt.subplot(236)
    plt.plot(axial[:-1]*1000, Klam)
    plt.hlines(y=3e-6, xmin=axial[0]*1000, xmax=axial[-1]*1000, linestyle = 'dashed', color = 'r')
    # plt.plot(timesteps, pref[0,:-1]/1e5)
    plt.xlabel('Distance from faceplate [mm]')
    plt.ylabel('Relaminsarisation Factor [-]')
    plt.grid()
    
    
