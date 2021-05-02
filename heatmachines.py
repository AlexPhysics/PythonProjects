#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 12:24:32 2020

@author: julienreveillon
"""
import cantera as ct
import sys

def pump_comp(fluid, p_o, eta):
    """pump or compress a fluid from p_i to p_o, using
    a pump or compressor with isentropic efficiency eta."""
    # initialization
    h_i         = fluid.h 
    s_i         = fluid.s     
    # isentropic enthalpy
    s_o_is      = s_i 
    fluid.SP    = s_o_is, p_o  #.SP : change of entropy and pressuire of fluid
    h_o_is      = fluid.h   
    # efficiency
    w_is        = h_o_is - h_i   # isentropic work
    w_actual    = w_is / eta     # actual work
    # actual fluid   
    h_o         = h_i + w_actual
    fluid.HP    = h_o, p_o #.HP : change of enthalpy and pressure of fluid
    return fluid

def turbine(fluid, p_o, eta):
    """expand a fluid from p_i to p_o, using
    a turbine with isentropic efficiency eta."""
    # initialization
    h_i         = fluid.h 
    s_i         = fluid.s     
    # isentropic enthalpy
    s_o_is      = s_i 
    fluid.SP    = s_o_is, p_o  #.SP : change of entropy and pressuire of fluid
    h_o_is      = fluid.h   
    # efficiency
    w_is        = h_o_is - h_i   # isentropic work
    w_actual    = w_is * eta     # actual work
    # actual fluid   
    h_o         = h_i + w_actual
    fluid.HP    = h_o, p_o #.HP : change of enthalpy and pressure of fluid
    return fluid    


def isclose(a,b,abs_tol=1e-4):
    """return True if a and b are close within absolute tolerance abs_tol"""
    #
    diff = abs(b - a)
    return(diff<=abs_tol)

def checkFirstLaw(check):
    """ check is first law is respected"""
    if(not isclose(check,0)):
        print('-----------------------------------')
        print('Fatal First Law Error')
        print('sum q + sum w  = ',check)
        print('-----------------------------------')   
        sys.exit()

def F2K(T_in_F):
    """ returns the F temperature in  K """
    return((T_in_F - 32) * 5/9 + 273.15)

def C2K(T_in_C):
    """ returns the C temperature in  K """
    return(T_in_C + 273.15)

def K2C(T_in_K):
    """ returns the K temperature in  C """
    return(T_in_K - 273.15)

def MPa2Pa(P_in_MPa):
    """ returns the megaPa pressure in Pa """
    return (P_in_MPa*1000000.0)

def kPa2Pa(P_in_MPa):
    """ returns the megaPa pressure in Pa """
    return (P_in_MPa*1000.0)

def bar2Pa(P_in_bar):
    """ returns the bar pressure in Pa """
    return (P_in_bar*100000.0)

def psi2Pa(P_in_psi):
    """ returns the psi pressure in Pa """
    return (P_in_psi*6894.7572931783)

def atm2Pa(P_in_atm):
    """ returns the atm pressure in Pa """
    return (P_in_atm*ct.one_atm)

def inch2m(L_in_inch):
    """return the inch lengh in m"""
    return(L_in_inch*0.0254)

def mm2m(L_in_mm):
    """return the mm lengh in m"""
    return(L_in_mm*0.001)

def MJ2J(E_in_MJ):
    """return the megaJ energy in J"""
    return(E_in_MJ*1000000)

def MW2W(P_in_MW):
    """return the megaW power in W"""
    return(P_in_MW*1000000)

def W2HP(P_in_W):
    """return the W power in HorsePower"""
    return(P_in_W*0.0013596216)

def min2sec(t_in_min):
    """return the min in sec"""
    return(t_in_min*60)    


