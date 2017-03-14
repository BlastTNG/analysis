#! /usr/bin/python

import numpy as np
import pylab as pl
from scipy.constants import c, h, k

def from_K_antenna_to_JySR(temp, nu):
    return 2*k*(nu/c)**2**temp/1.0e-26

def from_Jy_to_K_antenna(power, nu):
    return power / 1.0e26 / (2*k*(nu/c)**2)

def from_K_antenna_to_K_thermo(k_antenna, nu):
    t_cmb = 2.725 #K
    x = h*nu / (k*t_cmb)
    return k_antenna * (np.exp(x) - 1)**2 / (x**2 * np.exp(x))

def from_K_thermo_to_K_antenna(k_thermo, nu):
    t_cmb = 2.725 #K
    x = h*nu / (k*t_cmb)
    return k_thermo * (x**2 * np.exp(x)) / (np.exp(x) - 1)**2 

def from_JySR_to_K_thermo(power, nu):
    t_antenna = from_JySR_to_K_antenna(power, nu)
    t_thermo = from_K_antenna_to_K_thermo(t_antenna, nu)
    return t_thermo

def from_Jy_to_K_thermo_deriv_around_exact(power_in_JySR, nu):
    power_in_normal_units = power_in_JySR * 1.0e-26
    Tcmb = 2.725
    x = h*nu / (k*Tcmb)
    dI_over_dT = 2*h**2*nu**4 / (c**2*k*Tcmb**2) * np.exp(x) / (np.exp(x)-1)**2 
    return  power_in_normal_units / dI_over_dT

def from_JySR_to_K_thermo_exact(power_in_JySR, nu):
    power_in_normal_units = power_in_JySR * 1.0e-26
    one_over_T = np.log(2*h*nu**3/c**2 / power_in_normal_units + 1) * k / (h*nu)
    T = 1.0/one_over_T
    return T

def from_K_thermo_to_JySR_exact(T_thermo, nu):
    power_in_normal_units = 2*h*nu**3/c**2 / (np.exp(h*nu/(T_thermo * k)) - 1)
    power_in_JySR = power_in_normal_units * 1.0e26
    return power_in_JySR

def from_K_thermo_to_JySR(t_thermo, nu):
    t_antenna = from_K_thermo_to_K_antenna(t_thermo, nu)
    power = from_K_antenna_to_JySR(t_antenna, nu)
    return power

def black_body_Jy(nu, T_thermo):
    I_Jy = 2*h*nu**3 / c**2 * 1.0/(np.exp(h*nu/(k*T_thermo)) - 1) * 1.0e26
    return I_Jy

def scale_dust_RJ(nu_ref, nu, input_T_RJ):
    Td = 19.6
    beta = 1.57 # Planck 2016 polarization scaling
    output_T_RJ = (nu/nu_ref)**(beta-2) * black_body_Jy(nu, Td) / black_body_Jy(nu_ref, Td) * input_T_RJ
    return output_T_RJ




if __name__ == "__main__":
    #nu = 22.5e9
    power_Jy = 2.23e-1*1.0e6
    t_antenna_K = 1.0
    new_nu = 150e9
    nu = 600e9
    input_uK_RJ = 7.13
    output_uK_RJ = scale_dust_RJ(nu, new_nu, input_uK_RJ)
    print "input ", input_uK_RJ, " uK_RJ at 600 GHz, output ", output_uK_RJ, " uK_RJ at 150 GHz which corresponds to ", from_K_antenna_to_K_thermo(output_uK_RJ, new_nu), "uK CMB"
